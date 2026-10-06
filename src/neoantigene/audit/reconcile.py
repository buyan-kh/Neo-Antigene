"""Can a published epitope actually be derived from its published variant?

This is the audit's cheapest oracle and its sharpest one. A paper claims a
peptide, a gene, and a protein change. Either the reference proteome yields
that peptide from that change or it does not, and the answer is exact
arithmetic rather than judgment. No assay, no model, no opinion.

It matters because a claimed epitope that does not reconcile invalidates a row
in somebody's benchmark, and because the *real* epitope in that case is
something else — which is a correction with biological content rather than
bookkeeping. This repo already hit one instance in its own labels: CASP1
p.P172S matched no Ensembl 116 isoform, because no isoform in that release
carries proline at 172.

Three subtleties are why this belongs in the library rather than being left to
whoever is auditing:

**Isoform choice.** A gene has many isoforms and a protein coordinate means
nothing without one. The only defensible rule is to require the annotated
reference residue at the annotated position, which is what the pipeline's own
`ReferenceMismatch` guard enforces. Picking the longest isoform and hoping is
how a reconciliation silently succeeds against the wrong sequence.

**Retired symbols.** Published tables carry the symbols of their publication
year. `ACPP` is now `ACP3`, `MLL3` is `KMT2C`. Without following the rename a
gene simply vanishes and looks like a failed reconciliation. On the Ott
benchmark that step recovered 284 variants, one of them a validated positive,
where Ensembl's REST symbol lookup recovered none.

**Frameshifts cannot be checked without their tail.** A frameshift's novel
sequence is in no proteome, so the honest verdict for a published frameshift
epitope with no supplied tail is "unverifiable", not "wrong". Collapsing those
two would manufacture findings.

The generator is the pipeline's own `generate_for_variant`, deliberately. An
audit that reimplements the thing it audits measures the reimplementation.
"""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Iterable, Iterator, Sequence
from enum import StrEnum
from pathlib import Path

from pydantic import BaseModel, Field

from ..io.fasta import read_fasta
from ..io.ids import strip_version
from ..models import AMINO_ACIDS, Variant, VariantClass
from ..peptides.generate import PeptideGenerationError, generate_for_variant
from ..peptides.proteome import ProteomeIndex

_MISSENSE = re.compile(r"^p?\.?([A-Z])(\d+)([A-Z])$")
_FRAMESHIFT = re.compile(r"^p?\.?([A-Z])(\d+)(?:[A-Z]*)fs", re.IGNORECASE)
_DEL_ONE = re.compile(r"^p?\.?([A-Z])(\d+)del$")
_DEL_RANGE = re.compile(r"^p?\.?([A-Z])(\d+)_([A-Z])(\d+)del$")
_INS_RANGE = re.compile(r"^p?\.?([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$")

#: A deletion spanning more than this is almost certainly a parse error rather
#: than a real in-frame event, and resolving it would scan a huge span.
MAX_DELETION_SPAN = 30


class Verdict(StrEnum):
    DERIVABLE = "derivable"
    PEPTIDE_NOT_GENERATED = "peptide_not_generated"
    NO_CONSISTENT_ISOFORM = "no_consistent_isoform"
    GENE_NOT_FOUND = "gene_not_found"
    UNPARSEABLE_CHANGE = "unparseable_change"
    TAIL_REQUIRED = "tail_required"

    @property
    def is_verified(self) -> bool:
        return self is Verdict.DERIVABLE

    @property
    def is_unverifiable(self) -> bool:
        """Distinct from a failure: nothing here says the claim is wrong."""
        return self in (Verdict.TAIL_REQUIRED, Verdict.UNPARSEABLE_CHANGE)


class ProteinChange(BaseModel):
    start: int
    end: int
    aa_ref: str
    aa_alt: str
    variant_class: VariantClass


def parse_protein_change(change: str) -> ProteinChange | None:
    """`p.G261R` -> a missense at 261; `p.S754fs` -> a frameshift; None if unreadable.

    Frameshift is matched before missense because `p.S754Qfs*12` and friends
    would otherwise fall through to no match at all.
    """
    text = (change or "").strip()

    if match := _FRAMESHIFT.match(text):
        position = int(match.group(2))
        return ProteinChange(
            start=position,
            end=position,
            aa_ref=match.group(1),
            aa_alt="",
            variant_class=VariantClass.FRAMESHIFT,
        )
    if match := _MISSENSE.match(text):
        position = int(match.group(2))
        return ProteinChange(
            start=position,
            end=position,
            aa_ref=match.group(1),
            aa_alt=match.group(3),
            variant_class=VariantClass.MISSENSE,
        )
    if match := _DEL_RANGE.match(text):
        start, end = int(match.group(2)), int(match.group(4))
        if end < start or end - start > MAX_DELETION_SPAN:
            return None
        # Only the flanking residues are named; the interior is recovered from
        # the reference during isoform resolution.
        return ProteinChange(
            start=start, end=end, aa_ref="", aa_alt="", variant_class=VariantClass.INFRAME_DEL
        )
    if match := _DEL_ONE.match(text):
        position = int(match.group(2))
        return ProteinChange(
            start=position,
            end=position,
            aa_ref=match.group(1),
            aa_alt="",
            variant_class=VariantClass.INFRAME_DEL,
        )
    if match := _INS_RANGE.match(text):
        inserted = match.group(5)
        if not set(inserted) <= AMINO_ACIDS:
            return None
        position = int(match.group(2))
        return ProteinChange(
            start=position,
            end=position,
            aa_ref="",
            aa_alt=inserted,
            variant_class=VariantClass.INFRAME_INS,
        )
    return None


class SymbolIndex:
    """Every isoform of every gene symbol, for annotation-consistent lookup.

    `ProteomeIndex` keeps one sequence per symbol, which is right for ranking
    and wrong here: reconciliation has to try each isoform until one agrees
    with the annotated residue.
    """

    def __init__(self) -> None:
        self._by_symbol: dict[str, list[tuple[str, str, bool]]] = {}
        self._aliases: dict[str, str] = {}

    @classmethod
    def from_fasta(cls, path: Path) -> SymbolIndex:
        index = cls()
        for header, raw in read_fasta(Path(path)):
            sequence = raw.rstrip("*")
            if not sequence:
                continue
            tokens = header.split()
            attributes: dict[str, str] = {}
            for token in tokens[1:]:
                if ":" in token:
                    key, _, value = token.partition(":")
                    attributes.setdefault(key, value)
            transcript = strip_version(attributes.get("transcript"))
            symbol = attributes.get("gene_symbol")
            if not transcript or not symbol:
                continue
            # Primary assembly ahead of patch scaffolds, so a tie resolves to
            # the canonical-looking isoform rather than an arbitrary one.
            primary = len(tokens) > 2 and tokens[2].startswith("chromosome:")
            index._by_symbol.setdefault(symbol, []).append((transcript, sequence, primary))
        return index

    def load_aliases(self, path: Path) -> int:
        """Adopt a retired-symbol map, keeping only targets present in this FASTA."""
        mapping = json.loads(Path(path).read_text())
        self._aliases = {
            old: new
            for old, new in mapping.items()
            if new in self._by_symbol and old not in self._by_symbol
        }
        return len(self._aliases)

    def symbols(self) -> Iterable[str]:
        return self._by_symbol.keys()

    def resolve(self, symbol: str, change: ProteinChange) -> tuple[str, str] | None:
        """First isoform agreeing with the annotation, as `(transcript, sequence)`."""
        isoforms = self._by_symbol.get(symbol)
        if not isoforms:
            current = self._aliases.get(symbol)
            isoforms = self._by_symbol.get(current, []) if current else []
        if not isoforms:
            return None

        for transcript, sequence, _primary in sorted(
            isoforms, key=lambda item: (not item[2], -len(item[1]))
        ):
            if change.end > len(sequence):
                continue
            span = sequence[change.start - 1 : change.end]
            if change.aa_ref and span != change.aa_ref:
                continue
            if not span:
                continue
            return transcript, sequence
        return None


class EpitopeClaim(BaseModel):
    """One published claim: this peptide came from this change in this gene."""

    peptide: str
    gene: str
    protein_change: str
    source: str = ""
    downstream_protein: str | None = None


class Reconciliation(BaseModel):
    claim: EpitopeClaim
    verdict: Verdict
    transcript: str | None = None
    detail: str = ""

    def row(self) -> dict[str, str]:
        return {
            "peptide": self.claim.peptide,
            "gene": self.claim.gene,
            "protein_change": self.claim.protein_change,
            "source": self.claim.source,
            "verdict": self.verdict.value,
            "transcript": self.transcript or "",
            "detail": self.detail,
        }


class ReconciliationReport(BaseModel):
    results: list[Reconciliation] = Field(default_factory=list)

    @property
    def counts(self) -> dict[str, int]:
        counts = dict.fromkeys((v.value for v in Verdict), 0)
        for result in self.results:
            counts[result.verdict.value] += 1
        return counts

    def describe(self) -> str:
        total = len(self.results)
        if not total:
            return "no claims checked"
        verified = sum(1 for r in self.results if r.verdict.is_verified)
        unverifiable = sum(1 for r in self.results if r.verdict.is_unverifiable)
        contradicted = total - verified - unverifiable

        lines = [
            f"{total} claims checked",
            f"  derivable      {verified:>5}  ({verified / total:.1%})",
            f"  contradicted   {contradicted:>5}  ({contradicted / total:.1%})",
            f"  unverifiable   {unverifiable:>5}  ({unverifiable / total:.1%})",
            "",
        ]
        lines += [f"  {name:<24} {count}" for name, count in sorted(self.counts.items()) if count]
        if unverifiable:
            lines.append(
                "\nUnverifiable is not contradicted. A frameshift with no supplied tail "
                "cannot be checked against any proteome, and an unreadable change string "
                "says nothing about the claim."
            )
        return "\n".join(lines)


def reconcile_claim(
    claim: EpitopeClaim,
    index: SymbolIndex,
    lengths: Sequence[int] = (8, 9, 10, 11),
) -> Reconciliation:
    change = parse_protein_change(claim.protein_change)
    if change is None:
        return Reconciliation(
            claim=claim,
            verdict=Verdict.UNPARSEABLE_CHANGE,
            detail=f"could not read {claim.protein_change!r}",
        )

    if change.variant_class is VariantClass.FRAMESHIFT and not claim.downstream_protein:
        return Reconciliation(
            claim=claim,
            verdict=Verdict.TAIL_REQUIRED,
            detail="frameshift tail is in no proteome; supply downstream_protein to check",
        )

    resolved = index.resolve(claim.gene, change)
    if resolved is None:
        known = claim.gene in set(index.symbols())
        return Reconciliation(
            claim=claim,
            verdict=Verdict.NO_CONSISTENT_ISOFORM if known else Verdict.GENE_NOT_FOUND,
            detail=(
                f"no isoform of {claim.gene} carries {change.aa_ref!r} at {change.start}"
                if known
                else f"{claim.gene} is not in the proteome, before or after alias resolution"
            ),
        )

    transcript, sequence = resolved
    peptides = _generated_peptides(claim, change, transcript, sequence, lengths)
    if peptides is None:
        return Reconciliation(
            claim=claim,
            verdict=Verdict.NO_CONSISTENT_ISOFORM,
            transcript=transcript,
            detail="the resolved isoform could not produce peptides for this change",
        )

    if claim.peptide in peptides:
        return Reconciliation(claim=claim, verdict=Verdict.DERIVABLE, transcript=transcript)
    return Reconciliation(
        claim=claim,
        verdict=Verdict.PEPTIDE_NOT_GENERATED,
        transcript=transcript,
        detail=(
            f"{transcript} yields {len(peptides)} peptides overlapping "
            f"{claim.protein_change} and none is the claimed sequence"
        ),
    )


def _generated_peptides(
    claim: EpitopeClaim,
    change: ProteinChange,
    transcript: str,
    sequence: str,
    lengths: Sequence[int],
) -> set[str] | None:
    """Peptides the pipeline's own generator produces for this change."""
    index = ProteomeIndex()
    index._by_transcript = {transcript: sequence}
    index._by_protein = {transcript: sequence}

    variant = Variant(
        chrom="0",
        pos=1,
        ref="N",
        alt="N",
        gene=claim.gene,
        transcript=transcript,
        variant_class=change.variant_class,
        protein_start=change.start,
        protein_end=change.end,
        aa_ref=change.aa_ref or sequence[change.start - 1 : change.end],
        aa_alt=change.aa_alt,
        downstream_protein=claim.downstream_protein,
    )
    requested = sorted({len(claim.peptide), *lengths})
    try:
        candidates = generate_for_variant(variant, index, requested)
    except PeptideGenerationError:
        return None
    return {c.mutant_peptide for c in candidates}


def reconcile(
    claims: Iterable[EpitopeClaim],
    index: SymbolIndex,
    lengths: Sequence[int] = (8, 9, 10, 11),
) -> ReconciliationReport:
    return ReconciliationReport(
        results=[reconcile_claim(claim, index, lengths) for claim in claims]
    )


def read_claims(path: Path) -> Iterator[EpitopeClaim]:
    """Read claims from a TSV.

    Required columns: `peptide`, `gene`, `protein_change`.
    Optional: `source` (a citation or row id) and `downstream_protein`.
    """
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        available = {name.lower() for name in (reader.fieldnames or [])}
        missing = [c for c in ("peptide", "gene", "protein_change") if c not in available]
        if missing:
            raise ValueError(f"{path}: missing required columns {missing}")
        for raw in reader:
            row = {k.lower(): (v or "").strip() for k, v in raw.items() if k}
            if not row.get("peptide"):
                continue
            yield EpitopeClaim(
                peptide=row["peptide"],
                gene=row["gene"],
                protein_change=row["protein_change"],
                source=row.get("source", ""),
                downstream_protein=row.get("downstream_protein") or None,
            )


def write_report(report: ReconciliationReport, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = ("peptide", "gene", "protein_change", "source", "verdict", "transcript", "detail")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for result in report.results:
            writer.writerow(result.row())
    return path
