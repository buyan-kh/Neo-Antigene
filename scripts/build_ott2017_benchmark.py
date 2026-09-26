#!/usr/bin/env python
"""Build the Ott et al. 2017 benchmark case from the paper's own supplements.

    uv run --extra plots python scripts/build_ott2017_benchmark.py

Downloads three supplementary tables from Springer's CDN, verifies them against
pinned SHA-256 digests, and emits pipeline inputs plus a separate label file.

The separation is the point. Everything under `inputs/` is derived only from
Supplementary Table 2 (the full somatic mutation list) and Table 4 (HLA
typing). The ELISPOT columns of Table 5 are read only when writing
`validated.tsv`, which the pipeline never reads. A reviewer can check that
claim by grepping this file for the assay column names: they appear once.

Source
    Ott PA, Hu Z, Keskin DB, et al. "An immunogenic personal neoantigen
    vaccine for patients with melanoma." Nature 2017;547(7662):217-221.
    doi:10.1038/nature22991

Coordinate harmonization, stated plainly. The MAF is on NCBI build 37 with
UCSC transcript identifiers; peptides are built against Ensembl 116 / GRCh38.
For each variant this script selects an Ensembl isoform of the same gene whose
residue at the annotated position equals the annotated reference residue.
Variants with no consistent isoform are dropped and counted. That is a
resolution step on the inputs, not a filter on outcomes, and it cannot see
labels.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import sys
import urllib.request
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import time, timedelta
from pathlib import Path
from typing import Any, Final

REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[1]
OUT_ROOT: Final[Path] = REPO_ROOT / "data" / "benchmarks" / "ott2017"

DOI: Final[str] = "10.1038/nature22991"
PMCID: Final[str] = "PMC5577644"
CDN: Final[str] = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fnature22991/MediaObjects/41586_2017_BFnature22991_MOESM{n}_ESM.xlsx"
)

#: Verified on download, 2026-09-26. A changed digest means the publisher
#: replaced the file and the numbers in docs/BENCHMARK.md no longer describe it.
DIGESTS: Final[dict[int, str]] = {
    2: "bc537bca05393bad3c29afc04ea130cfe5b0d889dd9445985e3f3a4f6a83f00d",
    4: "c761b08f529ad21f1ff3c3a817a0187ee8913f5605687de4c38a893131e9b3b9",
    5: "9b10a6b6a430a331d88973889caed3cbe9328609774c936c7c2a851f41476b60",
}

#: Patients 1-6 received the vaccine and have ELISPOT data in Table 5.
VACCINATED: Final[tuple[int, ...]] = (1, 2, 3, 4, 5, 6)

#: MAF classes this pipeline can enumerate peptides for. Frameshift and
#: nonsense variants are excluded because neo-ORF enumeration is not
#: implemented; `docs/BENCHMARK.md` records what that costs.
SUPPORTED_CLASSES: Final[dict[str, str]] = {
    "Missense_Mutation": "missense_variant",
    "In_Frame_Del": "inframe_deletion",
    "In_Frame_Ins": "inframe_insertion",
}

AMINO_ACIDS: Final[frozenset[str]] = frozenset("ACDEFGHIKLMNPQRSTVWY")

_MISSENSE = re.compile(r"^p\.([A-Z])(\d+)([A-Z])$")
_DEL_ONE = re.compile(r"^p\.([A-Z])(\d+)del$")
_DEL_RANGE = re.compile(r"^p\.([A-Z])(\d+)_([A-Z])(\d+)del$")
_INS_RANGE = re.compile(r"^p\.([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$")


@dataclass
class Tally:
    """Every row this script touched, so nothing disappears unexplained."""

    counts: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    def hit(self, reason: str, n: int = 1) -> None:
        self.counts[reason] += n

    def report(self, title: str) -> str:
        lines = [title]
        for reason, count in sorted(self.counts.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {reason:<46} {count:>7}")
        return "\n".join(lines)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def download(number: int, cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / f"ott2017_S{number}.xlsx"
    if not target.exists():
        url = CDN.format(n=number)
        print(f"  downloading S{number} <- {url}")
        request = urllib.request.Request(url, headers={"User-Agent": "neoantigene-benchmark/0.1"})
        with urllib.request.urlopen(request, timeout=120) as response, open(target, "wb") as out:
            shutil.copyfileobj(response, out)

    actual = sha256(target)
    expected = DIGESTS[number]
    if actual != expected:
        raise SystemExit(
            f"S{number} digest mismatch\n  expected {expected}\n  actual   {actual}\n"
            "The publisher's file changed. Re-verify before trusting any benchmark number."
        )
    print(f"  S{number} verified sha256={actual[:16]}...")
    return target


def load_sheet(path: Path) -> list[list[Any]]:
    import openpyxl

    workbook = openpyxl.load_workbook(path, data_only=True, read_only=True)
    sheet = workbook[workbook.sheetnames[0]]
    return [list(row) for row in sheet.iter_rows(values_only=True)]


def decode_allele(value: Any) -> str | None:
    """Recover an HLA allele from Excel's time-coerced cells.

    Table 4 stores allotypes as times, so `A*24:02` arrives as a timedelta of
    one day plus two minutes. Total minutes recovers `field1:field2`. Cells
    that are already text (`04:01:01G`) are passed through.
    """
    if value is None:
        return None
    if isinstance(value, str):
        parts = value.strip().split(":")
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            return f"{int(parts[0]):02d}:{int(parts[1]):02d}"
        return None
    if isinstance(value, time):
        total = value.hour * 60 + value.minute
    elif isinstance(value, timedelta):
        total = int(value.total_seconds() // 60)
    else:
        return None
    return f"{total // 60:02d}:{total % 60:02d}"


def read_hla(path: Path) -> dict[int, list[str]]:
    """Class I HLA-A and HLA-B per patient. Table 4 has no HLA-C column."""
    rows = load_sheet(path)
    header = [str(c) if c is not None else "" for c in rows[1]]
    genes = [h for h in header[1:]]

    typing: dict[int, list[str]] = {}
    for row in rows[2:]:
        if not row or not isinstance(row[0], int | float):
            continue
        patient = int(row[0])
        alleles: list[str] = []
        for gene, cell in zip(genes, row[1:], strict=False):
            if gene not in {"HLA-A", "HLA-B"}:
                continue
            decoded = decode_allele(cell)
            if decoded:
                alleles.append(f"HLA-{gene.removeprefix('HLA-')}*{decoded}")
        # Homozygous patients list the same allele twice.
        typing[patient] = sorted(dict.fromkeys(alleles))
    return typing


def parse_protein_change(change: str) -> tuple[int, int, str, str] | None:
    """`p.G261R` -> (261, 261, 'G', 'R'); deletions -> alt ''; None if unparseable."""
    change = (change or "").strip()

    if match := _MISSENSE.match(change):
        ref, position, alt = match.group(1), int(match.group(2)), match.group(3)
        return position, position, ref, alt
    if match := _DEL_ONE.match(change):
        ref, position = match.group(1), int(match.group(2))
        return position, position, ref, ""
    if match := _DEL_RANGE.match(change):
        start, end = int(match.group(2)), int(match.group(4))
        if end < start or end - start > 30:
            return None
        # Only the flanking residues are named; the interior is unknown here
        # and is recovered from the reference during isoform resolution.
        return start, end, "", ""
    if match := _INS_RANGE.match(change):
        start = int(match.group(2))
        inserted = match.group(5)
        if not set(inserted) <= AMINO_ACIDS:
            return None
        return start, start, "", inserted
    return None


#: Ensembl's REST endpoint resolves retired symbols through its own synonym
#: table, so the alias map is fetched rather than typed out here. The paper's
#: MAF is annotated against hg19/UCSC and carries 2013-era symbols: `MLL3` is
#: now `KMT2C`, `LPHN2` is `ADGRL2`, `ODZ3` is `TENM3`, `ACPP` is `ACP3`.
#: Without this step those genes silently vanish from the candidate pool.
HGNC_URL: Final[str] = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)


class Proteome:
    """Gene symbol to Ensembl isoforms, for build-37 to build-38 remapping."""

    def __init__(self, path: Path) -> None:
        self._by_symbol: dict[str, list[tuple[str, str, bool]]] = defaultdict(list)
        self._aliases: dict[str, str] = {}
        opener = gzip.open if path.suffix == ".gz" else open
        header: str | None = None
        chunks: list[str] = []
        with opener(path, "rt") as handle:  # type: ignore[operator]
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    self._add(header, "".join(chunks))
                    header, chunks = line[1:], []
                else:
                    chunks.append(line.strip())
        self._add(header, "".join(chunks))

    def _add(self, header: str | None, sequence: str) -> None:
        if header is None:
            return
        sequence = sequence.rstrip("*")
        if not sequence:
            return
        attributes: dict[str, str] = {}
        for token in header.split()[1:]:
            if ":" in token:
                key, _, value = token.partition(":")
                attributes.setdefault(key, value)
        symbol = attributes.get("gene_symbol")
        transcript = attributes.get("transcript")
        if not transcript:
            return
        # Prefer the primary assembly over patch scaffolds.
        primary = header.split()[2].startswith("chromosome:") if len(header.split()) > 2 else False
        entry = (transcript.split(".")[0], sequence, primary)
        if symbol:
            self._by_symbol[symbol].append(entry)

    def learn_aliases(self, symbols: Iterable[str], cache: Path) -> int:
        """Map symbols retired since 2017 onto their current approved symbol.

        Uses HGNC's own `prev_symbol` and `alias_symbol` columns. Ensembl's
        `lookup/symbol` endpoint was tried first and is not suitable: it does
        not consult the synonym table, so it returns nothing for exactly the
        renames that matter here (ACPP, ASNA1, ADCK3) while resolving a few
        hundred pseudogenes that have no protein sequence at all.

        Only renames pointing at a gene actually present in the pinned FASTA
        are kept, so a rename can never invent a lookup target.
        """
        wanted = {s for s in symbols if s and s not in self._by_symbol}
        if not wanted:
            return 0

        table = self._hgnc_table(cache)
        aliases: dict[str, str] = {}
        with open(table, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                current = (row.get("symbol") or "").strip()
                if not current or (row.get("status") or "").strip() != "Approved":
                    continue
                if current not in self._by_symbol:
                    continue
                for column in ("prev_symbol", "alias_symbol"):
                    for old in (row.get(column) or "").strip('"').split("|"):
                        old = old.strip()
                        # An approved symbol never yields to another gene's alias.
                        if old in wanted and old not in self._by_symbol:
                            aliases.setdefault(old, current)

        self._aliases = aliases
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(aliases, indent=2, sort_keys=True))
        return len(aliases)

    @staticmethod
    def _hgnc_table(cache: Path) -> Path:
        target = cache.with_name("hgnc_complete_set.txt")
        if not target.exists():
            print(f"  downloading HGNC symbol table <- {HGNC_URL}")
            target.parent.mkdir(parents=True, exist_ok=True)
            request = urllib.request.Request(
                HGNC_URL, headers={"User-Agent": "neoantigene-benchmark/0.1"}
            )
            with (
                urllib.request.urlopen(request, timeout=300) as response,
                open(target, "wb") as out,
            ):
                shutil.copyfileobj(response, out)
        print(f"  HGNC table sha256={sha256(target)[:16]}...")
        return target

    def resolve(
        self, symbol: str, start: int, end: int, aa_ref: str
    ) -> tuple[str, str, str] | None:
        """Pick an isoform consistent with the annotation.

        Returns `(transcript, aa_ref, aa_alt_reference_span)`. Candidates are
        ordered primary-assembly first, then longest, so a tie resolves to the
        canonical-looking isoform rather than an arbitrary one.
        """
        isoforms = self._by_symbol.get(symbol)
        if not isoforms:
            # The MAF was annotated in 2017; follow the HGNC rename if there is one.
            current = self._aliases.get(symbol)
            isoforms = self._by_symbol.get(current, []) if current else []
        candidates = sorted(isoforms, key=lambda item: (not item[2], -len(item[1])))
        for transcript, sequence, _primary in candidates:
            if end > len(sequence):
                continue
            span = sequence[start - 1 : end]
            if aa_ref and span != aa_ref:
                continue
            if not aa_ref and not span:
                continue
            return transcript, span, span
        return None


def maf_symbols(path: Path) -> list[str]:
    """Every gene symbol the mutation table mentions, for alias resolution."""
    rows = load_sheet(path)
    header = [str(c).strip() if c is not None else "" for c in rows[2]]
    column = header.index("Hugo Symbol")
    return [str(row[column]).strip() for row in rows[3:] if row and row[column]]


def read_variants(path: Path, proteome: Proteome, tally: Tally) -> dict[int, list[dict[str, str]]]:
    """Per-patient pipeline inputs from the full somatic mutation table."""
    rows = load_sheet(path)
    header = [str(c).strip() if c is not None else "" for c in rows[2]]
    index = {name: position for position, name in enumerate(header)}

    per_patient: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in rows[3:]:
        if not row or not isinstance(row[index["Patient ID"]], int | float):
            continue
        tally.hit("maf rows read")
        patient = int(row[index["Patient ID"]])
        if patient not in VACCINATED:
            tally.hit("dropped: patient not vaccinated")
            continue

        classification = str(row[index["Variant Classification"]] or "")
        consequence = SUPPORTED_CLASSES.get(classification)
        if consequence is None:
            tally.hit(f"dropped: unsupported class {classification}")
            continue

        parsed = parse_protein_change(str(row[index["Protein Change"]] or ""))
        if parsed is None:
            tally.hit("dropped: unparseable protein change")
            continue
        start, end, aa_ref, aa_alt = parsed

        symbol = str(row[index["Hugo Symbol"]] or "").strip()
        resolved = proteome.resolve(symbol, start, end, aa_ref)
        if resolved is None:
            tally.hit("dropped: no Ensembl isoform matches annotation")
            continue
        transcript, span, _ = resolved

        if consequence == "inframe_deletion":
            amino_acids = f"{span}/-"
        elif consequence == "inframe_insertion":
            amino_acids = f"-/{aa_alt}"
        else:
            amino_acids = f"{aa_ref}/{aa_alt}"

        tally.hit("kept")
        per_patient[patient].append(
            {
                "chrom": str(row[index["Chromosome"]]),
                "pos": str(int(row[index["Start position"]])),
                "ref": str(row[index["Reference Allele"]] or "N"),
                "alt": str(row[index["Tumor Seq Allele2"]] or "N"),
                "gene": symbol,
                "transcript": transcript,
                "consequence": consequence,
                "protein_position": str(start) if start == end else f"{start}-{end}",
                "amino_acids": amino_acids,
                "filter": "PASS",
            }
        )
    return per_patient


#: The only place the ELISPOT columns are read. Nothing under `inputs/`
#: derives from this function.
ELISPOT_PRIMARY: Final[str] = "Peptide pulsed autologous APC"
ELISPOT_MINIGENE: Final[str] = "Minigene expressing autologous B cell"
ELISPOT_TUMOR: Final[str] = "Autoloogus tumor"  # sic, as published


def read_labels(path: Path, tally: Tally) -> list[dict[str, str]]:
    """Assay labels from Table 5, in the schema `neoantigene.assays` reads.

    The call is taken from the peptide-pulsed autologous APC ELISPOT, which is
    the readout available for every row. The minigene and autologous-tumor
    columns are carried into `notes` as secondary evidence rather than merged
    into the call, because merging would invent a composite endpoint the paper
    did not report.
    """
    rows = load_sheet(path)
    # Two header rows: group labels on one, sub-labels on the next.
    upper = [str(c).strip() if c is not None else "" for c in rows[2]]
    lower = [str(c).strip() if c is not None else "" for c in rows[3]]

    column: dict[str, int] = {}
    for position, (top, bottom) in enumerate(zip(upper, lower, strict=False)):
        for name in (bottom, top):
            if name and name not in column:
                column[name] = position

    required = ["Patient ID", "Gene", "Protein change", "HLA allele", ELISPOT_PRIMARY]
    missing = [name for name in required if name not in column]
    if missing:
        raise SystemExit(f"Table 5 is missing expected columns: {missing}")
    peptide_column = column["Sequence"]

    labels: list[dict[str, str]] = []
    for row in rows[4:]:
        if not row or not isinstance(row[column["Patient ID"]], int | float):
            continue
        tally.hit("table 5 rows read")
        peptide = str(row[peptide_column] or "").strip().upper()
        if not peptide or not set(peptide) <= AMINO_ACIDS:
            tally.hit("dropped: no clean mutant peptide sequence")
            continue

        allele = str(row[column["HLA allele"]] or "").strip()
        if not allele:
            tally.hit("dropped: no HLA restriction")
            continue

        raw = row[column[ELISPOT_PRIMARY]]
        call = {1: "positive", 0: "negative"}.get(int(raw) if isinstance(raw, int | float) else -1)
        if call is None:
            tally.hit("dropped: ELISPOT not determined (n.d.)")
            continue
        tally.hit(f"label: {call}")

        labels.append(
            {
                "sample_id": f"OTT2017-PT{int(row[column['Patient ID']])}",
                "peptide": peptide,
                "allele": f"HLA-{allele}" if not allele.startswith("HLA-") else allele,
                "assay": "ifng_elispot",
                "call": call,
                "notes": (
                    f"gene={str(row[column['Gene']] or '').strip()};"
                    f"protein_change={str(row[column['Protein change']] or '').strip()};"
                    f"minigene={_flag(row, column, ELISPOT_MINIGENE)};"
                    f"autologous_tumor={_flag(row, column, ELISPOT_TUMOR)};"
                    f"source=doi:{DOI} Supplementary Table 5"
                ),
            }
        )
    return labels


def _flag(row: list[Any], column: dict[str, int], name: str) -> str:
    if name not in column:
        return "na"
    value = row[column[name]]
    if isinstance(value, int | float):
        return {1: "positive", 0: "negative"}.get(int(value), "nd")
    return "nd"


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(path: Path, patient: int, alleles: list[str], proteome: str) -> None:
    hla = "\n".join(f"  - {allele}" for allele in alleles)
    path.write_text(
        f"""# Ott et al. 2017 (doi:{DOI}), Patient {patient}.
# Generated by scripts/build_ott2017_benchmark.py -- do not hand-edit.
#
# No RNA expression table: the paper publishes per-peptide TPM only for
# vaccine-selected peptides, so using it would advantage exactly the peptides
# that were assayed. The expression and tumor_selectivity features are
# therefore neutral for every candidate in this benchmark.
#
# No tumor_purity: not published per patient, so clonality is neutral too.

sample_id: OTT2017-PT{patient}
cancer_type: melanoma

hla:
{hla}

processed:
  variant_tsv: variants.tsv
  proteome_fasta: {proteome}
"""
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT_ROOT)
    parser.add_argument(
        "--proteome",
        type=Path,
        default=None,
        help="Ensembl peptide FASTA. Defaults to the fetch-proteome cache.",
    )
    args = parser.parse_args()

    from neoantigene.peptides import reference

    proteome_path = args.proteome or reference.proteome_path()
    if not proteome_path.exists():
        raise SystemExit(
            f"reference proteome not found at {proteome_path}\n"
            "Run `uv run neoantigene fetch-proteome` first."
        )

    cache = args.out / ".cache"
    print(f"Ott et al. 2017 benchmark  doi:{DOI}  {PMCID}")
    print("1. fetching supplements")
    tables = {n: download(n, cache) for n in (2, 4, 5)}

    print(f"2. indexing {proteome_path.name}")
    proteome = Proteome(proteome_path)

    print("3. reading HLA typing (Table 4)")
    typing = read_hla(tables[4])
    for patient in VACCINATED:
        print(f"   PT{patient}: {' '.join(typing.get(patient, [])) or 'MISSING'}")

    print("4. resolving gene symbols retired since 2017")
    learned = proteome.learn_aliases(maf_symbols(tables[2]), cache / "symbol_aliases.json")
    print(f"   mapped {learned} retired symbols onto current Ensembl gene IDs")

    print("5. converting the full mutation table (Table 2)")
    variant_tally = Tally()
    per_patient = read_variants(tables[2], proteome, variant_tally)
    print(variant_tally.report("   mutation accounting:"))

    print("6. reading ELISPOT labels (Table 5)")
    label_tally = Tally()
    labels = read_labels(tables[5], label_tally)
    print(label_tally.report("   label accounting:"))

    columns = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "transcript",
        "consequence",
        "protein_position",
        "amino_acids",
        "filter",
    ]
    for patient in VACCINATED:
        rows = per_patient.get(patient, [])
        alleles = typing.get(patient, [])
        if not rows or not alleles:
            print(f"   PT{patient}: skipped ({len(rows)} variants, {len(alleles)} alleles)")
            continue
        directory = args.out / "inputs" / f"patient-{patient}"
        write_tsv(directory / "variants.tsv", rows, columns)
        write_manifest(directory / "sample.yaml", patient, alleles, str(proteome_path))
        print(f"   PT{patient}: {len(rows)} variants, {len(alleles)} alleles -> {directory}")

    label_path = args.out / "validated.tsv"
    write_tsv(
        label_path,
        labels,
        ["sample_id", "peptide", "allele", "assay", "call", "notes"],
    )
    positives = sum(1 for row in labels if row["call"] == "positive")
    print(f"6. wrote {len(labels)} labels ({positives} positive) -> {label_path}")

    digests = {n: sha256(path) for n, path in tables.items()}
    (args.out / "SOURCES.md").write_text(_sources(digests, typing, variant_tally, label_tally))
    print(f"7. wrote provenance -> {args.out / 'SOURCES.md'}")
    return 0


def _sources(
    digests: dict[int, str],
    typing: dict[int, list[str]],
    variants: Tally,
    labels: Tally,
) -> str:
    rows = "\n".join(
        f"| Supplementary Table {n} | `{CDN.format(n=n)}` | `{digest}` |"
        for n, digest in sorted(digests.items())
    )
    hla = "\n".join(
        f"| Patient {p} | {' '.join(typing.get(p, [])) or 'not published'} |" for p in VACCINATED
    )
    accounting = "\n".join(
        f"| {reason} | {count} |" for reason, count in sorted(variants.counts.items())
    )
    label_rows = "\n".join(
        f"| {reason} | {count} |" for reason, count in sorted(labels.counts.items())
    )
    return f"""# Benchmark provenance

Every peptide and label in this directory traces to the supplementary tables of
one paper. Regenerate with:

```bash
uv run python scripts/build_ott2017_benchmark.py
```

## Citation

Ott PA, Hu Z, Keskin DB, Shukla SA, Sun J, Bozym DJ, et al. "An immunogenic
personal neoantigen vaccine for patients with melanoma." *Nature*
2017;547(7662):217-221. doi:[{DOI}](https://doi.org/{DOI}). PMID 28678778.
{PMCID}.

## Files downloaded, with digests verified at build time

| Table | URL | SHA-256 |
| --- | --- | --- |
{rows}

Table 2 is the complete somatic mutation list for Patients 1-10. Table 4 is
HLA typing. Table 5 is the immunizing peptides with IFN-gamma ELISPOT results.

## Gene symbols retired since 2017

The MAF was annotated in 2017 and names genes that current Ensembl no longer
uses, so a naive symbol lookup silently discards them. Renames are followed
using HGNC's approved symbol table (`prev_symbol` and `alias_symbol`),
downloaded from
`https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`
and cached in `.cache/`. A rename is only accepted when the current symbol is
actually present in the pinned Ensembl release, so it can never invent a
lookup target.

This matters: without it, 284 supported-class variants are lost, among them
`ACPP` p.E34K in Patient 1 (now `ACP3`), whose peptide is one of the
experimentally validated positives. Ensembl's `lookup/symbol` REST endpoint
was tried first and rejected -- it does not consult the synonym table, so it
misses exactly these renames.

## What went into the pipeline, and what did not

`inputs/` derives **only** from Tables 2 and 4. `validated.tsv` derives from
Table 5's ELISPOT columns and is read only after ranking is complete.

Table 5 also publishes each peptide's predicted MHC affinity and per-gene TPM.
Neither is used: the affinities are another predictor's output, and the TPM
values exist only for vaccine-selected peptides, so using them would advantage
precisely the peptides that were assayed.

## HLA typing used

| Patient | Class I alleles |
| --- | --- |
{hla}

Supplementary Table 4 has HLA-A and HLA-B columns only, with no HLA-C, so
C-restricted candidates cannot be scored for this cohort. The allotypes are
stored as Excel time values and are decoded by total minutes; the decoded
values were cross-checked against the plain-text allele column in Table 5.

## Mutation accounting

| outcome | rows |
| --- | --- |
{accounting}

Frameshift and nonsense variants are excluded because neo-ORF peptide
enumeration is not implemented in this package. That exclusion is visible in
the benchmark as assayed peptides the ranker never scored.

## Label accounting

| outcome | rows |
| --- | --- |
{label_rows}

The call comes from the **peptide-pulsed autologous APC** IFN-gamma ELISPOT.
Rows marked `n.d.` are dropped rather than treated as negative. The minigene
and autologous-tumor readouts are preserved in the `notes` column as secondary
evidence, not merged into the call.

Peptides that were never included in a vaccine pool have no ELISPOT result and
are therefore **unlabeled, not negative**. Metrics are computed only over the
assayed intersection.
"""


if __name__ == "__main__":
    sys.exit(main())
