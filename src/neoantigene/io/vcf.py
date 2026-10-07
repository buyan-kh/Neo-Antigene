"""Somatic VCF reader for VEP-annotated calls.

Neo Antigene consumes protein-level consequences rather than re-deriving them,
so the input VCF must carry a VEP `CSQ` INFO field. Run VEP with at least:

    vep --input_file somatic.vcf --format vcf --vcf --symbol --terms SO \
        --canonical --biotype --transcript_version --cache --offline

A plain TSV path is also supported for teams that annotate elsewhere; see
`read_variant_tsv`.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import IO

from ..models import Variant, VariantClass

_CSQ_FORMAT_RE = re.compile(r"Format:\s*([^\"']+)")
#: `ensembl=114` is the release branch Downstream.pm was taken from. Fall back
#: to the major of `##VEP="v114.2"` when the cache line is absent. A word
#: boundary keeps `ensembl-io=114` from matching.
_VEP_ENSEMBL_RE = re.compile(r"\bensembl=(\d+)")
_VEP_VERSION_RE = re.compile(r'##VEP="v(\d+)')

_CONSEQUENCE_MAP = {
    "missense_variant": VariantClass.MISSENSE,
    "inframe_insertion": VariantClass.INFRAME_INS,
    "inframe_deletion": VariantClass.INFRAME_DEL,
    "frameshift_variant": VariantClass.FRAMESHIFT,
    "protein_altering_variant": VariantClass.OTHER,
    "stop_lost": VariantClass.OTHER,
    "start_lost": VariantClass.OTHER,
}

CODING_CLASSES = frozenset(
    {
        VariantClass.MISSENSE,
        VariantClass.INFRAME_INS,
        VariantClass.INFRAME_DEL,
        VariantClass.FRAMESHIFT,
    }
)

_TSV_REQUIRED = ("chrom", "pos", "ref", "alt", "transcript", "protein_position", "amino_acids")

#: CSQ keys that may carry a frameshift's novel tail, most specific first.
#: `DownstreamProtein` is VEP's Downstream plugin, which ships in
#: Ensembl/VEP_plugins. `FrameshiftSequence` is pVACtools' Frameshift.pm, which
#: does not — it ships inside pVACtools — and which emits the whole mutant
#: protein from residue 1 rather than the tail, so it is sliced on read.
_DOWNSTREAM_KEYS = ("DownstreamProtein", "FrameshiftSequence")
_WHOLE_PROTEIN_KEYS = frozenset({"FrameshiftSequence"})


class VCFError(ValueError):
    pass


@contextmanager
def _open(path: Path) -> Iterator[IO[str]]:
    handle: IO[str] = (
        gzip.open(path, "rt")  # noqa: SIM115 - closed by this context manager
        if str(path).endswith(".gz")
        else open(path)  # noqa: SIM115 - closed by this context manager
    )
    try:
        yield handle
    finally:
        handle.close()


def _classify(consequence_field: str) -> VariantClass:
    terms = consequence_field.split("&")
    for term in terms:
        mapped = _CONSEQUENCE_MAP.get(term)
        if mapped is not None and mapped in CODING_CLASSES:
            return mapped
    for term in terms:
        if term in _CONSEQUENCE_MAP:
            return _CONSEQUENCE_MAP[term]
    return VariantClass.OTHER


def _parse_protein_position(value: str) -> tuple[int, int] | None:
    if not value or value == "-":
        return None
    head = value.split("/", 1)[0]
    parts = head.split("-")
    try:
        start = int(parts[0])
    except ValueError:
        return None
    if len(parts) == 1:
        return start, start
    try:
        end = int(parts[1])
    except ValueError:
        return start, start
    return start, max(start, end)


def _parse_amino_acids(value: str) -> tuple[str, str]:
    if not value or value == "-":
        return "", ""
    if "/" not in value:
        return value, value
    ref, alt = value.split("/", 1)
    return ("" if ref == "-" else ref), ("" if alt == "-" else alt)


def _parse_info(info_field: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for entry in info_field.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            info[key] = value
        else:
            info[entry] = "true"
    return info


def _extract_vaf(fmt: str, sample: str) -> tuple[float | None, int | None]:
    record = dict(zip(fmt.split(":"), sample.split(":"), strict=False))
    depth = _maybe_int(record.get("DP"))
    raw_ad = record.get("AD")
    if raw_ad and raw_ad != ".":
        try:
            counts = [int(x) for x in raw_ad.split(",")]
        except ValueError:
            counts = []
        total = sum(counts)
        if total > 0 and len(counts) >= 2:
            return counts[1] / total, depth if depth is not None else total
    affinity = _maybe_float(record.get("AF"))
    return affinity, depth


def _population_af(info: dict[str, str]) -> float | None:
    for key in ("gnomAD_AF", "MAX_AF", "AF_popmax", "POPAF"):
        value = _maybe_float(info.get(key))
        if value is not None:
            return value
    return None


def _csq_fields(header_lines: Sequence[str]) -> list[str]:
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            match = _CSQ_FORMAT_RE.search(line)
            if match:
                return [f.strip() for f in match.group(1).rstrip('">').split("|")]
    raise VCFError(
        "no CSQ INFO header found; annotate the VCF with VEP before ranking "
        "(see neoantigene.io.vcf docstring)"
    )


def _pick_transcript(entries: list[dict[str, str]]) -> dict[str, str] | None:
    coding = [e for e in entries if _classify(e.get("Consequence", "")) in CODING_CLASSES]
    if not coding:
        return None
    for entry in coding:
        if entry.get("CANONICAL") == "YES":
            return entry
    protein_coding = [e for e in coding if e.get("BIOTYPE") == "protein_coding"]
    return (protein_coding or coding)[0]


def read_vep_vcf(
    path: Path,
    tumor_sample: str | None = None,
    all_transcripts: bool = False,
) -> Iterator[Variant]:
    """Yield protein-altering `Variant` records from a VEP-annotated somatic VCF."""
    header_lines: list[str] = []
    columns: list[str] = []
    csq_keys: list[str] = []
    vep_release: int | None = None

    with _open(Path(path)) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                columns = line.lstrip("#").split("\t")
                csq_keys = _csq_fields(header_lines)
                vep_release = _vep_release(header_lines)
                continue
            if not line:
                continue
            if not columns:
                raise VCFError(f"{path}: data line before #CHROM header")

            fields = line.split("\t")
            row = dict(zip(columns, fields, strict=False))
            info = _parse_info(row.get("INFO", ""))
            raw_csq = info.get("CSQ")
            if not raw_csq:
                continue

            entries = [
                dict(zip(csq_keys, block.split("|"), strict=False)) for block in raw_csq.split(",")
            ]
            selected = entries if all_transcripts else [_pick_transcript(entries)]

            vaf, depth = _sample_metrics(columns, fields, row, tumor_sample)
            population_af = _population_af(info)
            filters = tuple(f for f in (row.get("FILTER") or "").split(";") if f)

            for entry in selected:
                if not entry:
                    continue
                variant = _build_variant(
                    row, entry, vaf, depth, population_af, filters, vep_release
                )
                if variant is not None:
                    yield variant


def _vep_release(header_lines: Sequence[str]) -> int | None:
    """Ensembl release recorded on the VEP header, if the VCF still has one."""
    for line in header_lines:
        if not line.startswith("##VEP="):
            continue
        ensembl = _VEP_ENSEMBL_RE.search(line)
        if ensembl:
            return int(ensembl.group(1))
        version = _VEP_VERSION_RE.search(line)
        if version:
            return int(version.group(1))
    return None


def _build_variant(
    row: dict[str, str],
    entry: dict[str, str],
    vaf: float | None,
    depth: int | None,
    population_af: float | None,
    filters: tuple[str, ...],
    vep_release: int | None,
) -> Variant | None:
    variant_class = _classify(entry.get("Consequence", ""))
    if variant_class not in CODING_CLASSES:
        return None
    protein_position = _parse_protein_position(entry.get("Protein_position", ""))
    if protein_position is None:
        return None
    aa_ref, aa_alt = _parse_amino_acids(entry.get("Amino_acids", ""))
    tail, length_change = _downstream_annotation(entry, protein_position[0])
    return Variant(
        chrom=row["CHROM"],
        pos=int(row["POS"]),
        ref=row["REF"],
        alt=row["ALT"].split(",")[0],
        gene=entry.get("SYMBOL") or entry.get("Gene") or "",
        transcript=entry.get("Feature") or "",
        variant_class=variant_class,
        protein_start=protein_position[0],
        protein_end=protein_position[1],
        aa_ref=aa_ref,
        aa_alt=aa_alt,
        dna_vaf=vaf,
        tumor_depth=depth,
        population_af=population_af,
        filters=filters,
        downstream_protein=tail,
        protein_length_change=length_change,
        vep_release=vep_release,
    )


def _downstream_annotation(
    entry: dict[str, str], protein_start: int
) -> tuple[str | None, int | None]:
    """Novel frameshift tail, and a length change only when the tail is VEP's.

    A whole-protein field (`FrameshiftSequence`) is sliced to start at
    `protein_start`. Its length change is dropped: pVACtools does not emit
    one, and VEP's number describes `DownstreamProtein`, not this slice.
    `ProteinLengthChange` is kept only for `DownstreamProtein`, which is the
    field it was computed from.
    """
    for key in _DOWNSTREAM_KEYS:
        value = (entry.get(key) or "").strip()
        if not value or value == "-":
            continue
        if key in _WHOLE_PROTEIN_KEYS:
            return value[protein_start - 1 :] or None, None
        return value, _maybe_int(entry.get("ProteinLengthChange"))
    return None, None


def _sample_metrics(
    columns: list[str],
    fields: list[str],
    row: dict[str, str],
    tumor_sample: str | None,
) -> tuple[float | None, int | None]:
    if "FORMAT" not in columns:
        return None, None
    sample_columns = columns[columns.index("FORMAT") + 1 :]
    if not sample_columns:
        return None, None
    if tumor_sample is None:
        name = sample_columns[-1]
    elif tumor_sample in sample_columns:
        name = tumor_sample
    else:
        raise VCFError(f"tumor sample {tumor_sample!r} not in VCF samples {sample_columns}")
    return _extract_vaf(row["FORMAT"], fields[columns.index(name)])


def read_variant_tsv(path: Path) -> Iterator[Variant]:
    """Read pre-annotated variants from a TSV.

    Required columns: chrom, pos, ref, alt, transcript, protein_position,
    amino_acids (`G/D`), consequence.
    Optional: gene, dna_vaf, rna_vaf, tumor_depth, copy_number, population_af, filter,
    downstream_protein, protein_length_change, vep_release.

    `downstream_protein` is the novel tail a frameshift produces. Frameshift
    rows without it are read and reported but yield no peptides. Supply
    `vep_release` of 114 or later together with `protein_length_change` when
    the tail is VEP `DownstreamProtein` from that release; the generator uses
    the pair to place the tail. Either column alone leaves the join at
    `protein_position`.
    """
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        available = {name.lower() for name in (reader.fieldnames or [])}
        missing = [c for c in _TSV_REQUIRED if c not in available]
        if missing:
            raise VCFError(f"{path}: missing required columns {missing}")
        for raw_row in reader:
            record = {k.lower(): (v or "").strip() for k, v in raw_row.items()}
            protein_position = _parse_protein_position(record["protein_position"])
            if protein_position is None:
                continue
            aa_ref, aa_alt = _parse_amino_acids(record["amino_acids"])
            yield Variant(
                chrom=record["chrom"],
                pos=int(record["pos"]),
                ref=record["ref"],
                alt=record["alt"],
                gene=record.get("gene", ""),
                transcript=record["transcript"],
                variant_class=_classify(record.get("consequence", "")),
                protein_start=protein_position[0],
                protein_end=protein_position[1],
                aa_ref=aa_ref,
                aa_alt=aa_alt,
                dna_vaf=_maybe_float(record.get("dna_vaf")),
                rna_vaf=_maybe_float(record.get("rna_vaf")),
                tumor_depth=_maybe_int(record.get("tumor_depth")),
                copy_number=_maybe_float(record.get("copy_number")),
                population_af=_maybe_float(record.get("population_af")),
                filters=tuple(f for f in (record.get("filter") or "").split(";") if f),
                downstream_protein=(record.get("downstream_protein") or "").strip() or None,
                protein_length_change=_maybe_int(record.get("protein_length_change")),
                vep_release=_maybe_int(record.get("vep_release")),
            )


def _maybe_float(value: str | None) -> float | None:
    if value is None or value in ("", "."):
        return None
    try:
        return float(value.split(",")[0])
    except ValueError:
        return None


def _maybe_int(value: str | None) -> int | None:
    parsed = _maybe_float(value)
    return None if parsed is None else int(parsed)
