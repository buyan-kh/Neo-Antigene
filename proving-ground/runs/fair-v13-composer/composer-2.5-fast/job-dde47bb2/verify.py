#!/usr/bin/env python3
import csv
import json
from pathlib import Path

from cds_offsets import CDS_SPAN

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

ROOT = Path(__file__).resolve().parent


def read_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def parse_gff(path, origin):
    """Return per-transcript CDS and 3' UTR spans in 1-based transcript coords."""
    spans = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in ("CDS", "three_prime_UTR"):
                continue
            gstart, gend = int(fields[3]), int(fields[4])
            parent = fields[8]
            tid = None
            for part in parent.split(";"):
                if part.startswith("Parent=transcript:"):
                    tid = part.split(":", 1)[1]
            if not tid:
                continue
            tstart = gstart - origin
            tend = gend - origin
            spans.setdefault(tid, {})
            spans[tid][feature] = (tstart, tend)
    return spans


def translate_until_stop(nucleotides):
    residues = []
    for index in range(0, len(nucleotides) - 2, 3):
        residue = CODON_TABLE.get(nucleotides[index:index + 3], "X")
        residues.append(residue)
        if residue == "*":
            break
    return "".join(residues)


def read_variants(path):
    variants = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(item.split("=", 1) for item in fields[7].split(";"))
            variants.append(
                {
                    "variant_id": fields[2],
                    "ref": fields[3],
                    "alt": fields[4],
                    "transcript": info["TRANSCRIPT"],
                    "cds_pos": int(info["CDSPOS"]),
                }
            )
    return variants


def correct_downstream(transcript_seq, cds_span, utr3_span, cds_pos_1based, ref, alt):
    cds_start, cds_end = cds_span
    cds = transcript_seq[cds_start - 1 : cds_end]
    zero_based = cds_pos_1based - 1
    mutant_cds = cds[:zero_based] + alt + cds[zero_based + len(ref) :]
    first_codon = zero_based // 3
    tail_cds = mutant_cds[first_codon * 3 :]
    utr3 = ""
    if utr3_span:
        ustart, uend = utr3_span
        utr3 = transcript_seq[ustart - 1 : uend]
    tail = tail_cds + utr3
    return translate_until_stop(tail), first_codon


def pipeline_downstream(transcript_seq, transcript_id, cds_pos_1based, ref, alt):
    start, end = CDS_SPAN[transcript_id]
    cds = transcript_seq[start - 1 : end]
    zero_based = cds_pos_1based - 1
    mutant = cds[:zero_based] + alt + cds[zero_based + len(ref) :]
    first_codon = zero_based // 3
    tail = mutant[first_codon * 3 :]
    return translate_until_stop(tail), first_codon


def main():
    origin = 1_000_000
    for line in open(ROOT / "annotation.gff3"):
        if line.startswith("##transcript-origin"):
            origin = int(line.strip().split()[-1])
            break

    gff_spans = parse_gff(ROOT / "annotation.gff3", origin)
    transcripts = read_fasta(ROOT / "transcripts.fasta")
    variants = read_variants(ROOT / "variants.vcf")

    reported = {}
    with open(ROOT / "reported_neoorfs.csv") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            reported[row["variant_id"]] = (
                row["downstream_peptide"],
                int(row["start_offset"]),
            )

    span_diffs = []
    for tid in sorted(gff_spans):
        gff_cds = gff_spans[tid].get("CDS")
        cached = CDS_SPAN.get(tid)
        if gff_cds != cached:
            span_diffs.append((tid, gff_cds, cached))

    affected = []
    corrected = []
    codon_shift = None

    for v in variants:
        tid = v["transcript"]
        gff_cds = gff_spans[tid]["CDS"]
        utr3 = gff_spans[tid].get("three_prime_UTR")
        seq = transcripts[tid]
        pep, off = correct_downstream(
            seq, gff_cds, utr3, v["cds_pos"], v["ref"], v["alt"]
        )
        rep_pep, rep_off = reported[v["variant_id"]]
        if pep != rep_pep or off != rep_off:
            affected.append(v["variant_id"])
            corrected.append(
                {
                    "variant_id": v["variant_id"],
                    "downstream_peptide": pep,
                    "start_offset": off,
                }
            )
            cached_start = CDS_SPAN[tid][0]
            ann_start = gff_cds[0]
            delta = (ann_start - cached_start) // 3
            if codon_shift is None:
                codon_shift = delta
            elif codon_shift != delta:
                print(f"WARNING inconsistent shift for {v['variant_id']}: {delta} vs {codon_shift}")

    drift = len(affected) > 0
    if not drift:
        mechanism = "none"
        magnitude = {"offset": 0, "direction": "none"}
    else:
        mechanism = "annotation_release_mismatch"
        if codon_shift is None:
            codon_shift = 0
        if codon_shift > 0:
            magnitude = {"offset": abs(codon_shift), "direction": "downstream"}
        elif codon_shift < 0:
            magnitude = {"offset": abs(codon_shift), "direction": "upstream"}
        else:
            magnitude = {"offset": 0, "direction": "none"}

    out = {
        "drift_detected": drift,
        "mechanism": mechanism,
        "magnitude": magnitude,
        "affected_variants": sorted(affected),
        "corrected": sorted(corrected, key=lambda x: x["variant_id"]),
    }

    print("Span diffs (gff vs cache):", span_diffs)
    print(json.dumps(out, indent=2))

    out_path = ROOT / "output" / "submission.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
