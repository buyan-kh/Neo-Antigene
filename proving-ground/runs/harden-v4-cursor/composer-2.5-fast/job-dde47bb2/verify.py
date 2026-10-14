#!/usr/bin/env python3
import csv
import json

from enumerate_downstream import (
    CODON_TABLE,
    read_fasta,
    read_variants,
    translate,
    downstream_peptide,
)
from cds_offsets import CDS_SPAN

CODON_TABLE_STD = CODON_TABLE


def parse_gff3_cds(path):
    spans = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if fields[2] != "CDS":
                continue
            tid = fields[0]
            start, end = int(fields[3]), int(fields[4])
            spans[tid] = (start, end)
    return spans


def coding_sequence_gff(transcript_id, transcript_seq, spans):
    start, end = spans[transcript_id]
    return transcript_seq[start - 1 : end]


def main():
    gff_spans = parse_gff3_cds("annotation.gff3")
    transcripts = read_fasta("transcripts.fasta")
    variants = read_variants("variants.vcf")

    print("CDS span mismatches (cached vs gff3):")
    for tid in sorted(gff_spans):
        g = gff_spans[tid]
        c = CDS_SPAN[tid]
        if g != c:
            print(f"  {tid}: cached={c} gff={g} start_delta={c[0]-g[0]} end_delta={c[1]-g[1]}")

    reported = {}
    with open("reported_neoorfs.csv") as f:
        for row in csv.DictReader(f):
            reported[row["variant_id"]] = (
                row["downstream_peptide"],
                int(row["start_offset"]),
            )

    corrected = []
    affected = []
    for v in variants:
        tid = v["transcript"]
        cds = coding_sequence_gff(tid, transcripts[tid], gff_spans)
        pep, off = downstream_peptide(cds, v["cds_pos"], v["ref"], v["alt"])
        vid = v["variant_id"]
        rep_pep, rep_off = reported[vid]
        ok = pep == rep_pep and off == rep_off
        if not ok:
            affected.append(vid)
            corrected.append(
                {"variant_id": vid, "downstream_peptide": pep, "start_offset": off}
            )
            print(f"\nMISMATCH {vid}")
            print(f"  reported: off={rep_off} pep={rep_pep[:60]}...")
            print(f"  correct:  off={off} pep={pep[:60]}...")

    # Check if errors match using cached vs gff cds
    print("\n--- Pipeline reproduce (cached spans) ---")
    for v in variants:
        tid = v["transcript"]
        cds = transcripts[tid][CDS_SPAN[tid][0] - 1 : CDS_SPAN[tid][1]]
        pep, off = downstream_peptide(cds, v["cds_pos"], v["ref"], v["alt"])
        vid = v["variant_id"]
        rep_pep, rep_off = reported[vid]
        print(f"{vid}: match_reported={pep==rep_pep and off==rep_off}")

    out = {
        "drift_detected": len(affected) > 0,
        "mechanism": "annotation_release_mismatch" if affected else "none",
        "magnitude": {"offset": 0, "direction": "none"},
        "affected_variants": sorted(affected),
        "corrected": sorted(corrected, key=lambda x: x["variant_id"]),
    }

    # Compute codon offset magnitude for start mismatch pattern
    if affected:
        deltas = []
        for v in variants:
            tid = v["transcript"]
            if gff_spans[tid] != CDS_SPAN[tid]:
                d = CDS_SPAN[tid][0] - gff_spans[tid][0]
                if d:
                    deltas.append(d // 3)
        if deltas and len(set(deltas)) == 1:
            codons = deltas[0]
            out["magnitude"] = {
                "offset": abs(codons),
                "direction": "upstream" if codons < 0 else "downstream",
            }

    with open("output/submission.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nWrote output/submission.json")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
