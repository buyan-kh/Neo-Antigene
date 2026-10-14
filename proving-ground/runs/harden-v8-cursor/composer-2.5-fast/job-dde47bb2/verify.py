import csv
import json
from pathlib import Path

from enumerate_downstream import (
    CODON_TABLE,
    read_fasta,
    read_variants,
    coding_sequence,
    downstream_peptide,
    translate,
)


def correct_downstream_peptide(cds, cds_pos_1based, ref, alt):
    zero_based = cds_pos_1based - 1
    mutant = cds[:zero_based] + alt + cds[zero_based + len(ref) :]
    codon_index = zero_based // 3
    tail = mutant[codon_index * 3 :]
    return translate(tail), codon_index


def parse_gff3_cds(path):
    spans = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if fields[2] != "CDS":
                continue
            attrs = dict(
                item.split("=", 1) for item in fields[8].split(";") if "=" in item
            )
            tid = attrs["Parent"].replace("transcript:", "")
            genomic_start, genomic_end = int(fields[3]), int(fields[4])
            spans[tid] = (genomic_start, genomic_end)
    return spans


def cds_from_gff3(transcript_id, transcript_seq, genomic_origin):
    gstart, gend = parse_gff3_cds("annotation.gff3")[transcript_id]
    tstart = gstart - genomic_origin
    tend = gend - genomic_origin
    return transcript_seq[tstart - 1 : tend]


def main():
    base = Path(__file__).parent
    transcripts = read_fasta(str(base / "transcripts.fasta"))
    variants = read_variants(str(base / "variants.vcf"))
    reported = {}
    with open(base / "reported_neoorfs.csv") as f:
        for row in csv.DictReader(f):
            reported[row["variant_id"]] = (
                row["downstream_peptide"],
                int(row["start_offset"]),
            )

    affected = []
    corrected = []
    for v in variants:
        tid = v["transcript"]
        cds_cached = coding_sequence(tid, transcripts[tid])
        cds_gff = cds_from_gff3(tid, transcripts[tid], 1000000)
        assert cds_cached == cds_gff, (tid, len(cds_cached), len(cds_gff))

        pipe_pep, pipe_off = downstream_peptide(
            cds_cached, v["cds_pos"], v["ref"], v["alt"]
        )
        good_pep, good_off = correct_downstream_peptide(
            cds_cached, v["cds_pos"], v["ref"], v["alt"]
        )
        rep_pep, rep_off = reported[v["variant_id"]]

        match_pipe = pipe_pep == rep_pep and pipe_off == rep_off
        match_good = good_pep == rep_pep and good_off == rep_off
        if not match_good:
            affected.append(v["variant_id"])
            corrected.append(
                {
                    "variant_id": v["variant_id"],
                    "downstream_peptide": good_pep,
                    "start_offset": good_off,
                }
            )
        print(
            v["variant_id"],
            "cds_pos",
            v["cds_pos"],
            "rep_off",
            rep_off,
            "pipe_off",
            pipe_off,
            "good_off",
            good_off,
            "rep==pipe",
            match_pipe,
            "rep==good",
            match_good,
        )
        if rep_pep != good_pep:
            print("  rep:", rep_pep[:60], "..." if len(rep_pep) > 60 else "")
            print("  good:", good_pep[:60], "..." if len(good_pep) > 60 else "")

    out = {
        "drift_detected": bool(affected),
        "mechanism": "downstream_start_skips_variant_codon"
        if affected
        else "none",
        "magnitude": {"offset": 1, "direction": "downstream"} if affected else {"offset": 0, "direction": "none"},
        "affected_variants": sorted(affected),
        "corrected": sorted(corrected, key=lambda x: x["variant_id"]),
    }
    out_path = base / "output" / "submission.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print("written", out_path)


if __name__ == "__main__":
    main()
