#!/usr/bin/env python3
import csv
import json
from pathlib import Path

ROOT = Path(__file__).parent

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


def parse_gff_cds(path, origin_by_tx):
    spans = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "CDS":
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            parent = attrs.get("Parent", "")
            if parent.startswith("transcript:"):
                tx = parent.split(":")[1]
            else:
                continue
            gstart, gend = int(fields[3]), int(fields[4])
            origin = origin_by_tx[tx]
            tstart = gstart - origin
            tend = gend - origin
            spans[tx] = (tstart, tend)
    return spans


def read_origin(path):
    origin = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            origin[row["transcript_id"]] = int(row["genomic_origin"])
    return origin


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


def translate(nucleotides):
    residues = []
    for index in range(0, len(nucleotides) - 2, 3):
        residue = CODON_TABLE.get(nucleotides[index : index + 3], "X")
        residues.append(residue)
        if residue == "*":
            break
    return "".join(residues)


def downstream_peptide(cds, cds_pos_1based, ref, alt):
    zero_based = cds_pos_1based - 1
    mutant = cds[:zero_based] + alt + cds[zero_based + len(ref) :]
    first_codon = zero_based // 3
    tail = mutant[first_codon * 3 :]
    return translate(tail), first_codon


def read_reported(path):
    reported = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            reported[row["variant_id"]] = (
                row["downstream_peptide"],
                int(row["start_offset"]),
            )
    return reported


def main():
    from cds_offsets import CDS_SPAN as CACHED

    origin = read_origin(ROOT / "transcript_origin.tsv")
    gff_spans = parse_gff_cds(ROOT / "annotation.gff3", origin)
    transcripts = read_fasta(ROOT / "transcripts.fasta")
    variants = read_variants(ROOT / "variants.vcf")
    reported = read_reported(ROOT / "reported_neoorfs.csv")

    print("CDS span comparison (gff vs cached):")
    for tx in sorted(gff_spans):
        g, c = gff_spans[tx], CACHED[tx]
        if g != c:
            print(f"  {tx}: gff={g} cached={c} start_delta={g[0]-c[0]} end_delta={g[1]-c[1]}")

    corrected = []
    affected = []
    for v in variants:
        tx = v["transcript"]
        start, end = gff_spans[tx]
        cds = transcripts[tx][start - 1 : end]
        pep, off = downstream_peptide(cds, v["cds_pos"], v["ref"], v["alt"])
        vid = v["variant_id"]
        rpep, roff = reported[vid]
        if pep != rpep or off != roff:
            affected.append(vid)
            corrected.append(
                {"variant_id": vid, "downstream_peptide": pep, "start_offset": off}
            )
            print(f"MISMATCH {vid}: reported=({roff},{rpep!r}) correct=({off},{pep!r})")
        else:
            print(f"OK {vid}")

    # magnitude: how many codons separate code CDS start from annotation
    # Check if consistent pattern
    start_deltas = []
    for tx in gff_spans:
        gstart = gff_spans[tx][0]
        cstart = CACHED[tx][0]
        start_deltas.append((gstart - cstart) // 3)

    print("Start codon deltas (gff - cached):", sorted(set(start_deltas)))

    drift = len(affected) > 0
    result = {
        "drift_detected": drift,
        "mechanism": "annotation_release_mismatch" if drift else "none",
        "magnitude": {"offset": 0, "direction": "none"},
        "affected_variants": sorted(affected),
        "corrected": sorted(corrected, key=lambda x: x["variant_id"]),
    }

    if drift:
        # Primary issue is cached CDS spans; quantify dominant start offset mismatch
        deltas_bp = [gff_spans[tx][0] - CACHED[tx][0] for tx in gff_spans]
        # use mode of non-zero if mixed
        from collections import Counter

        codon_deltas = Counter((gff_spans[tx][0] - CACHED[tx][0]) // 3 for tx in gff_spans)
        most_common = codon_deltas.most_common()
        print("Codon start deltas:", most_common)
        # annotation start relative to cached: if gff start > cached, annotation is downstream in transcript
        # magnitude.offset = codons between starts
        if len(set(deltas_bp)) == 1 and deltas_bp[0] != 0:
            off_codons = abs(deltas_bp[0]) // 3
            direction = "downstream" if deltas_bp[0] > 0 else "upstream"
            result["magnitude"] = {"offset": off_codons, "direction": direction}
        else:
            # mixed end/start issues - still annotation mismatch; use 1 codon as README pattern for starts
            plus_one = sum(1 for d in deltas_bp if d == 3)
            minus_one = sum(1 for d in deltas_bp if d == -3)
            zero = sum(1 for d in deltas_bp if d == 0)
            print(f"Mixed: +1codon={plus_one} -1codon={minus_one} zero={zero}")
            # For magnitude field when mixed, pick the most common start codon offset
            cd = codon_deltas.most_common(1)[0][0]
            if cd == 0:
                result["magnitude"] = {"offset": 0, "direction": "none"}
            elif cd > 0:
                result["magnitude"] = {"offset": cd, "direction": "downstream"}
            else:
                result["magnitude"] = {"offset": -cd, "direction": "upstream"}

    out = ROOT / "output" / "submission.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
