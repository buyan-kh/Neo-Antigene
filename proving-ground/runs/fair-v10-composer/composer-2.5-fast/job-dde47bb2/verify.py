#!/usr/bin/env python3
import csv
import re

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

GENOMIC_ORIGIN = 1_000_000


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


def parse_gff3_cds(path):
    spans = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = dict(
                item.split("=", 1) for item in fields[8].split(";") if "=" in item
            )
            parent = attrs.get("Parent", "")
            m = re.search(r"ENST\d+", parent)
            if not m:
                continue
            tid = m.group(0)
            gstart, gend = int(fields[3]), int(fields[4])
            tstart = gstart - GENOMIC_ORIGIN
            tend = gend - GENOMIC_ORIGIN
            spans[tid] = (tstart, tend)
    return spans


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


def coding_sequence_from_span(transcript_seq, start, end):
    return transcript_seq[start - 1 : end]


def downstream_peptide(cds, cds_pos_1based, ref, alt):
    zero_based = cds_pos_1based - 1
    mutant = cds[:zero_based] + alt + cds[zero_based + len(ref) :]
    first_codon = zero_based // 3
    tail = mutant[first_codon * 3 :]
    return translate(tail), first_codon


def main():
    from cds_offsets import CDS_SPAN as CACHED

    transcripts = read_fasta("transcripts.fasta")
    gff_spans = parse_gff3_cds("annotation.gff3")
    variants = read_variants("variants.vcf")

    reported = {}
    with open("reported_neoorfs.csv") as f:
        r = csv.DictReader(f)
        for row in r:
            reported[row["variant_id"]] = (
                row["downstream_peptide"],
                int(row["start_offset"]),
            )

    print("CDS span comparison (cached vs gff3):")
    for tid in sorted(gff_spans):
        g = gff_spans[tid]
        c = CACHED[tid]
        diff = g[0] - c[0]
        print(f"  {tid}: cached={c} gff={g} start_diff={diff}")

    corrected = []
    affected = []
    for v in variants:
        tid = v["transcript"]
        seq = transcripts[tid]
        gstart, gend = gff_spans[tid]
        cds = coding_sequence_from_span(seq, gstart, gend)
        pep, off = downstream_peptide(cds, v["cds_pos"], v["ref"], v["alt"])
        vid = v["variant_id"]
        rep_pep, rep_off = reported[vid]
        ok = pep == rep_pep and off == rep_off
        print(f"{vid}: match={ok}")
        if not ok:
            affected.append(vid)
            print(f"  reported: off={rep_off} pep={rep_pep[:60]}...")
            print(f"  correct:  off={off} pep={pep}")
            corrected.append(
                {"variant_id": vid, "downstream_peptide": pep, "start_offset": off}
            )

    print("\naffected:", affected)


if __name__ == "__main__":
    main()
