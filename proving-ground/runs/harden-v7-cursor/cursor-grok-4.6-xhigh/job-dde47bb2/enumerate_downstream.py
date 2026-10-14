"""Downstream frameshift peptide enumeration.

Reduced reimplementation of the VEP Downstream plugin logic, operating in
transcript space. CDS spans are cached in cds_offsets.py so the annotation
does not have to be parsed on every run.
"""

import csv
import sys

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
        residue = CODON_TABLE.get(nucleotides[index:index + 3], "X")
        residues.append(residue)
        if residue == "*":
            break
    return "".join(residues)


def coding_sequence(transcript_id, transcript_seq):
    start, end = CDS_SPAN[transcript_id]
    return transcript_seq[start - 1:end]


def downstream_peptide(cds, cds_pos_1based, ref, alt):
    """Return the frameshifted tail and the zero-based codon it starts at."""
    zero_based = cds_pos_1based - 1
    mutant = cds[:zero_based] + alt + cds[zero_based + len(ref):]
    first_codon = zero_based // 3 + 1
    tail = mutant[first_codon * 3:]
    return translate(tail), first_codon


def main(fasta_path, vcf_path, out_path):
    transcripts = read_fasta(fasta_path)
    rows = []
    for variant in read_variants(vcf_path):
        cds = coding_sequence(variant["transcript"], transcripts[variant["transcript"]])
        peptide, offset = downstream_peptide(
            cds, variant["cds_pos"], variant["ref"], variant["alt"]
        )
        rows.append((variant["variant_id"], peptide, offset))
    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant_id", "downstream_peptide", "start_offset"])
        writer.writerows(sorted(rows))


if __name__ == "__main__":
    main(*sys.argv[1:4])
