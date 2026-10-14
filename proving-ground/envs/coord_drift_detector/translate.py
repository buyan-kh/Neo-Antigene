from __future__ import annotations

from pathlib import Path
from typing import Dict

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


def read_fasta(path: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    name = None
    chunks = []
    with open(path, "r", encoding="utf-8") as handle:
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


def translate(nucleotides: str, stop_at_stop: bool = True) -> str:
    residues = []
    for index in range(0, len(nucleotides) - 2, 3):
        codon = nucleotides[index : index + 3]
        residue = CODON_TABLE.get(codon, "X")
        residues.append(residue)
        if residue == "*" and stop_at_stop:
            break
    return "".join(residues)


def apply_indel(cds: str, zero_based_pos: int, ref: str, alt: str) -> str:
    return cds[:zero_based_pos] + alt + cds[zero_based_pos + len(ref) :]


def downstream_peptide(cds: str, variant_cds_pos: int, ref: str, alt: str) -> tuple:
    mutant = apply_indel(cds, variant_cds_pos, ref, alt)
    first_changed_codon = variant_cds_pos // 3
    tail = mutant[first_changed_codon * 3 :]
    peptide = translate(tail)
    return peptide, first_changed_codon
