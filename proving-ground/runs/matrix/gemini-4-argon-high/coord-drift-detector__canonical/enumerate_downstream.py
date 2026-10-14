"""Downstream frameshift peptide enumeration.

Reduced reimplementation of the VEP Downstream plugin logic.

IMPORTANT: the coordinate convention here targets Ensembl release 113.
Under that convention the downstream tail is taken starting one codon after
the codon containing the variant.
"""

ENSEMBL_RELEASE_TARGET = 113

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


def translate(nucleotides):
    residues = []
    for index in range(0, len(nucleotides) - 2, 3):
        residue = CODON_TABLE.get(nucleotides[index:index + 3], "X")
        residues.append(residue)
        if residue == "*":
            break
    return "".join(residues)


def apply_indel(cds, zero_based_pos, ref, alt):
    return cds[:zero_based_pos] + alt + cds[zero_based_pos + len(ref):]


def downstream_peptide(cds, variant_cds_pos, ref, alt):
    """Return the frameshifted tail and where it starts.

    Release 113 convention: the codon containing the variant is
    dropped from the tail, so enumeration begins at the following codon.
    """
    mutant = apply_indel(cds, variant_cds_pos, ref, alt)
    first_codon = variant_cds_pos // 3 + 1
    tail = mutant[first_codon * 3:]
    return translate(tail), first_codon
