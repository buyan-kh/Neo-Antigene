#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# NCBI standard table, written independently
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
BASES = "TCAG"
TABLE = {}
i = 0
for a in BASES:
    for b in BASES:
        for c in BASES:
            TABLE[a + b + c] = AA[i]
            i += 1


def read_fasta(path):
    seqs = {}
    name, chunks = None, []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.upper())
    seqs[name] = "".join(chunks)
    return seqs


def translate(seq):
    out = []
    for i in range(0, (len(seq) // 3) * 3, 3):
        aa = TABLE[seq[i:i+3]]
        out.append(aa)
        if aa == "*":
            break
    return "".join(out)


transcripts = read_fasta(ROOT / "transcripts.fasta")

# GFF CDS in transcript coords (origin 1000000)
gff_cds = {}
for line in open(ROOT / "annotation.gff3"):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip().split("\t")
    if f[2] != "CDS":
        continue
    start = int(f[3]) - 1000000
    end = int(f[4]) - 1000000
    parent = dict(x.split("=") for x in f[8].split(";"))["Parent"].split(":")[-1]
    gff_cds[parent] = (start, end)

variants = []
for line in open(ROOT / "variants.vcf"):
    if line.startswith("#"):
        continue
    f = line.rstrip().split("\t")
    info = dict(x.split("=", 1) for x in f[7].split(";"))
    variants.append((f[2], f[3], f[4], info["TRANSCRIPT"], int(info["CDSPOS"])))

print("NCBI table TGA/TAA/TAG:", TABLE["TGA"], TABLE["TAA"], TABLE["TAG"])
print()

for vid, ref, alt, tid, cds_pos in variants:
    start, end = gff_cds[tid]
    cds = transcripts[tid][start - 1:end]
    assert cds[cds_pos - 1: cds_pos - 1 + len(ref)] == ref, (vid, cds[cds_pos-1:cds_pos-1+len(ref)], ref)
    z = cds_pos - 1
    wt_codon_i = z // 3
    mutant = cds[:z] + alt + cds[z + len(ref):]
    wt_aa = translate(cds[wt_codon_i * 3:])
    mut_from_var = translate(mutant[wt_codon_i * 3:])
    mut_from_next = translate(mutant[(wt_codon_i + 1) * 3:])
    # also translate full WT protein for first differing AA
    wt_prot = translate(cds)
    mut_prot = translate(mutant)
    first_diff = None
    for i, (a, b) in enumerate(zip(wt_prot, mut_prot)):
        if a != b:
            first_diff = i
            break
    else:
        if len(mut_prot) != len(wt_prot):
            first_diff = min(len(wt_prot), len(mut_prot))
    print(vid)
    print(f"  pos_in_codon={z%3} insdel={len(alt)-len(ref)} codon={wt_codon_i}")
    print(f"  wt_codon={cds[wt_codon_i*3:wt_codon_i*3+3]} -> {TABLE[cds[wt_codon_i*3:wt_codon_i*3+3]]}")
    print(f"  mut_codon={mutant[wt_codon_i*3:wt_codon_i*3+3]} -> {TABLE[mutant[wt_codon_i*3:wt_codon_i*3+3]]}")
    print(f"  first_diff_aa_index={first_diff} wt_aa_at_var={wt_prot[wt_codon_i] if wt_codon_i < len(wt_prot) else 'NA'} mut_aa={mut_prot[wt_codon_i] if wt_codon_i < len(mut_prot) else 'NA'}")
    print(f"  pep_from_var={mut_from_var}")
    print(f"  pep_from_next={mut_from_next}")
    print(f"  full_mut_from_diff={mut_prot[first_diff:] if first_diff is not None else mut_prot}")
