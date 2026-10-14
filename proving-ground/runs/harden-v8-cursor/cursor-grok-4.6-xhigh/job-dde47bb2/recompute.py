#!/usr/bin/env python3
"""Independently recompute frameshift downstream peptides from sequences + annotation."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent

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
                chunks.append(line)
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def read_gff_cds(path, origin):
    spans = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "CDS":
                continue
            start = int(fields[3]) - origin
            end = int(fields[4]) - origin
            attrs = dict(item.split("=", 1) for item in fields[8].split(";"))
            parent = attrs["Parent"].split(":")[-1]
            spans[parent] = (start, end, fields[6])
    return spans


def translate(nucleotides, stop=True):
    residues = []
    for index in range(0, len(nucleotides) - 2, 3):
        residue = CODON_TABLE.get(nucleotides[index:index + 3], "X")
        residues.append(residue)
        if stop and residue == "*":
            break
    return "".join(residues)


def rc(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def main():
    origin_path = ROOT / "transcript_origin.tsv"
    origins = {}
    with open(origin_path) as handle:
        next(handle)
        for line in handle:
            tid, chrom, orig = line.strip().split("\t")
            origins[tid] = int(orig)

    origin = origins["ENST00000000001"]
    transcripts_raw = read_fasta(ROOT / "transcripts.fasta")
    transcripts = {k: v.upper() for k, v in transcripts_raw.items()}
    spans = read_gff_cds(ROOT / "annotation.gff3", origin)

    print("=== CDS spans from GFF vs transcript lengths ===")
    from cds_offsets import CDS_SPAN
    for tid, (start, end, strand) in spans.items():
        cached = CDS_SPAN[tid]
        seq = transcripts[tid]
        print(
            f"{tid}: gff_tx=({start},{end}) strand={strand} cached={cached} "
            f"tx_len={len(seq)} cds_len={end-start+1} match={cached==(start,end)}"
        )

    variants = []
    with open(ROOT / "variants.vcf") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info = dict(item.split("=", 1) for item in fields[7].split(";"))
            variants.append({
                "variant_id": fields[2],
                "chrom": fields[0],
                "pos": int(fields[1]),
                "ref": fields[3],
                "alt": fields[4],
                "transcript": info["TRANSCRIPT"],
                "cds_pos": int(info["CDSPOS"]),
            })

    reported = {}
    with open(ROOT / "reported_neoorfs.csv") as handle:
        next(handle)
        for line in handle:
            vid, pep, off = line.strip().split(",")
            reported[vid] = (pep, int(off))

    print("\n=== REF allele checks ===")
    for v in variants:
        tid = v["transcript"]
        start, end, strand = spans[tid]
        seq = transcripts[tid]
        cds = seq[start - 1:end]
        cds_pos = v["cds_pos"]
        tx_pos = v["pos"]
        ref = v["ref"]
        print(f"\n{v['variant_id']} REF={ref} ALT={v['alt']} CDSPOS={cds_pos} VCF_POS={tx_pos}")
        cds_slice = cds[cds_pos-1:cds_pos-1+len(ref)]
        tx_slice = seq[tx_pos-1:tx_pos-1+len(ref)]
        print(f"  CDS len={len(cds)} CDS at {cds_pos-1}={cds_slice}")
        print(f"  TX at {tx_pos-1}={tx_slice}")
        # CDS-relative vs transcript-relative
        # If POS is transcript position, CDS position would be POS - (start-1)
        inferred_cds_from_tx = tx_pos - (start - 1)
        print(f"  if VCF POS is transcript-space, implied CDSPOS={inferred_cds_from_tx}")
        print(f"  TX at CDS start+CDSPOS-1: tx_idx={start-1+cds_pos-1} base={seq[start-1+cds_pos-1:start-1+cds_pos-1+len(ref)]}")
        # codon context
        z = cds_pos - 1
        codon_idx = z // 3
        codon_start = codon_idx * 3
        print(f"  codon 0-based={codon_idx} nt-in-codon={z%3} wt_codon={cds[codon_start:codon_start+3]}")

    print("\n=== Independent translations (include variant codon) ===")
    for v in variants:
        tid = v["transcript"]
        start, end, strand = spans[tid]
        seq = transcripts[tid]
        cds = seq[start - 1:end]
        z = v["cds_pos"] - 1
        mutant = cds[:z] + v["alt"] + cds[z + len(v["ref"]):]
        codon0 = z // 3
        # include variant codon
        tail_incl = mutant[codon0 * 3:]
        pep_incl = translate(tail_incl)
        # skip variant codon (pipeline)
        tail_skip = mutant[(codon0 + 1) * 3:]
        pep_skip = translate(tail_skip)
        # 0-based pos without -1
        z0 = v["cds_pos"]
        if z0 < len(cds):
            mutant0 = cds[:z0] + v["alt"] + cds[z0 + len(v["ref"]):]
            codon0b = z0 // 3
            pep_0based = translate(mutant0[codon0b * 3:])
        else:
            pep_0based = "OOB"
        # apply to full transcript then extract CDS-like region
        seq_mut = seq[: v["pos"] - 1] + v["alt"] + seq[v["pos"] - 1 + len(v["ref"]):]
        # if CDS coordinates still used as slice on mutant transcript
        cds_wrong_level = seq_mut[start - 1:end]
        pep_wrong_level = translate(cds_wrong_level[codon0 * 3:])
        # apply indel to transcript, then take from CDS start to end of mutant seq (include 3'UTR)
        cds_plus_utr = seq_mut[start - 1:]
        pep_tx_from_cds = translate(cds_plus_utr[codon0 * 3:])
        # skip + include 3' UTR because pipeline only uses CDS slice without 3'UTR
        # wait pipeline uses CDS_SPAN end so 3'UTR is excluded... unless stop is after CDS
        pep_incl_utr = translate((mutant + seq[end:])[codon0 * 3:])

        rep_pep, rep_off = reported[v["variant_id"]]
        print(f"{v['variant_id']}")
        print(f"  reported: offset={rep_off} pep={rep_pep}")
        print(f"  incl_var codon={codon0} pep={pep_incl} MATCH_PEP={pep_incl==rep_pep} MATCH_OFF={rep_off==codon0}")
        print(f"  skip_var codon={codon0+1} pep={pep_skip} MATCH_PEP={pep_skip==rep_pep} MATCH_OFF={rep_off==codon0+1}")
        print(f"  0based_pos pep={pep_0based}")
        print(f"  wrong_level pep={pep_wrong_level}")
        print(f"  tx_from_cds pep={pep_tx_from_cds}")
        print(f"  incl_3utr pep={pep_incl_utr}")
        print(f"  mutant_len={len(mutant)} wt_len={len(cds)} delta={len(mutant)-len(cds)}")
        print(f"  first mutant codon at {codon0}: {mutant[codon0*3:codon0*3+3]} -> {translate(mutant[codon0*3:codon0*3+3], stop=False)}")
        print(f"  next codon: {mutant[(codon0+1)*3:(codon0+1)*3+3]} -> {translate(mutant[(codon0+1)*3:(codon0+1)*3+3], stop=False)}")


if __name__ == "__main__":
    main()
