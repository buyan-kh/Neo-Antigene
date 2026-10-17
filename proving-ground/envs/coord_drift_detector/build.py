from __future__ import annotations

import csv
import json
import random
import sys
from pathlib import Path

ENV_DIR = Path(__file__).resolve().parent
ROOT = ENV_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from envs.coord_drift_detector.translate import downstream_peptide  # noqa: E402

SEED = 20261014
TRANSCRIPT_COUNT = 30
AFFECTED_COUNT = 4
END_DECOY_COUNT = 6
END_DECOY_NT = 6
CDS_CODONS = 180
ANNOTATION_RELEASE = 114
CODE_RELEASE = 113
DRIFT_CODONS = 1
GENOMIC_ORIGIN = 1000000
MECHANISM = "annotation_release_mismatch"
MECHANISMS = [
    "none",
    "annotation_release_mismatch",
    "zero_vs_one_based_position",
    "downstream_start_skips_variant_codon",
    "wrong_strand",
    "wrong_codon_table",
    "indel_applied_at_wrong_sequence_level",
    "translation_not_stopped_at_stop_codon",
]

BASES = "ACGT"
STOPS = ("TAA", "TAG", "TGA")

ENUMERATE_SOURCE = '''"""Downstream frameshift peptide enumeration.

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
            fields = line.rstrip("\\n").split("\\t")
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
    first_codon = zero_based // 3
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
'''

OFFSETS_HEADER = '''"""CDS span per transcript: 1-based, inclusive, transcript coordinates.

Cached so the annotation does not have to be parsed on every run.
"""

CDS_SPAN = {
'''

README_TEXT = """# Frameshift downstream enumeration

This pipeline reads frameshift variants and works out the new protein tail
each one produces.

```
variants.vcf              the frameshift variants (transcript space)
transcripts.fasta         full transcript sequences, soft-masked
annotation.gff3           UTR and CDS coordinates on the reference chromosome
transcript_origin.tsv     genomic origin of each transcript sequence
enumerate_downstream.py   the enumeration code
cds_offsets.py            cached CDS spans used by the code
reported_neoorfs.csv      the current output
```

Run with `python enumerate_downstream.py transcripts.fasta variants.vcf out.csv`.

`annotation.gff3` is the current reference release. `cds_offsets.py` was
cached from the previous release and was not regenerated when the annotation
was refreshed.
"""


def random_codon(rng: random.Random, allow_stop: bool = False) -> str:
    while True:
        codon = "".join(rng.choice(BASES) for _ in range(3))
        if allow_stop or codon not in STOPS:
            return codon


def soft_mask(rng: random.Random, sequence: str) -> str:
    masked = list(sequence)
    for _ in range(rng.randint(2, 4)):
        length = rng.randint(15, 45)
        start = rng.randint(0, len(sequence) - length)
        for index in range(start, start + length):
            masked[index] = masked[index].lower()
    return "".join(masked)


def make_transcript(rng: random.Random, upstream_atg: bool) -> dict:
    utr5_len = rng.randint(21, 90)
    utr5 = "".join(rng.choice(BASES) for _ in range(utr5_len))
    if upstream_atg:
        utr5 = utr5[:-3] + "ATG"
    else:
        while utr5[-3:] == "ATG":
            utr5 = utr5[:-3] + random_codon(rng, allow_stop=True)
    codons = ["ATG"] + [random_codon(rng) for _ in range(CDS_CODONS - 2)] + ["TAA"]
    cds = "".join(codons)
    utr3 = "".join(rng.choice(BASES) for _ in range(rng.randint(90, 180)))
    return {
        "sequence": utr5 + cds + utr3,
        "cds_start": utr5_len + 1,
        "cds_end": utr5_len + len(cds),
    }


def build_fixture(name: str) -> dict:
    rng = random.Random(SEED)
    fixture_dir = ENV_DIR / "fixtures" / name
    inputs_dir = fixture_dir / "inputs"
    key_dir = fixture_dir / "key"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    key_dir.mkdir(parents=True, exist_ok=True)

    transcript_ids = ["ENST{:011d}".format(i + 1) for i in range(TRANSCRIPT_COUNT)]
    drift_ids = set(rng.sample(transcript_ids, AFFECTED_COUNT))
    decoy_ids = set(rng.sample([item for item in transcript_ids if item not in drift_ids], END_DECOY_COUNT))

    transcripts = {}
    for transcript_id in transcript_ids:
        transcripts[transcript_id] = make_transcript(rng, transcript_id in drift_ids)

    with open(inputs_dir / "transcripts.fasta", "w", encoding="utf-8") as handle:
        for transcript_id in transcript_ids:
            sequence = soft_mask(rng, transcripts[transcript_id]["sequence"])
            handle.write(">{} transcript, soft-masked\n".format(transcript_id))
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    origin = GENOMIC_ORIGIN
    with open(inputs_dir / "annotation.gff3", "w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write("##coordinates 1-based genomic, inclusive\n")
        handle.write("##sequence-region chr1 {} {}\n".format(origin + 1, origin + 200000))
        handle.write("##transcript-origin {}\n".format(origin))
        for transcript_id in transcript_ids:
            record = transcripts[transcript_id]
            length = len(record["sequence"])
            handle.write(
                "chr1\tensembl\tfive_prime_UTR\t{start}\t{end}\t.\t+\t.\tParent=transcript:{tx}\n".format(
                    tx=transcript_id,
                    start=origin + 1,
                    end=origin + record["cds_start"] - 1,
                )
            )
            handle.write(
                "chr1\tensembl\tCDS\t{start}\t{end}\t.\t+\t0\tID=CDS:{tx};Parent=transcript:{tx}\n".format(
                    tx=transcript_id,
                    start=origin + record["cds_start"],
                    end=origin + record["cds_end"],
                )
            )
            handle.write(
                "chr1\tensembl\tthree_prime_UTR\t{start}\t{end}\t.\t+\t.\tParent=transcript:{tx}\n".format(
                    tx=transcript_id,
                    start=origin + record["cds_end"] + 1,
                    end=origin + length,
                )
            )

    with open(inputs_dir / "transcript_origin.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", "chrom", "genomic_origin"])
        for transcript_id in transcript_ids:
            writer.writerow([transcript_id, "chr1", origin])

    stale_spans = {}
    for transcript_id in transcript_ids:
        record = transcripts[transcript_id]
        start = record["cds_start"]
        end = record["cds_end"]
        if transcript_id in drift_ids:
            start -= 3 * DRIFT_CODONS
            assert record["sequence"][start - 1 : record["cds_start"] - 1] == "ATG"
        elif transcript_id in decoy_ids:
            end += END_DECOY_NT
            assert end <= len(record["sequence"])
        stale_spans[transcript_id] = (start, end)

    with open(inputs_dir / "cds_offsets.py", "w", encoding="utf-8") as handle:
        handle.write(OFFSETS_HEADER)
        for transcript_id in transcript_ids:
            start, end = stale_spans[transcript_id]
            handle.write('    "{}": ({}, {}),\n'.format(transcript_id, start, end))
        handle.write("}\n")

    variants = []
    for index, transcript_id in enumerate(transcript_ids):
        record = transcripts[transcript_id]
        cds = record["sequence"][record["cds_start"] - 1 : record["cds_end"]]
        annotated_tail = record["sequence"][record["cds_start"] - 1 :]
        placed = False
        for _ in range(4000):
            if transcript_id in drift_ids:
                codon_index = rng.randint(CDS_CODONS - 30, CDS_CODONS - 8)
            else:
                codon_index = rng.randint(20, CDS_CODONS - 60)
            cds_pos = codon_index * 3 + rng.randint(0, 2)
            if index % 2 == 0:
                ref = cds[cds_pos]
                alt = ref + rng.choice(BASES)
            else:
                ref = cds[cds_pos : cds_pos + 2]
                alt = ref[0]
            cds_pep, _ = downstream_peptide(cds, cds_pos, ref, alt)
            correct_pep, correct_off = downstream_peptide(annotated_tail, cds_pos, ref, alt)
            stale_start, stale_end = stale_spans[transcript_id]
            stale_cds = record["sequence"][stale_start - 1 : stale_end]
            stale_pep, stale_off = downstream_peptide(stale_cds, cds_pos, ref, alt)
            differs = stale_pep != correct_pep or stale_off != correct_off
            if transcript_id in drift_ids:
                escapes = correct_pep.endswith("*") and correct_pep != cds_pep
                if escapes and differs:
                    placed = True
                    break
            elif cds_pep.endswith("*") and cds_pep == correct_pep and not differs:
                placed = True
                break
        if not placed:
            raise RuntimeError("could not place a variant on " + transcript_id)
        variants.append(
            {
                "variant_id": "{}:{}".format(transcript_id, cds_pos + 1),
                "transcript": transcript_id,
                "cds_pos": cds_pos,
                "ref": ref,
                "alt": alt,
            }
        )

    with open(inputs_dir / "variants.vcf", "w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=transcript sequences\n")
        handle.write("##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description=\"Transcript\">\n")
        handle.write(
            "##INFO=<ID=CDSPOS,Number=1,Type=Integer,Description=\"1-based position within the coding sequence in annotation.gff3\">\n"
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant in variants:
            handle.write(
                "{tx}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\tTRANSCRIPT={tx};CDSPOS={pos}\n".format(
                    tx=variant["transcript"],
                    pos=variant["cds_pos"] + 1,
                    vid=variant["variant_id"],
                    ref=variant["ref"],
                    alt=variant["alt"],
                )
            )

    (inputs_dir / "enumerate_downstream.py").write_text(ENUMERATE_SOURCE, encoding="utf-8")
    (inputs_dir / "README.md").write_text(README_TEXT, encoding="utf-8")

    correct = {}
    reported = {}
    affected = []
    for variant in variants:
        record = transcripts[variant["transcript"]]
        cds = record["sequence"][record["cds_start"] - 1 : record["cds_end"]]
        assert cds[variant["cds_pos"] : variant["cds_pos"] + len(variant["ref"])] == variant["ref"]
        peptide, offset = downstream_peptide(
            record["sequence"][record["cds_start"] - 1 :],
            variant["cds_pos"],
            variant["ref"],
            variant["alt"],
        )
        correct[variant["variant_id"]] = {
            "downstream_peptide": peptide,
            "start_offset": offset,
        }
        stale_start, stale_end = stale_spans[variant["transcript"]]
        stale_cds = record["sequence"][stale_start - 1 : stale_end]
        stale_peptide, stale_offset = downstream_peptide(
            stale_cds, variant["cds_pos"], variant["ref"], variant["alt"]
        )
        reported[variant["variant_id"]] = {
            "downstream_peptide": stale_peptide,
            "start_offset": stale_offset,
        }
        if stale_peptide != peptide or stale_offset != offset:
            affected.append(variant["variant_id"])

    assert {v.split(":")[0] for v in affected} == drift_ids

    with open(inputs_dir / "reported_neoorfs.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["variant_id", "downstream_peptide", "start_offset"])
        for variant_id in sorted(reported):
            writer.writerow(
                [
                    variant_id,
                    reported[variant_id]["downstream_peptide"],
                    reported[variant_id]["start_offset"],
                ]
            )

    key = {
        "fixture": name,
        "mechanism": MECHANISM,
        "mechanism_options": MECHANISMS,
        "magnitude": {"offset": DRIFT_CODONS, "direction": "downstream"},
        "annotation_release": ANNOTATION_RELEASE,
        "code_release": CODE_RELEASE,
        "affected_variants": sorted(affected),
        "affected_transcripts": sorted({v.split(":")[0] for v in affected}),
        "corrected": correct,
        "reported": reported,
        "seed": SEED,
    }
    (key_dir / "meta.json").write_text(
        json.dumps(key, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    write_reference_submissions(name, key)
    return {
        "fixture": name,
        "variants": len(variants),
        "affected": len(affected),
        "annotation_release": ANNOTATION_RELEASE,
        "code_release": CODE_RELEASE,
    }


def write_reference_submissions(name: str, key: dict) -> None:
    reference_dir = ENV_DIR / "reference" / name
    reference_dir.mkdir(parents=True, exist_ok=True)

    def rows(variant_ids, source):
        return [
            {
                "variant_id": variant_id,
                "downstream_peptide": key[source][variant_id]["downstream_peptide"],
                "start_offset": key[source][variant_id]["start_offset"],
            }
            for variant_id in variant_ids
        ]

    def write(filename: str, payload: dict) -> None:
        (reference_dir / filename).write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    write(
        "good.json",
        {
            "drift_detected": True,
            "mechanism": MECHANISM,
            "magnitude": key["magnitude"],
            "affected_variants": key["affected_variants"],
            "corrected": rows(key["affected_variants"], "corrected"),
        },
    )
    write(
        "bad_no_drift.json",
        {
            "drift_detected": False,
            "mechanism": "none",
            "magnitude": {"offset": 0, "direction": "none"},
            "affected_variants": [],
            "corrected": [],
        },
    )
    write(
        "bad_echoed_buggy_output.json",
        {
            "drift_detected": True,
            "mechanism": MECHANISM,
            "magnitude": key["magnitude"],
            "affected_variants": key["affected_variants"],
            "corrected": rows(key["affected_variants"], "reported"),
        },
    )
    write(
        "bad_wrong_mechanism.json",
        {
            "drift_detected": True,
            "mechanism": "zero_vs_one_based_position",
            "magnitude": key["magnitude"],
            "affected_variants": key["affected_variants"],
            "corrected": rows(key["affected_variants"], "corrected"),
        },
    )
    write(
        "bad_trusted_matching_tables.json",
        {
            "drift_detected": False,
            "mechanism": "none",
            "magnitude": {"offset": 0, "direction": "none"},
            "affected_variants": [],
            "corrected": rows(sorted(key["corrected"]), "reported"),
        },
    )


def main() -> None:
    print(json.dumps(build_fixture("canonical"), indent=2))


if __name__ == "__main__":
    main()
