# Hardening log

Goal: each environment hard enough that a flagship model fails it, while a
careful expert can still solve it. Only the shortcut gets closed.

Loop per environment: run gemini-3.7-flash-high; if it passes, harden and
re-run; once flash fails, escalate gpt-5.6-sol-max -> gemini-4-argon-high ->
grok-4.7-xhigh -> claude-opus-5-5-max. After every change: rebuild fixtures,
re-run `verify_determinism.py --runs 100`.

## Baseline (v0, before hardening)

| env | version | model | pass | score | failed checks | time |
| --- | --- | --- | --- | --- | --- | --- |
| contam-sentinel | v0 | gemini-3.7-flash-high | PASS | 1.0 | - | 210.0s |
| power-gate | v0 | gemini-3.7-flash-high | PASS | 1.0 | - | 170.7s |
| pvalue-police | v0 | gemini-3.7-flash-high | PASS | 1.0 | - | 130.3s |
| allele-normalizer | v0 | gemini-3.7-flash-high | PASS | 1.0 | - | 107.6s |
| coord-drift-detector | v0 | gemini-3.7-flash-high | PASS | 1.0 | - | 159.9s |
| all five | v0 | gpt-5.6-sol-max | PASS 5/5 | 1.0 | - | 61-133s each |
| all five | v0 | gemini-4-argon-high | PASS 5/5 | 1.0 | - | 67-140s each |
| all five | v0 | grok-4.7-xhigh | n/a | - | cursor-agent spawn aborted | - |
| power-gate | v0 | claude-opus-5-5-max | PASS | 1.0 | - | - |
| 4 bio envs | v0 | claude-opus-5-5-max | n/a | - | blocked by Anthropic usage policy | - |

## Changes

### v1 (all environments) — strip tells, close the "copy the team's number" shortcuts

Cross-cutting:
- Every `task.md`: removed the paragraph that told the solver what the trap
  was; replaced the filled-in example JSON (which contained the actual
  answers, e.g. `"leaked_pair_count": 1486`, `"power": 0.2545`) with a typed
  skeleton. Free-text answer fields (`mechanism`, `bug_class`) became a closed
  list with plausible distractors so the exact string is still fair.
- `run_models.py`: sandbox directories renamed from `<env>__<fixture>` to
  `job-<hash>` so the working-directory name no longer names the bug.
- `verify_determinism.py`: now also asserts `good*.json` passes and
  `bad*.json` fails. This exposed two latent reference bugs (power-gate
  `bad_hedging` was the right answer on the adequately-powered mirror;
  pvalue-police `bad_wrong_location` never failed). Both fixed.

contam-sentinel:
- canonical is now a partial leak: 594 of 1486 benchmark pairs in training
  (was 1486/1486, which made "uncomputable" a constant answer).
- 300 decoy training rows with the same peptide but a different allele, and
  200 with a one-residue-different peptide and the same allele. Joining on
  peptide alone or truncating alleles now over-counts.
- Leaked scores less bimodal (0.80/0.20 +-0.15). Clean AUROC 0.743.
- Old full-leak fixture kept as `saturated` mirror.

power-gate:
- `analysis.py` now contains a correct noncentral-hypergeometric power
  routine that the team evaluates at the *observed* odds ratio (5.64) and
  prints "achieved power = 0.851, adequately powered". `claim.md` cites it.
  True power at the pre-specified effect (OR 2.5) is 0.2545.
- observed_overlap 5 -> 7 so the post-hoc number clears 0.80.
- Removed the `notes`/`power_definition` prose from study_design.json; the
  conventions live only in task.md, phrased as definitions.

pvalue-police:
- Code now really calls a correct Benjamini-Hochberg routine, but only on the
  features that already passed `p < 0.05`. With m=30 that keeps all 30.
  README says "FDR controlled at 5% using Benjamini-Hochberg".
- bug_class answer is `correction_applied_after_prefiltering`;
  `no_multiple_testing_correction` is a distractor and fails.
- localization now required to pass; span widened to the whole faulty block
  (lines 39-46). Deviation from the pasted spec, flagged.

allele-normalizer:
- 14 of 30 shared alleles are spelled identically in both files, so the
  naive join returns 14 rows and looks healthy (was 0, a giant tell).
- Decoys: `Cw7` (C*07:01) only in A vs `C*0702` only in B; `A*0104` shared
  vs `HLA-A*01:04N` (null allele) only in B. Collapsing either over-merges
  and fails join_correctness (30/31 = 0.968 < 0.99).
- Added lowercase and trailing-space spellings. Grader strips mapping keys.

coord-drift-detector:
- The drift is now in the data, not the prose. FASTA holds full transcripts,
  `annotation.gff3` (release 114) gives CDS spans, and the code reads CDS
  spans from a cached `cds_offsets.py` generated from release 113. For 7 of
  12 transcripts release 114 moved the CDS start one codon downstream (113
  had picked an upstream in-frame ATG). The other 5 are unaffected.
- Removed the old `+1` bug; the team's code now reproduces its own output
  exactly from the inputs.
- `affected_set` now required to pass (was only detection + peptides, so
  "fix everything" would have passed). Deviation from spec, flagged.
- mechanism answer is `annotation_release_mismatch`;
  `zero_vs_one_based_position` is the main distractor.

Determinism after v1: 32 submissions x 100 runs, all identical, verdicts OK.

### v2 — flash still passed all five on v1; add realistic complexity where careless-but-plausible methods diverge from careful ones

power-gate:
- Our own key was wrong. Exact-test power is a sawtooth in n (drops
  0.827 -> 0.736 at n=232 when the critical value steps 10 -> 11). Flash's
  answer 225 was right; our bisection gave 247. `stats.py` now scans and
  caches; `min_sample_size` is now required to pass. New reference
  `bad_bisected_min_n.json` (247) fails.

pvalue-police:
- Data regenerated (seed 20261023, effect 1.4, 6 features with one missing
  replicate) so that keeping the project's z-test gives BH = 8 features but
  "upgrading" to a Welch t-test gives 6 (Bonferroni 4 vs 1). The task says
  keep the project's statistic; the data now enforces it.
- Reported list is 34; still equals the uncorrected set.

allele-normalizer:
- Added three-field spellings (`HLA-B*27:05:01`) that must truncate to two
  fields, class II genes (DRB1, DQB1) with 4-character gene names, and
  serological aliases DR4/DR15/DQ8.

coord-drift-detector:
- Two more transcripts now differ in the cached table at the CDS *end* only
  (+6 nt). Their results are unchanged, so "flag every transcript whose table
  row differs" gives 9 affected instead of 7 and fails. Variants are
  resampled until every frameshift tail terminates at a stop inside the CDS,
  which removes the "translate into the UTR or not" ambiguity.
- New reference `bad_flagged_every_table_difference.json` fails.

contam-sentinel:
- Training set is now a directory of three panels as received: `panel_a.csv`
  (peptide, mhc_allele), `panel_b.tsv` (hla, sequence, ic50_nm), and
  `panel_c.csv` keyed by `allele_id` through `allele_ids.tsv` (compact
  names). Leaks and decoys are spread across all three.
- `predictions.csv` spells alleles compactly (`A*2402`) while `benchmark.csv`
  uses full names, so the clean AUROC also needs a normalised join.

Determinism after v2: 35 submissions x 100 runs, all identical, verdicts OK.

### v3 — flash still 5/5 on v2; add silent-failure traps (plausible method, confident wrong answer, no error)

power-gate:
- Raw `population.tsv` (5000 items, in_group flag) and `sample.tsv` (60
  items) added. `study_design.json` still says group size 120 (as
  pre-registered); the actual data has 115. README states which is which.
  Trusting the summary gives power 0.2545 / min_n 225; truth is 0.4181 / 235
  (the critical value steps). Bisection on the true K gives 258. All three
  are required checks. Grader now counts from the TSVs.
- New reference `bad_trusted_design_summary.json` fails.

pvalue-police:
- Pre-filter moved out of `stats.py` into `run.py` as
  `drop_flat_features()` ("skip features with no visible difference so the
  FDR step stays fast"). `stats.py` now reads as a clean BH-over-everything
  module. Bug span is `analysis/run.py` lines 22-34.

coord-drift-detector:
- FASTA is soft-masked (random lower-case runs, header says so). A parser
  that forgets to upper-case gets `X` residues.

allele-normalizer:
- Two exact duplicate rows in source_a (same string, same payload).
  `source_a_rows` in the key is distinct alleles. Counting rows gives 32.
- Four-field spellings (`HLA-A*02:01:01:01`).

contam-sentinel:
- Four-field allele spellings in training (`HLA-A*02:01:01:01`), 75 rows in
  panel_a alone.
- `panel_b.tsv` sequences carry three lower-case flanking residues on each
  side (`qkhYQPAILLRHkip`); README documents the convention.

Determinism after v3: 37 submissions x 100 runs, all identical, verdicts OK.

## Run log

| env | version | model | pass | score | failed checks | flags | time |
| --- | --- | --- | --- | --- | --- | --- | --- |
| contam-sentinel | v1 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 272.3s |
| power-gate | v1 | gemini-3.7-flash-high | PASS | 0.8 | min_sample_size (225; our key said 247 and was wrong) | - | 97.1s |
| pvalue-police | v1 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 116.8s |
| allele-normalizer | v1 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 127.4s |
| coord-drift-detector | v1 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 153.9s |
| contam-sentinel | v2 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 290.0s |
| power-gate | v2 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 83.9s |
| pvalue-police | v2 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 94.8s |
| allele-normalizer | v2 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 187.3s |
| coord-drift-detector | v2 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 134.9s |
| power-gate | v3 | gemini-3.7-flash-high | FAIL | 0.6 | power_value (0.2545 vs 0.4181), min_sample_size (225 vs 235) | - | 101.7s |
| pvalue-police | v3 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 143.2s |
| coord-drift-detector | v3 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 200.8s |
| allele-normalizer | v3 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 197.7s |
| contam-sentinel | v3 | gemini-3.7-flash-high | PASS | 1.0 | - | - | 279.7s |
| power-gate | v3 | gpt-5.6-sol-max | **FAIL** | 0.6 | power_value (0.2545 vs 0.4181), min_sample_size (225 vs 235) — trusted study_design.json, never recounted population.tsv | - | 72.9s |
| power-gate | v3 | gemini-4-argon-high | PASS | 1.0 | - | - | 685.5s |
| power-gate | v3 | grok-4.7-xhigh | PASS | 1.0 | - | - | 454.1s |
| power-gate | v3 | claude-opus-5-5-max | PASS | 1.0 | - | - | 199.0s |

**power-gate: flagship broken at v3 (gpt-5.6-sol-max). What broke it: the design summary (`study_design.json`, group size 120) disagrees with the raw population list (115). GPT-5.6 computed everything correctly from the summary and never reconciled it with the data. Hardening of power-gate stops here.**

### v4 — four remaining environments; move the trap out of the obvious file

pvalue-police:
- Analysis code is now a clean BH-over-everything pair (`stats.py` + `run.py`).
  The pre-filter lives in `analysis/prepare.py`, which writes `data/features.csv`
  (34 rows) from `data/all_measured.csv` (400). README describes the pipeline.
  Bug span is `prepare.py` lines 33-47. Grader recomputes from `all_measured.csv`.
  A solver that only reads `analysis/stats.py` sees a correct BH routine.

coord-drift-detector:
- Stripped every release number (FASTA headers, GFF, VCF, `cds_offsets.py`
  docstring). README says the sequences were refreshed and the cached table
  was not. The mismatch is only visible by comparing the two tables.

allele-normalizer:
- `HLA-A*02:01G` (G-group) only in source_b, while `HLA-A*02:01` is shared.
  Stripping the G letter over-merges. Canonical form keeps the suffix.

contam-sentinel:
- `held_out_from_fit.csv` lists 80 pairs that appear in the training panels
  but were not used to fit the predictor. True leak is 514, not 594. Naive
  overlap AUROC on the remaining 892 is 0.743; honoring the hold-out gives
  0.780 on 972. Jaccard 514/594 = 0.865 fails the 0.98 threshold.
- New reference `bad_ignored_holdout.json` fails.

Determinism after v4: 39 submissions x 100 runs, all identical, verdicts OK.

Note: gemini-3.7-flash-high hit the Ultra usage cap mid-v4 (`usage limits reset 11/9/2026`). Cheap-probe stand-in from here is `grok-4.7-low`.

| pvalue-police | v4 | grok-4.7-low | PASS | 1.0 | - | - | 29.7s |
| coord-drift-detector | v4 | grok-4.7-low | PASS | 1.0 | - | - | 88.3s |
| allele-normalizer | v4 | grok-4.7-low | PASS | 1.0 | - | - | 133.8s |
| contam-sentinel | v4 | grok-4.7-low | PASS | 1.0 | - | - | 181.7s |
| 4 remaining | v4 | gpt-5.6-sol-max | n/a | - | Ultra usage cap (same as Gemini; resets 11/9/2026) | - | 4s |
| pvalue-police | v4 | grok-4.7-xhigh | PASS | 1.0 | - | - | 195.1s |
| coord-drift-detector | v4 | grok-4.7-xhigh | PASS | 1.0 | - | - | 487.4s |
| allele-normalizer | v4 | grok-4.7-xhigh | PASS | 1.0 | - | - | 666.7s |
| contam-sentinel | v4 | grok-4.7-xhigh | PASS | 1.0 | - | - | 395.4s |

Cursor-only window: cheap probe `composer-2.5-fast`, flagship `cursor-grok-4.6-xhigh`.

| pvalue-police | v4 | composer-2.5-fast | PASS | 1.0 | - | - | 25.6s |
| coord-drift-detector | v4 | composer-2.5-fast | PASS | 1.0 | - | - | 30.5s |
| allele-normalizer | v4 | composer-2.5-fast | PASS | 1.0 | - | - | 65.4s |
| contam-sentinel | v4 | composer-2.5-fast | FAIL | 0.6 | leaked_set_match (594 vs 514), recomputation (0.743/892 vs 0.780/972) | 1 malformed leaked_pairs | 68.2s |
| pvalue-police | v4 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 116.2s |
| coord-drift-detector | v4 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 221.6s |
| allele-normalizer | v4 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 309.0s |
| contam-sentinel | v4 | cursor-grok-4.6-xhigh | **FAIL** | 0.6 | leaked_set_match (594 vs 514), recomputation (0.743/892 vs 0.780/972) — joined every panel overlap, ignored `held_out_from_fit.csv` | - | 437.1s |

**contam-sentinel: flagship broken at v4 (cursor-grok-4.6-xhigh). What broke it: 80 overlapping pairs were in the training panels but listed in `held_out_from_fit.csv` as not used in the fit. Both Composer and Cursor Grok 4.6 Extra High counted all 594 overlaps. Hardening of contam-sentinel stops here.**

### v5 — remaining three; Cursor-only models

pvalue-police:
- Raw table moved to `raw/run_2024-11-02_plate_reader.csv`. Filter script
  moved to `bin/select_candidates.py`. `analysis/` now contains only a clean
  BH driver. Bug span is `bin/select_candidates.py` 34-48.

coord-drift-detector:
- `annotation.gff3` is genomic (origin 1_000_000). `transcript_origin.tsv`
  gives the lift. Cache stays in transcript coordinates. Using GFF numbers
  as transcript coords slices off the end of the FASTA.

allele-normalizer:
- `HLA-Cw*dd:dd` spelling (historical Cw gene name). Canonical gene is C.

Determinism after v5: 39 submissions x 100 runs, all identical, verdicts OK.

| pvalue-police | v5 | composer-2.5-fast | PASS | 1.0 | - | - | 21.9s |
| coord-drift-detector | v5 | composer-2.5-fast | PASS | 1.0 | - | - | 37.8s |
| allele-normalizer | v5 | composer-2.5-fast | PASS | 1.0 | - | - | 39.9s |
| pvalue-police | v5 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 93.6s |
| coord-drift-detector | v5 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 232.6s |
| allele-normalizer | v5 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 244.5s |

### v6 — remaining three; close the "compare the two tables / BH everything" shortcuts

pvalue-police:
- Raw export now has 80 spike-in controls (`spike_####`, role=spike_in), 8
  with real effects. Team BH is run on p<0.05 across screen+spike-ins
  (reported 41). Correct BH is screen analytes only (4). BH on all 480
  rows is a different set. Task says spike-ins are not hypotheses.

coord-drift-detector:
- Cache now matches the GFF exactly. The bug is only in
  `enumerate_downstream.py`: `first_codon = pos//3 + 1` skips the variant
  codon. All 12 reported peptides are wrong. Comparing the two tables says
  "no drift" and fails. Rerunning the team's code echoes the bug.
  Mechanism is `downstream_start_skips_variant_codon`.

allele-normalizer:
- `DQA1*01:02-DQB1*05:01` heterodimer only in source_b. Splitting on `-`
  and joining to shared `DQB1*05:01` over-merges.

Determinism after v6: 38 submissions x 100 runs, all identical, verdicts OK.

| pvalue-police | v6 | composer-2.5-fast | PASS | 1.0 | - | - | 30.9s |
| coord-drift-detector | v6 | composer-2.5-fast | PASS | 1.0 | - | - | 30.5s |
| allele-normalizer | v6 | composer-2.5-fast | PASS | 0.8 | mapping_validity (1 incorrect); not required then | - | 46.0s |
| pvalue-police | v6 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 149.0s |
| coord-drift-detector | v6 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 167.8s |
| allele-normalizer | v6 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 341.3s |

### v7 — hide spike-in identity; mapping now required

pvalue-police:
- Spike-ins renamed `feat_0400`–`feat_0479` so they look like screen
  analytes. Role lives only in `raw/sample_sheet.csv`.
- Grader BH/Bonferroni uses only `role=screen`.

allele-normalizer:
- `mapping_validity` added to `required_to_pass` (Composer v6 had a wrong
  mapping and still passed).
- G-group only-b spelling is now `HLA-A*02:01:01:01G`.

Determinism after v7: 40 submissions x 100 runs, all identical, verdicts OK.

| pvalue-police | v7 | composer-2.5-fast | PASS | 1.0 | - | - | 26.6s |
| coord-drift-detector | v7 | composer-2.5-fast | PASS | 1.0 | - | - | 44.4s |
| allele-normalizer | v7 | composer-2.5-fast | FAIL | 0.8 | mapping_validity | - | 59.0s |
| pvalue-police | v7 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 128.1s |
| coord-drift-detector | v7 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 126.3s |
| allele-normalizer | v7 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 346.4s |

### v8 — strip remaining tells; alias table can be wrong

allele-normalizer:
- Canonicalize now parses a two-field name before consulting the alias
  table. The alias file maps `A*0206` / `HLA-A*02:06` → `HLA-A*02:01`
  (stale). Trusting the table over-merges A*02:06 into A*02:01.
- `HLA-A*02:06` is forced into the shared set.

coord-drift-detector:
- Removed the sentence that said the peptide starts at the variant codon.
  Softened the mechanism blurb so it does not name the +1 direction.

pvalue-police:
- Removed the sentence that said spike-ins are not hypotheses.

Determinism after v8: 41 submissions x 100 runs, all identical, verdicts OK.

| pvalue-police | v8 | composer-2.5-fast | FAIL | 0.8 | corrected_set | - | 25.7s |
| coord-drift-detector | v8 | composer-2.5-fast | PASS | 1.0 | - | - | 23.8s |
| allele-normalizer | v8 | composer-2.5-fast | FAIL | 0.0 | join_correctness, payload_integrity, mapping_validity | - | 65.5s |
| pvalue-police | v8 | cursor-grok-4.6-xhigh | FAIL | 0.8 | corrected_set | - | 114.7s |
| coord-drift-detector | v8 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 377.7s |
| allele-normalizer | v8 | cursor-grok-4.6-xhigh | PASS | 1.0 | - | - | 479.8s |
| pvalue-police | v8 | grok-4.7-xhigh | **FAIL** | 0.8 | corrected_set (8 features including spike-ins feat_0405/0433/0451/0472; screen-only BH is 4) | - | 231.2s |
| coord-drift-detector | v8 | grok-4.7-xhigh | PASS | 1.0 | - | - | 782.6s |
| allele-normalizer | v8 | grok-4.7-xhigh | PASS | 1.0 | - | - | 659.8s |

**pvalue-police: flagship broken at v8 (grok-4.7-xhigh). What broke it: BH over every feat_* in the plate-reader file, including spike-ins listed as role=spike_in in sample_sheet.csv. Hardening of pvalue-police stops here.**

### v9 — last two environments

coord-drift-detector:
- Shipped code no longer has a visible off-by-one. It matches its own output.
- 7 of 12 CDS features have phase 1 in `annotation.gff3`. The code ignores
  phase. Mechanism is `cds_phase_ignored`.

allele-normalizer:
- README says `payload_b` is a percentage of `payload_a`'s scale and must be
  divided by 100. Raw values fail `payload_integrity`.

Determinism after v9: 42 submissions x 100 runs, all identical, verdicts OK.

| coord-drift-detector | v9 | composer-2.5-fast | FAIL | 0.6 | corrected_peptides, offsets | - | 62.6s |
| allele-normalizer | v9 | composer-2.5-fast | FAIL | 0.0 | join_correctness, payload_integrity, counts, mapping_validity, no_silent_drops | - | 82.4s |
| coord-drift-detector | v9 | grok-4.7-xhigh | **FAIL** | 0.0 | detection (drift_detected=false, mechanism=none) — code reproduces its own output, phase column never applied | - | 2173.2s |
| allele-normalizer | v9 | grok-4.7-xhigh | **FAIL** | 0.8 | mapping_validity — heterodimer mapped to HLA-DQA1*01:02-DQB1*05:01; expected HLA-DQA1*01:02-HLA-DQB1*05:01. Join and the /100 scale were correct. | - | 1138.5s |

**coord-drift-detector: flagship broken at v9 (grok-4.7-xhigh). What broke it: the enumeration code matches its reported output, so the model reported no drift. Seven CDS features have phase 1. Hardening stops here.**

**allele-normalizer: flagship broken at v9 (grok-4.7-xhigh). What broke it: the DQA1-DQB1 heterodimer was canonicalized without repeating HLA- on the second gene. Hardening stops here.**
