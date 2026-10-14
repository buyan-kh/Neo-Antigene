# Proving Ground

**We sell the exam, not the student.**

AI labs can already generate endless answers. What they cannot easily build is
a *test that cannot be cheated* on real scientific work. That is what this is:
a set of self-grading tests where a model is handed messy real files, asked to
do a piece of careful scientific checking, and then marked by plain code that
either reproduces the right answer or does not.

Every test here comes from a mistake that actually shipped in a real tool and
was caught by hand. That is the part a labeling vendor cannot copy: they can
hire people to write questions, but they do not know which mistake to look
for.

## The rules every test follows

- **Marked by code, never by a person or another model.** Same answer in, same
  mark out, every time.
- **Cannot be gamed.** The mark is based on numbers the grader works out for
  itself and on lists the model had to actually derive. Not on guessable
  words.
- **The answer sheet is out of reach.** The model sees only the input files.
  The answer key lives in a separate folder the model never gets.
- **Guessing the same thing every time fails.** Several tests ship in mirrored
  pairs, where the lazy constant answer is right in one and wrong in the
  other.

## How a test runs

1. The model gets a folder of input files and the task description
   (`task.md`). No network.
2. The model writes its answer to one file: `/output/submission.json`.
3. A separate grader process reads that answer, reads the hidden key, and
   produces a score.

The grader always returns the same shape:

```json
{
  "score": 0.8,
  "passed": false,
  "checks": [{"name": "recomputation", "passed": false, "detail": "..."}],
  "reward_hack_flags": ["fabricated a metric with 0 clean pairs"]
}
```

Some checks are gates. Miss a gate and the score is zero no matter what else
is right.

## What is built so far

| Test | The mistake it catches |
| --- | --- |
| `contam-sentinel` | The test set was secretly inside the training data, hidden because the same thing was spelled two different ways. |
| `power-gate` | A small "significant" number mistaken for a study that was actually big enough. |

Three more are specified and not yet built: `pvalue-police`, `allele-normalizer`,
and `coord-drift-detector`.

### contam-sentinel

A tool is reported to score brilliantly on a held-out test. In the canonical
version, **every single one of the 1,486 test examples is already in the
training data** — invisible to anyone who compares the text exactly as
written, because the names are spelled differently in the two files. The model
has to normalise the names, find the overlap, and then refuse to report a
clean score, because there is nothing clean left to score.

The mirrored version leaks only half. There, a clean score *can* be computed,
and it comes out at **0.518 — a coin flip.** The brilliant tool was memorising.

### power-gate

A study reports a genuinely significant result (p = 0.014, under the 0.05
bar). It is also badly underpowered: a study that size would only detect the
claimed effect about **25%** of the time, against a stated target of 80%. It
would have needed a sample of **247**, not 60.

The mirrored version is properly powered, so a model that reflexively cries
"too small" without computing gets zero.

The definition of "big enough" is pinned exactly in `task.md` (Fisher's
noncentral hypergeometric at the stated odds ratio). Without pinning it, two
careful statisticians would get two different numbers and both would fail.

## Running it

Needs Python 3.9+. No packages to install.

```bash
# regenerate the input files and answer keys
python3 envs/contam_sentinel/build.py
python3 envs/power_gate/build.py

# mark one answer
python3 envs/contam_sentinel/grade.py --submission path/to/submission.json --fixture canonical
python3 envs/power_gate/grade.py --submission path/to/submission.json --fixture canonical

# prove the marking never drifts
python3 verify_determinism.py --runs 100
```

## Layout

```
contract/              the shared rules every test obeys
  core.py              answer loading, checks, scoring, grader output
  determinism.py       runs a grader many times and compares
envs/<test>/
  task.md              what the model is asked to do
  build.py             makes the input files and the answer key
  grade.py             marks an answer
  fixtures/<name>/
    inputs/            what the model sees
    key/               the answer sheet; never shown to the model
  reference/<name>/    known-good and known-bad answers, used to test the grader
verify_determinism.py  runs every grader 100 times and checks nothing moves
```
