# Feature screen

A screening analysis over a panel of candidate features. Each feature was
measured in two groups with 8 replicates per group (a few wells
failed QC and were dropped, so some features have one fewer). The measurement
noise is well characterised, so a z-test against a known standard deviation of
1.0 is used. Significance calls are false-discovery-rate controlled at 5%
using Benjamini-Hochberg.

Layout:

```
bin/select_candidates.py  builds the analysis table from the raw export
analysis/stats.py         scoring and the significance call
analysis/run.py           driver that prints the significant features
raw/run_2024-11-02_plate_reader.csv   the instrument export
data/features.csv                     the table the screen is run on
results/reported_significant.csv      the current output
```

Rebuild the analysis table with `python bin/select_candidates.py`, then run
`python analysis/run.py`.
