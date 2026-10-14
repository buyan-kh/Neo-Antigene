# Feature screen

A screening analysis over a panel of candidate features. Each feature was
measured in two groups with 8 replicates per group (a few wells
failed QC and were dropped, so some features have one fewer). The measurement
noise is well characterised, so a z-test against a known standard deviation of
1.0 is used. Significance calls are false-discovery-rate controlled at 5%
using Benjamini-Hochberg.

Layout:

```
analysis/prepare.py  builds the analysis table from the raw measurements
analysis/stats.py    scoring and the significance call
analysis/run.py      driver that prints the significant features
data/all_measured.csv   every measurement that was taken
data/features.csv       the table the screen is run on
results/reported_significant.csv   the current output
```

Rebuild the analysis table with `python analysis/prepare.py`, then run
`python analysis/run.py`.
