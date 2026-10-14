# Peptide-MHC benchmark

This directory holds a trained predictor's training records, a held-out
benchmark, and the predictor's scores on that benchmark.

The benchmark was assembled from a public epitope resource and kept separate
from the training records, so no leakage is expected. Reported accuracy on
this benchmark has been strong and the team considers it ready to publish.

Files:
  predictor_training_set/     records the predictor was built on, as received
    panel_a.csv               binding panel, peptide + allele per row
    panel_b.tsv               legacy panel export (tab separated); `sequence` is the
                              measured peptide in upper case with three residues of
                              flanking context on each side in lower case
    panel_c.csv               panel keyed by allele_id
    allele_ids.tsv            allele_id -> allele name for panel_c
    held_out_from_fit.csv     pairs that appear in the panels above but were
                              held out of the fit and were not used to build
                              the predictor
  benchmark.csv               held-out evaluation pairs
  predictions.csv             predictor score for each benchmark pair
