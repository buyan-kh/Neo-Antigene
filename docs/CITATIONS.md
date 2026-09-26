# Citations behind the ranking weights

Every scored feature and its shipped weight, with the primary source. All
metadata below was resolved against Crossref, PubMed or PMC rather than
recalled, and the "what it does **not** say" column exists because three of
the attributions commonly made for these papers — including ones I made in an
earlier draft of this repo — are not supported by the papers themselves.

Weights are priors. They were frozen before the benchmark in
[`docs/BENCHMARK.md`](BENCHMARK.md) was run, and were not adjusted afterwards.

## The weights

| feature | weight | primary source | what it supports | what it does **not** support |
| --- | --- | --- | --- | --- |
| `presentation` | 2.5 | Wells 2020 (TESLA) | Binding affinity had the smallest p-value of any presentation feature tested (p = 4×10⁻⁶); their model filters on presentation first | That presentation is *measurably* necessary. The sequential architecture is inherited from Schreiber 2011 as a modelling assumption, not a result TESLA quantified |
| `clonality` | 1.5 | McGranahan 2016 | Clonal neoantigen burden tracks overall survival in lung adenocarcinoma and enhanced sensitivity to PD-1/CTLA-4 blockade; subclonal neoantigens did the opposite | Anything about vaccine-induced responses specifically |
| `expression` | 1.2 | Wells 2020 (TESLA) | Immunogenic pMHC had significantly higher tumor abundance (p = 0.01); their optimized filter used > 33 TPM | A safe threshold. 50% (5/10) of the immunogenic peptides their filter wrongly dropped were lost to the abundance cut |
| `agretopicity` | 0.5 | Duan 2014 (definition), Ghorani 2018 (outcomes) | Duan coined the differential agretopic index; Ghorani found mean DAI correlated with overall survival in 3/5 cohorts, **only for clonal peptides** | A per-peptide immunogenicity call. Both results are cohort-level survival associations |
| `tumor_selectivity` | 0.4 | — | Mechanistic: on-target/off-tumor risk. Central-tolerance logic | Nothing empirical. This is reasoning, not a measured effect |
| `self_dissimilarity` | 0.25 | Richman 2019 | Dissimilarity to the non-mutated proteome predicted peptide immunogenicity and correlated with survival after PD-1 therapy in NSCLC, independent of predicted MHC affinity | Our implementation. Richman used BLAST + Smith-Waterman (BLOSUM62) feeding the Łuksza 2017 partition function; ours is seed-and-extend BLOSUM62 without the partition function. Related measure, not the same one |
| `wt_dissimilarity` | 0.15 | — | A cheap stand-in for the above, over the positionally matched wild-type peptide | Nothing directly. Kept below `self_dissimilarity` for that reason |
| `mutation_exposure` | 0.1 | Capietto 2020 | Mutation position is an important determinant, and determines *which* affinity metric predicts immunogenicity: absolute affinity for non-anchor, relative-to-wild-type for anchor mutations | **That TCR-facing mutations are more immunogenic than anchor mutations.** The paper states "both anchor and nonanchor mutated peptides contain cases that show CD8 responses". TESLA found mutation position not useful for filtering |
| `hydrophobicity` | **0.0** | Chowell 2015 *vs* Wells 2020 | Chowell: strong bias toward hydrophobic TCR-contact residues in immunogenic epitopes, blind-validated in vivo on 364 HIV-1 Gag peptides | The sign, in tumors. TESLA found immunogenic pMHC significantly **less** hydrophobic (p = 0.04) and hydrophobicity unhelpful for filtering. The two disagree in the setting we care about, so the weight is zero |

Two features carry no empirical citation and say so. `tumor_selectivity` is
mechanistic reasoning; `wt_dissimilarity` is an implementation convenience.
Both are reported and refittable.

## Full references

1. **Wells DK, van Buuren MM, Dang KK, et al.** (Tumor Neoantigen Selection
   Alliance). "Key Parameters of Tumor Epitope Immunogenicity Revealed Through
   a Consortium Approach Improve Neoantigen Prediction." *Cell*
   2020;183(3):818-834.e13. DOI [10.1016/j.cell.2020.09.015](https://doi.org/10.1016/j.cell.2020.09.015).
   PMID 33038342. PMC7652061.
   Of 608 top-ranked peptides tested by pMHC multimer, 37 (6%) were
   immunogenic; median 3/51 per team.

2. **McGranahan N, Furness AJS, Rosenthal R, et al.** "Clonal neoantigens
   elicit T cell immunoreactivity and sensitivity to immune checkpoint
   blockade." *Science* 2016;351(6280):1463-1469.
   DOI [10.1126/science.aaf1490](https://doi.org/10.1126/science.aaf1490).
   PMID 26940869. PMC4984254.

3. **Duan F, Duitama J, Al Seesi S, et al.** "Genomic and bioinformatic
   profiling of mutational neoepitopes reveals new rules to predict anticancer
   immunogenicity." *J Exp Med* 2014;211(11):2231-2248.
   DOI [10.1084/jem.20141308](https://doi.org/10.1084/jem.20141308).
   PMID 25245761. PMC4203949 (open access).
   Source of the term "differential agretopic index".

4. **Ghorani E, Rosenthal R, McGranahan N, et al.** "Differential binding
   affinity of mutated peptides for MHC class I is a predictor of survival in
   advanced lung cancer and melanoma." *Ann Oncol* 2018;29(1):271-279.
   DOI [10.1093/annonc/mdx687](https://doi.org/10.1093/annonc/mdx687).
   PMID 29361136. PMC5834109.

5. **Richman LP, Vonderheide RH, Rech AJ.** "Neoantigen Dissimilarity to the
   Self-Proteome Predicts Immunogenicity and Response to Immune Checkpoint
   Blockade." *Cell Systems* 2019;9(4):375-382.e4.
   DOI [10.1016/j.cels.2019.08.009](https://doi.org/10.1016/j.cels.2019.08.009).
   PMID 31606370. PMC6813910.

6. **Capietto AH, Jhunjhunwala S, Pollock SB, et al.** "Mutation position is an
   important determinant for predicting cancer neoantigens." *J Exp Med*
   2020;217(4):e20190179.
   DOI [10.1084/jem.20190179](https://doi.org/10.1084/jem.20190179).
   PMID 31940002. PMC7144530.

7. **Chowell D, Krishna S, Becker PD, et al.** "TCR contact residue
   hydrophobicity is a hallmark of immunogenic CD8+ T cell epitopes." *PNAS*
   2015;112(14):E1754-62.
   DOI [10.1073/pnas.1500973112](https://doi.org/10.1073/pnas.1500973112).
   PMID 25831525. PMC4394253.

8. **O'Donnell TJ, Rubinsteyn A, Laserson U.** "MHCflurry 2.0: Improved
   Pan-Allele Prediction of MHC Class I-Presented Peptides by Incorporating
   Antigen Processing." *Cell Systems* 2020;11(1):42-48.e7.
   DOI [10.1016/j.cels.2020.06.010](https://doi.org/10.1016/j.cels.2020.06.010).
   PMID 32711842.
   Note: PMID 33091335 has the same title but is the published erratum
   (*Cell Systems* 11(4):418-419). Cite the former.

9. **Ott PA, Hu Z, Keskin DB, et al.** "An immunogenic personal neoantigen
   vaccine for patients with melanoma." *Nature* 2017;547(7662):217-221.
   DOI [10.1038/nature22991](https://doi.org/10.1038/nature22991).
   PMID 28678778. PMC5577644.
   The benchmark case. See [`docs/BENCHMARK.md`](BENCHMARK.md).

10. **Calis JJA, Maybeno M, Greenbaum JA, et al.** "Properties of MHC class I
    presented peptides that enhance immunogenicity." *PLoS Comput Biol*
    2013;9(10):e1003266.
    DOI [10.1371/journal.pcbi.1003266](https://doi.org/10.1371/journal.pcbi.1003266).
    PMID 24204222. PMC3808449.
    Cited here only to record what it is *not* evidence for: it supports a
    positional weighting (P4-6 most important) and an association with large
    aromatic side chains, which is not the same claim as hydrophobicity.

11. **Łuksza M, Riaz N, Makarov V, et al.** "A neoantigen fitness model
    predicts tumour response to checkpoint blockade immunotherapy." *Nature*
    2017;551(7681):517-520.
    DOI [10.1038/nature24473](https://doi.org/10.1038/nature24473).
    Referenced because Richman 2019 builds on its partition function, which
    this implementation does not reproduce.

## On label provenance in IEDB

The project brief cited a preprint reporting that ~55.8% of assessable IEDB
entries were labeled by prediction rather than experiment. That checks out, and
it is worth stating precisely because the precise version is narrower than the
headline.

**Preibisch G, Tyrolski M, Kucharski P, et al.** "Resolution of recursive data
corruption to transform T-cell epitope discovery." bioRxiv, posted April 2026.
DOI [10.64898/2026.03.30.710191](https://doi.org/10.64898/2026.03.30.710191).
The abstract states verbatim: "as of January 2025, 55.8% of assessable data are
labeled by computational models rather than verified experimentally."

Three caveats that matter before repeating it:

- **Scope.** The audit covers mass-spectrometry eluted-ligand entries only
  (3,971,663 records), not T-cell or B-cell assays. "Assessable" means
  Clean + Biased (3,417,839), excluding ~554k entries that were multi-allelic
  or had insufficient metadata.
- **What "labeled by models" means.** The *peptides* are real MS observations.
  It is the **HLA allele assignment** that is predictor-derived, via
  computational deconvolution or rank-threshold filtering.
- **Reproducibility.** Roughly half the "biased" count came from manual review
  of 37 publications rather than structured metadata, and it is a non-peer
  -reviewed preprint from a company whose product is the proposed alternative.

The underlying concern is independently verifiable from IEDB's own
documentation and API, which is the part worth relying on. IEDB explicitly
tags prediction-derived restrictions: its
[curation manual](https://curationwiki.iedb.org/wiki/index.php/Data_Field_Descriptions)
defines the MHC evidence code **"MHC binding prediction"** as used "when
predictive analysis of the epitope-MHC interaction... without experimental
assay, is used to infer a restriction", and **"Inferred by motif or alleles
present"** notes outright that "MHC binding predictions will also use this
evidence code".

The sharpest finding is a curation *rule*, not a count. The
[Curation Manual 2.0](https://curationwiki.iedb.org/wiki/index.php/Curation_Manual2.0)
instructs curators to record the most **precise** restriction rather than the
most direct one: if an epitope was eluted with a class-I-specific antibody
(class-level evidence only) but "was predicted to be an HLA-A\*02:01 binder,
then this 4 digit restriction represents a higher level of precision and
should be entered as the epitope's restriction". A prediction can therefore
supersede a coarser experimental call, which means prediction-derived labels
are not confined to obviously low-confidence records.

### The evidence-code counts do not reproduce 55.8%, and that is worth saying

Live counts from `https://query-api.iedb.org` (measured 2026-09-26, with
`Prefer: count=exact`):

| Table | Evidence code | Count | Share |
| --- | --- | --- | --- |
| `tcell_search` (n=578,145) | MHC binding prediction | 67,017 | 11.6% (15.1% of non-null) |
| `mhc_search` (n=5,805,610) | Inferred by motif or alleles present | 953,634 | 16.4% (19.2% of non-null) |
| `mhc_search` | Allele/locus-specific antibody | 3,425,442 | 59.0% |
| `mhc_search` | Single allele present | 600,168 | 10.3% |

So IEDB's own provenance field puts prediction-derived restriction around
12-19%, not 55.8%. This is **not** a refutation, because the two numbers
measure different things: the preprint locates the contamination in
prediction-based deconvolution and filtering applied to immunopeptidomics
datasets *upstream* of submission, and IEDB has no field that records whether
a submitted peptide's allele came from a deconvolution tool. That provenance
would land silently in the 59.0% antibody-purification bucket, and 85.7% of
MHC-ligand rows mention mass spectrometry, so the bucket is large enough for
it to hide in.

The operational consequence: **filtering on `mhc_allele_evidence` is a lower
bound on label contamination, not a sufficient cleaning step.** Neither figure
can be audited from IEDB metadata alone, and both drift as the database grows.

**How this project responds.** Benchmark labels come from a paper's own
supplementary ELISPOT table, not from an IEDB pull, so the provenance of every
label is a specific assay in a specific figure. See
[`docs/BENCHMARK.md`](BENCHMARK.md) for the exact filtering.
