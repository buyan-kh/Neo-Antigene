# Public exams

A fair exam of a picker needs four things at once: the full list of
mutations that were considered, the person's HLA, the short peptides, and a
yes or no on each peptide from an immune-cell test. The rule has to be
locked before looking at the yeses. Several famous datasets have some of
these and lack the others. They answer different questions.

## Parkhurst 2019, the colon hits

Parkhurst and colleagues, *Cancer Discovery* 2019. PMID 31164343. PMC7138461.
Seventy-five people with microsatellite-stable gastrointestinal cancers.
124 immunogenic neoantigens, about 1.6% of the variant transcripts they
screened.

Colorectal subset: 51 people, 45 with at least one hit, 94 immunogenic
neoantigens, 47 seen by CD8 cells and 47 by CD4 cells. The free main table
names the gene and the amino-acid change, for example KRAS G12D or TP53
R175H. It does not name the short 8–11 letter peptide.

Of those 94 colorectal names, 92 are simple missense changes. A 25-letter
stretch can be rebuilt as twelve residues on each side of the change. Two
are not simple: SMAD5 (P268delinsPKH) and COPS2 (P42L, P308L).

Moderna's 64 are an unnamed subset of those 94, from 29 of the 51
colorectal patients. The paper's supplements are supposed to hold HLA
(supplementary table 4), the variants that were screened (supplementary
table 2), and NetMHCpan 4.0 minimal epitopes for the CD8 hits (supplementary
table 6). Downloads of those spreadsheets from this environment returned a
captcha page, not the tables. They are not in this repository.

What the public names can be used for: positive controls. Rebuild the
25-letter stretch for a simple missense, and check known peptides such as
the KRAS G12D CD8 peptides restricted by HLA-C*08:02. What they cannot be
used for: a score against Moderna, or a training set that is then quoted as
a win. Half the hits are CD4, so a class-I-only picker misses half of them
before it starts. The screened mutations that were silent are not in the
main table, so there is no denominator.

The NCI often tested a 25-letter peptide or a minigene, not the 8–11 letter
piece that goes in a vaccine. A "hit" in that paper is not automatically a
vaccine peptide.

## TESLA

Wells and colleagues, *Cell* 2020. DOI 10.1016/j.cell.2020.09.015.
PMC7652061. Six people (three melanoma, three lung). Twenty-five teams.
608 peptides tested, 37 immunogenic, about 6%.

The peptide list is Table S4. The matched tumor sequencing is on Synapse
(`syn21048999`) and in dbGaP (`phs000452`), behind a data-use agreement this
project does not have. Forty-one percent recall and six percent precision
are different scores. TESLA can teach the hit rate of a prospective screen.
It cannot crown a win over Moderna, because Moderna's 64 were never
published and TESLA's sequencing is not in hand.

## Ott 2017 and CU04

These two were actually run. The numbers and the limits are in
`docs/BENCHMARK.md`. Short version: both are too small to separate methods.
CU04 cannot reach a significance line at any outcome. Ott is a re-rank of a
list the authors already chose, with clonality filled in as a constant
because variant allele fraction, depth, and purity were not published.

## What a later picker would be

A universal picker means one system that, for a new person, takes that
person's mutations and HLA and proposes pieces, including mutations it has
never seen. It does not know that person's T cells. TESLA's prospective
rate was about 6% (37 of 608). A doubling of that rate would be a product,
and it would still be wrong often.

The stronger product shape is the one KEYNOTE-603 already uses as a
manufacturing format, plus a grade: the model nominates a few hundred, the
contract lab tests blood, and the vaccine gets the yeses. The first forty
patients are one cancer and the HLA types those people happen to have. A
list copied from one person does not transfer. The bits that work in one
tumor are usually absent from the next person.

"Tissue type" is not an input. HLA is. The same mutation in a different HLA
is a different peptide. The organ decides which study the person is in, not
the recipe.
