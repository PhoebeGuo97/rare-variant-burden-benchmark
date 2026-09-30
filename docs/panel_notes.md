# Panel notes: what the gene list does to the study design

The panel is 56 AD/PD candidate genes: the 51 shared modifiers from the
dissertation, plus 5 additional 17q21.31 candidates (KANSL1, LRRC37A,
LRRC37A2, SPPL2C, STH). The `list` column records which set each gene came
from. Both sets are analyzed identically.

## What is and is not a control here

Phenotypes in the calibration arm are simulated at random. **Under a null
phenotype every gene in the panel is null**, whether it is a validated shared
modifier, a single-disease gene, a gene with opposite effects in AD and PD, or a
gene that was never tested. Biology is therefore irrelevant to type I error, and
no biologically null gene is needed anywhere in this design. The control for
calibration is the phenotype simulation itself.

An earlier version of this file labeled the five 17q21.31 genes as "controls".
That was wrong on both counts: no biological control is required, and those
genes are dissertation candidates, so the label misstated the underlying work in
a public repository. They are panel genes.

The property actually worth stratifying on is **mappability**, which belongs to
the genomic region and says nothing about candidate status. It is measured, not
asserted: see `config/gene_mappability.tsv`.

## Finding 1: the panel has two clusters, and 16p11.2 is the bigger one

**16p11.2, chr16:29.6-31.5 Mb, 14 genes.** CDIPT, SETD1A, HSD3B7, STX1B, STX4,
PRSS53, VKORC1, BCKDK, KAT8, PRSS8, PRSS36, FUS, ITGAM, ITGAX. Nine fall in a
200 kb window from SETD1A to PRSS36. STX1B and HSD3B7 are 109 bp apart.

**17q21.31, chr17:45.3-46.9 Mb.** ARHGAP27, PLEKHM1, CRHR1, MAPT, ARL17B,
ARL17A, NSF, WNT3, plus the five added candidates. FMNL1 (45.22 Mb), HIGD1B
(44.85 Mb) and GOSR2 (46.92 Mb) sit just outside the boundary used here and are
treated as unclustered, which is a judgment call, not a fact.

So the effective number of independent units is around 31, not 56. The two
clusters contribute closer to one unit each than to fourteen and thirteen.

## Consequence 1: confidence intervals go over replicates, not over genes

The original plan was to pool all gene tests and put a Wilson interval on the
pooled type I error. That interval would be far too narrow: rare variants in
adjacent genes within a haplotype block are correlated, and gene level tests
inherit the correlation.

Corrected: for each gene and mask, type I error is the proportion of the N
simulated replicates with p below alpha. Replicates are independent by
construction, so that is an honest binomial. Report the **distribution** of per
gene type I error across genes, stratified by cluster and by segdup flag, rather
than one pooled number. The distribution is also the more useful result. A mask
that is calibrated on average while being anticonservative at three genes is a
mask you would reject, and pooling hides exactly that.

## Consequence 2: mappability is measured, twice

1. **A priori.** Fraction of each gene body overlapping the UCSC hg38 segmental
   duplication track (`scripts/01b_mappability.py`). Computed for all 56 genes,
   continuous, no gene singled out by assumption.
2. **Empirical, at the slicing stage.** Per variant call rate, HWE departure and
   excess heterozygosity in the actual 1000G callset. Collapsed paralogs in
   duplicated regions show up as heterozygote excess. This is the better
   measurement because it comes from the data being analyzed.

Where the two disagree, the empirical one wins. If the 17q21.31 genes do not
come out artifact prone, the prior was wrong and the stratification follows the
measurement.

ZNF717 (chr3, highly paralogous zinc finger) is a candidate for the same flag on
the same objective basis, not on my say so.

## Consequence 3: power is indexed by variant count, and NOT because of gene size

Power curves are reported against the **number of qualifying variants in the
mask**, not against gene identity.

The original justification for this was wrong and is corrected here. It said
gene bodies span 1.6 kb to 959 kb, so a per-gene curve would mostly be a plot of
gene length. Measured, the correlation between gene body length and rare coding
variant count is **r = 0.118**. Gene body length is mostly intron; coding length
is what produces coding variants. DLGAP1 is 959 kb with 73 rare coding variants,
ARAP1 is 108 kb with 163.

The conclusion stands for a different reason: qualifying variant count is the
quantity the test actually consumes, so indexing by it normalizes by information
content rather than by size. Those are different claims and only the second is
supportable here.

Practical consequence for Day 1: small genes will have zero or one qualifying
variant under the strict masks (M1_lof at AAF < 0.001 in 2504 samples) and
REGENIE will emit no test for them. The count of testable genes per mask by AAF
cell is itself a result and gets recorded.

## The panel is GWAS derived, and what that does and does not cost

These 56 genes were selected on common variant evidence. Two consequences.

It does **not** limit rare variant content. Every gene carries rare coding
variation, because in 2504 people that follows from mutation rate and sample
size rather than from how the gene was chosen: median 51 rare coding variants
per gene at MAF < 0.01, minimum 3, none at zero. The sparse category is LoF
(median 2, ten genes at zero), and that is LoF being rare everywhere, not an
artifact of GWAS selection.

It does mean the panel would be a **poor choice for a rare variant discovery
study**, and there is no reason to expect real rare signal in these genes. This
project never asks that question, since the phenotypes are simulated, but the
limitation belongs in the write-up so the panel is not misread as a discovery
attempt.

Incidental: panel median gene length is 37.8 kb against a genome-wide median
near 25 kb, the mild large-gene skew you would expect from GWAS discovery.
Given r = 0.118 above, it barely affects variant counts.

## An unused connection worth keeping

Reason (2) for excluding a gene from the shared modifier list, opposite effects
in AD and PD, is the exact scenario the burden versus SKAT contrast is about.
Burden collapses signal by summing, so opposing effects cancel; SKAT tests the
variance and does not care about sign. That is a cleaner motivation for the
mixed-direction arm of the power grid than a generic simulation, and it is worth
saying that way in the write-up.
