# Measured mappability: what the segdup track actually says

Source: UCSC hg38 `genomicSuperDups`, overlap fraction of each gene body,
computed identically for all 56 panel genes. Full table in
`config/gene_mappability.tsv`.

Eight genes exceed 50% segmental duplication overlap:

| Gene | List | Cluster | segdup frac |
| --- | --- | --- | --- |
| ZNF717 | modifier | unclustered | 1.000 |
| PLEKHM1 | modifier | 17q21.31 | 1.000 |
| ARL17B | modifier | 17q21.31 | 1.000 |
| ARL17A | modifier | 17q21.31 | 1.000 |
| LRRC37A | candidate_17q21 | 17q21.31 | 1.000 |
| LRRC37A2 | candidate_17q21 | 17q21.31 | 1.000 |
| NSF | modifier | 17q21.31 | 0.697 |
| ARHGAP27 | modifier | 17q21.31 | 0.521 |

## Two priors this falsifies

**"16p11.2 is a segmental duplication region."** Wrong for these genes. All 14
are at 0.000 except ITGAX at 0.049. The recurrent 16p11.2 CNV breakpoints sit at
roughly 29.6 to 30.2 Mb; only CDIPT is anywhere near them, and the other 13 are
distal to the BP4-BP5 interval in ordinary sequence. So 16p11.2 in this panel is
a **gene density and LD problem, not a mappability problem**. Those are different
problems with different fixes: LD correlation is handled by stratified reporting
and by taking intervals over replicates, and neither is helped by a mappability
filter.

**"17q21.31 genes have bad mappability."** Too coarse. MAPT is 0.000. CRHR1,
WNT3, SPPL2C and STH are 0.000. KANSL1 is only 0.172. The duplication is in the
low copy repeat blocks that flank and punctuate the inversion, which is where
ARL17A/B, LRRC37A/A2, PLEKHM1 and the distal part of NSF sit. MAPT lies between
the blocks, in unique sequence.

So the stratum is not "17q21.31". It is the LCR blocks inside 17q21.31, plus
ZNF717 on chr3, which was the one gene I guessed right.

## LRRC37A and LRRC37A2 are a specific problem

Both are about 45.4 kb and both are 100% duplicated. They are paralogs of each
other. Short reads from one will map to the other, so their rare variant calls
are not merely noisy, they are **correlated with each other in particular**.
Treat any agreement between those two gene tests as one observation, not two.

## Falsifiable prediction for the slicing QC step

Stated before running it, so the check is a real check:

If the segdup annotation is capturing what it should, then in the 1000G callset
the eight flagged genes should show **heterozygote excess and HWE departure**,
because collapsed paralogs produce spurious heterozygotes, while MAPT, CRHR1 and
all of 16p11.2 should look ordinary.

If that does not happen, the segdup flag is not predicting anything about this
callset and the stratification should be built from the callset QC alone.

## Small gene warning

STH is 445 bp and SPPL2C is 2187 bp. At AAF < 0.001 in 2504 samples they will
likely carry zero qualifying variants and REGENIE will emit no test. Same risk
for HSPB1 (1623 bp), HSD3B7 (3941 bp), HBB (3932 bp), PRSS8 (4295 bp) and
HIGD1B (4124 bp). Expected, not a bug, and the count of testable genes per mask
by AAF cell is recorded as a result.
