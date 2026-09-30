# Rare variant QC: the segdup flag predicts genotype quality, but not in the direction assumed

Measured on the 169,249 variants with MAF < 0.01, the set that actually enters
the masks. Full table in `results/gene_rare_qc.tsv`.

## Group medians

| | segdup high (n=8) | segdup low (n=47) |
| --- | --- | --- |
| median MQ | 48.6 | 60.0 |
| % sites MQ < 50 | 53.2 | 4.9 |
| % sites with Mendelian error | 12.0 | 2.0 |
| median VQSLOD | -2.59 | +15.52 |
| rare alleles per sample per kb | 0.051 | 0.087 |
| rare variants per kb | 23.7 | 35.5 |

Three independent metrics agree and the separation is large. MQ 60 is the GATK
ceiling, so the unflagged genes sit at perfect unique mapping while the flagged
ones fall well below. Mendelian error is six fold higher, and it is real
genotype error measured in the 1000G trios, not an annotation artifact.
Median VQSLOD is **negative** in the flagged genes, meaning VQSR's own model
scores them as more likely false than true, despite every one of these variants
having passed the PASS filter.

`MQ0` was 0 for every gene and carried no information. GATK deprecated it; drop
it from the metric set rather than reporting it as a failed criterion.

## The criterion that failed, and why it matters more than the ones that passed

The prediction was that mismapping would **inflate** the burden quantity. It
does the opposite: rare alleles per sample per kb is 0.051 in flagged genes
against 0.087 in the rest, and rare variants per kb is 23.7 against 35.5.
The duplicated regions are **under-called, not over-called**.

The mechanism is straightforward in hindsight. Ambiguous mapping drives MQ down,
GATK either declines to call or calls with low confidence, VQSR assigns a poor
score, and the PASS filter removes much of what is left. What survives is a
thinner, noisier set of genotypes, not an excess of spurious ones.

**Consequence for the calibration arm.** Genotype error that is independent of
phenotype does not create false positives. The phenotypes here are simulated at
random, so error in these genes is non-differential by construction, and it
should show up as **lost power, not inflated type I error**. That is a testable
statement and the benchmark will test it: if the flagged genes show ordinary
type I error but visibly worse power curves, that is the expected result and the
right one to report.

Where this kind of error does break a study is when it is differential, for
example cases and controls sequenced on different platforms or at different
depths. This design has no such structure, which is worth saying plainly rather
than implying the benchmark covers it.

## The 0.5 segdup threshold is too generous

Two genes above it behave like ordinary genes:

- **ARHGAP27** (segdup 0.521): MQ 60.0, MQ<50 3.9%, ME 2.8%, VQSLOD +14.91.
- **KANSL1** (0.172, already below the threshold) likewise: MQ 60.0, VQSLOD +13.63.

And one gene separates on error without separating on mapping quality:

- **ZNF717** (1.000): MQ 57.4 and only 14.5% of sites below MQ 50, but 17.7%
  Mendelian error and VQSLOD -2.13. Paralogy producing genotype error without
  a mapping quality collapse.

So the a priori annotation was directionally right and quantitatively crude.

**Revised stratum definition: median VQSLOD < 0**, which is a direct statement
about the callset rather than about the reference. That gives seven genes:

PLEKHM1, ARL17B, ARL17A, LRRC37A, LRRC37A2, NSF, ZNF717.

ARHGAP27 and KANSL1 move into the ordinary stratum. This is the flag used from
here on; `segdup_frac` stays in the table as the a priori annotation it was.

## Worst genes, for reference

| Gene | segdup | median MQ | MQ<50% | ME>0% | VQSLOD |
| --- | --- | --- | --- | --- | --- |
| LRRC37A | 1.000 | 36.4 | 85.9 | 26.2 | -7.45 |
| LRRC37A2 | 1.000 | 35.7 | 88.0 | 11.5 | -6.56 |
| ARL17A | 1.000 | 36.2 | 84.0 | 12.4 | -6.18 |
| ARL17B | 1.000 | 45.6 | 57.3 | 14.1 | -3.06 |
| ZNF717 | 1.000 | 57.4 | 14.5 | 17.7 | -2.13 |
| PLEKHM1 | 1.000 | 56.9 | 23.7 | 1.4 | -1.56 |
| NSF | 0.697 | 51.5 | 49.2 | 7.9 | -1.10 |

LRRC37A and LRRC37A2 remain the pair to treat as one observation: same size,
both fully duplicated, paralogs of each other, and the two worst on every metric.

## Genes predicted to look ordinary, and that do

MAPT (MQ 60.0, VQSLOD 15.52), CRHR1 (60.0, 15.60), WNT3 (60.0, 15.35),
FUS (60.0, 15.58), INPP5D (60.0, 15.56), DLGAP1 (60.0, 15.95). The claim that
MAPT sits in unique sequence between the duplication blocks holds in the data.

## Method note: an approach that was discarded before it produced any result

Everything above comes from `scripts/05_rare_qc.sh` and
`scripts/05b_rare_qc_summary.py`, on the MAF < 0.01 set.

An earlier attempt (`scripts/04b_qc_summary.py`) tried to answer the same
question with heterozygote excess at MAF >= 0.05. It was abandoned before any
number from it was used, and `results/gene_qc.tsv` is its stale output.

It was wrong for two compounding reasons. That set is disjoint from the
variants under test and holds only 11% of the slice. And excess heterozygosity
is detected through missing homozygotes, whose expected count is N*p^2, so the
test has no power below roughly 1/sqrt(N) and could never have reached a
MAF < 0.01 analysis set at any threshold. Mapping quality, Mendelian error and
VQSLOD are per-site and frequency independent, so they apply directly to the
variants being tested.

`scripts/04_variant_qc.sh` is still in the pipeline: `05b` uses its AC, AN and
F_MISSING columns to define the rare set. Its ExcHet and HWE columns are not
used for the stratum.
