# Type I error calibration: results

Null phenotypes on real 1000G genotypes. 56 genes, 3 masks, 5 AAF bins, 6 tests.
500 quantitative replicates, 200 per binary arm, 300 structured.
Full table in `results/calibration.tsv`.

## 1. The baseline is clean

Arm A1, quantitative null with PC covariates, burden test, all 15 mask by AAF
bin cells:

- lambda between 0.979 and 1.019
- type I error at alpha 0.05 between 0.0489 and 0.0523
- worst single gene 0.082

Nothing in the pipeline is inflating on its own. Every failure below is caused
by the thing that arm was built to stress.

## 2. The bcftools/boost warning was benign, and we know rather than assume

Every REGENIE set printed `boost::math::tgamma: Result of tgamma is too large
to represent`, from inside the SKAT null distribution. Burden does not use that
code path, so comparing the two on identical null data settles it:

| test | cells | median lambda | median T1E(0.05) |
| --- | --- | --- | --- |
| ADD (burden) | 15 | 0.999 | 0.0504 |
| ADD-SKAT | 6 | 0.978 | 0.0503 |
| ADD-SKATO | 3 | 1.016 | 0.0501 |
| ADD-SKATO-ACAT | 3 | 0.989 | 0.0527 |
| ADD-ACATV | 3 | 0.989 | 0.0527 |
| ADD-ACATO | 3 | 0.989 | 0.0527 |

SKAT is as calibrated as burden. The warning does not corrupt the p-values.

## 3. Lambda is a lying diagnostic for sparse binary traits

The headline result. Arm A3, 1:99 case:control, no Firth, M1_lof mask:

| AAF bin | lambda | T1E at 0.001 | vs nominal |
| --- | --- | --- | --- |
| 0.001 | **0.056** | 0.0265 | **26.5x** |
| 0.005 | 0.071 | 0.0248 | 24.8x |
| 0.01 | 0.078 | 0.0233 | 23.3x |
| all | 0.086 | 0.0222 | 22.2x |
| singleton | **0.039** | 0.0221 | 22.1x |

Lambda near 0.05 says "strongly conservative". The tail says the false positive
rate is more than twenty times nominal. Both are true, because the failure is
bimodal: with 25 cases and a sparse LoF mask most genes have zero or one case
carrier, so the statistic is degenerate and p sits near 1, dragging the median
chi-square to almost nothing, while the handful of genes with a case carrier
produce absurd tail p-values.

**A pipeline that reported only lambda would have concluded this test was safe.**
This is the strongest single thing to be able to say in an interview about
calibration, and it is a measurement, not a citation.

It is also what crashed the first analysis script: p-values below 1e-16, under a
null where case labels were assigned at random, are not representable when you
reconstruct the tail probability in double precision.

## 4. Firth repairs the tail and leaves everything else alone

Same phenotypes, `--firth --approx --pThresh 0.01`:

| arm | mask / bin | lambda | T1E(0.05) | T1E(0.001) |
| --- | --- | --- | --- | --- |
| A3 no Firth | M1_lof / 0.001 | 0.056 | 0.0348 | 0.02651 |
| A5 Firth | M1_lof / 0.001 | 0.056 | 0.0348 | **0.00186** |
| A3 no Firth | M2_lof_mis / singleton | 0.537 | 0.0655 | 0.01518 |
| A5 Firth | M2_lof_mis / singleton | 0.537 | 0.0655 | **0.00055** |

Lambda and T1E at 0.05 are identical between the two arms, which is correct and
worth understanding rather than treating as a bug: `--pThresh 0.01` applies
Firth only where the score test p is below 0.01, so the bulk of the distribution
is untouched by construction. The correction is a tail correction.

It slightly overcorrects in places: A5 M4_all/all gives 0.00009 against a
nominal 0.001, about 11x conservative.

**Severity gradient.** At alpha 0.001, uncorrected:

- 1:9, M2_lof_mis: ~2x inflated
- 1:9, M1_lof: ~10x
- 1:99, M1_lof: ~26x

Worse with a more unbalanced design, and worse with a sparser mask. Both
directions point the same way: fewer effective carriers, worse asymptotics.

## 5. The prediction I got backwards

I predicted that because rare variants carry sharper population structure than
common ones, and the PCs are built from common variants, PC correction would be
**least** complete in the rarest bins.

The opposite is true. Arm B, most confounded tertile, burden:

| mask / bin | B1 lambda (no PCs) | B1 T1E(0.05) | B2 lambda (PCs) | B2 T1E(0.05) |
| --- | --- | --- | --- | --- |
| M4_all / all | **33.12** | 0.7266 | 1.300 | 0.0959 |
| M4_all / 0.01 | 6.65 | 0.4520 | 1.056 | 0.0580 |
| M4_all / 0.005 | 3.93 | 0.3270 | 1.026 | 0.0530 |
| M4_all / 0.001 | 1.44 | 0.0977 | 1.009 | 0.0464 |
| M4_all / singleton | **0.92** | 0.0389 | 0.987 | 0.0453 |

Without any correction at all, the singleton bin is **already calibrated**
(lambda 0.92, T1E 0.039) while the unrestricted bin is catastrophic (lambda 33,
73% of tests significant at 0.05). The same ordering holds for M2_lof_mis
(17.40 down to 0.98) and M1_lof (1.76 down to 1.03).

**Mechanism.** Confounding needs the gene's burden score to differ
systematically between populations *with precision*. Common variants are both
strongly differentiated and precisely measured, so a mask containing them
inherits ancestry directly. A singleton burden is a small, high-variance count
of private variants; whatever between-population difference exists is swamped by
sampling noise, so there is little for the confounder to act on.

**What the bins are.** `--aaf-bins 0.001,0.005,0.01` made REGENIE build
cutoffs `[0.001 0.005 0.01 1] + singletons`. The fourth cutoff of 1 is the
`all` bin: no frequency restriction at all. The annotation file was never
MAF-filtered, so `all` includes common variants, and 17.5% of the slice sits
above MAF 0.01. The bins are cumulative, not disjoint, so the 0.01 bin also
contains the 0.005, 0.001 and singleton variants.

**Two separate conclusions, and they are not the same claim.**

1. Uncorrected confounding scales monotonically with the AAF ceiling across the
   whole frequency range, not only at the common boundary. The MAF < 0.01 bin,
   containing nothing but rare variants, still reaches lambda 6.65 and a 45%
   false positive rate at alpha 0.05 with no correction. Rare variants confound
   badly on their own. A lower AAF ceiling means less bias to correct, and the
   ceiling is therefore a bias-control choice as well as a power choice.

2. PC correction was **sufficient for every rare bin** here: lambda 1.009 to
   1.082 across singleton through MAF < 0.01, with type I error 0.045 to 0.058.
   It was insufficient only for the unrestricted bin, lambda 1.300 and type I
   error 0.096. That is the operationally useful result: with PCs in the model,
   a properly capped rare variant analysis was calibrated, and the residual
   problem lives entirely in the bin that should not have been run without a cap.

The residual in the unrestricted bin is consistent with these PCs being built
from only 673 pruned SNPs inside 4.98 Mb of panel regions, which is thin for
absorbing common-variant ancestry structure.

**Dose response.** Within M2_lof_mis, AAF < 0.01, no PCs: lambda 1.375 in the
least confounded tertile, 1.958 in the middle, 4.113 in the most. With PCs:
1.042, 1.002, 1.082. The correction tracks the confounding.

## What this does not cover

- Genotyping error here is non-differential by construction, since phenotypes
  are random. Differential error, the case where cases and controls are
  sequenced differently, is the failure mode that actually breaks real studies
  and this design says nothing about it.
- Exome-wide alpha 2.5e-6 is not estimable: ~28,000 tests per cell means an
  expected count below one. Reported alphas stop at 1e-4.
- Arm B tertiles are defined post hoc from the realized R2, not pre-registered.
