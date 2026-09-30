# Burden vs SKAT: the contrast is real, gene-dependent, and the REGENIE default is a defensible trade-off

Supersedes `power_result.md`, which generalized from a single gene and was
wrong to. Four causal genes, 3 effect sizes, 2 direction arms, 50 replicates,
M2_lof_mis mask, MAF < 0.01, alpha 0.05. Standard error about 0.07 per cell.

Burden false positive rate across the 55 non-causal genes in the same runs:
0.0516 (same arm) and 0.0490 (mixed arm) against nominal 0.05. The simulation
is calibrated; nothing below is an artifact of the setup.

## The result

Power at h2 = 0.005:

| gene | causal vars | dir | MACthr 10 burden | SKAT | MACthr 1 burden | SKAT |
| --- | --- | --- | --- | --- | --- | --- |
| ABCE1 | 5 | same | 0.92 | 0.92 | 0.92 | 0.92 |
| ABCE1 | 5 | mixed | 0.12 | 0.12 | 0.12 | 0.12 |
| ERO1A | 23 | same | 0.94 | 0.94 | 0.94 | **0.06** |
| ERO1A | 23 | mixed | 0.30 | 0.28 | 0.30 | **0.10** |
| DLGAP1 | 40 | same | 0.94 | 0.94 | 0.94 | 0.80 |
| DLGAP1 | 40 | mixed | 0.36 | 0.46 | 0.36 | **0.60** |
| FBN1 | 116 | same | 0.86 | 0.84 | 0.86 | 0.56 |
| FBN1 | 116 | mixed | 0.14 | 0.20 | 0.14 | 0.20 |

Burden is identical across the two MACthr columns in every row. That is the
internal control: the flag touches only variance-component tests, and it did.

## Four genes, four behaviours

**DLGAP1 (40 variants) is the textbook case.** Turning off collapsing gives
SKAT 0.60 against burden 0.36 under mixed effects, and costs SKAT power under
aligned effects (0.80 against burden's 0.94). Burden wins when signs agree,
SKAT wins when they do not. This is the pilot result, and it replicates at 50
replicates.

**FBN1 (116 variants) shows a weak version.** SKAT gains 0.06 in the mixed arm
and loses 0.30 in the same arm. Directionally consistent, but the mixed-arm
advantage is within noise at this replicate count.

**ERO1A (23 variants) breaks.** At MACthr 1, SKAT has power 0.06 where burden
has 0.94, and 0.10 where burden has 0.30. That is not a trade-off, it is SKAT
having essentially no power at all. Most plausible reading: 23 variants that
are nearly all ultra-rare, so without collapsing the variance component is
estimated from almost nothing. This is precisely the instability the default
exists to prevent. Stated as a hypothesis; the allele frequency spectrum of
ERO1A's causal set was not checked.

**ABCE1 (5 variants) is identical in both arms.** MACthr 1 still collapses
MAC <= 1, and a 5-variant set that is essentially all singletons stays collapsed
either way. The flag cannot help a gene this sparse.

## Corrected conclusion

`--vc-MACthr` is a bias-variance trade-off, not a mistake to switch off.

- At the default (10), SKAT behaves like burden at n = 2504, because MAC 10 is
  MAF 0.002 and about 86% of this panel's rare variants sit below it. The
  sign-sensitivity advantage is lost, and the test is stable.
- At MACthr 1, that advantage is recovered **for genes with enough
  non-singleton variants** (DLGAP1) and SKAT becomes badly underpowered for
  genes dominated by ultra-rare ones (ERO1A).

Which side you want depends on the allele frequency spectrum of the gene, and
in a panel it differs gene to gene, so no single setting is right for all of
them. That is a more useful thing to know than either textbook statement.

## What the pilot got wrong, and why

The pilot ran one gene, DLGAP1, at 25 replicates, and produced a clean
reversal. I wrote it up as the finding. Three of the four genes here do not
behave that way, and one behaves in a way that contradicts the recommendation
the pilot implied.

DLGAP1 was not a bad choice; it happens to have the variant count where the
effect is clearest. But a single gene cannot distinguish "this is how the
method behaves" from "this is how this gene behaves", and the pilot's job was
to locate the effect size grid, not to establish a result. Reading a finding
off it was the error.

## Still unresolved

- SKATO and ACAT are emitted for some genes and not others (ABCE1 and FBN1 have
  them, DLGAP1 does not). No explanation. All numbers above are plain SKAT.
- MACthr 1 chunks ran 380s against ~840s for MACthr 10, better than twice as
  fast. Test availability is the same or greater at MACthr 1, so it is not
  skipping work. Unexplained.
- Only four genes. The gene-dependence is established; its relationship to the
  allele frequency spectrum is asserted from one case and not measured.
