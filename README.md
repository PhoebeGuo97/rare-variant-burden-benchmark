# Rare variant burden pipeline: calibration and power benchmarking

A gene-based rare variant association pipeline (REGENIE burden, SKAT, SKAT-O,
ACAT-V, ACAT-O) built on 1000 Genomes 30x WGS, benchmarked for type I error
calibration and power across mask definitions, MAF thresholds and effect
direction.

**This is a methods benchmark on simulated phenotypes. It is not a rare variant
discovery study and it identifies no genes.** Real genotypes, real LD, real
population structure, simulated outcomes.

---

## Headline findings

### 1. Genomic inflation factor is a misleading diagnostic for sparse binary traits

Uncorrected score test, 1:99 case:control, LoF mask, **null phenotypes with
randomly assigned case labels**:

| AAF bin | lambda | type I error at alpha 0.001 | vs nominal |
| --- | --- | --- | --- |
| 0.001 | **0.056** | 0.0265 | **26.5x** |
| singleton | **0.039** | 0.0221 | 22.1x |

Lambda near 0.05 reads as strongly conservative. The false positive rate is
more than twenty times nominal. Both are true: with 25 cases and a sparse mask
most genes have no case carrier, so the statistic is degenerate and p sits near
1, dragging the median chi-square toward zero, while the few genes with a
carrier throw extreme tail p-values.

A pipeline reporting only lambda would have called this test safe.

Firth correction (`--firth --approx --pThresh 0.01`) takes that 0.0265 to
0.00186 and leaves lambda and the alpha 0.05 rate untouched, which is correct
rather than suspicious: `--pThresh` makes it a tail correction by construction.

### 2. Confounding scales with the AAF ceiling, and PCs were sufficient for the rare bins

Structured null, most confounded tertile, M4_all mask, burden:

| bin (MAF <=) | lambda, no PCs | lambda, with PCs |
| --- | --- | --- |
| singleton | 0.92 | 0.99 |
| 0.001 | 1.44 | 1.01 |
| 0.005 | 3.93 | 1.03 |
| 0.01 | 6.65 | 1.06 |
| uncapped | **33.12** | 1.30 |

Two separate conclusions. Uncorrected confounding rises monotonically with the
AAF ceiling across the whole frequency range, so the ceiling is a bias-control
choice and not only a power choice. And PC correction was sufficient for every
rare bin (lambda 1.01 to 1.06) and insufficient only for the uncapped bin.

### 3. The burden vs SKAT distinction is erased by a REGENIE default at this sample size

DLGAP1, 40 causal variants, h2 = 0.004, alpha 0.05, effects scaled to equal
variance explained in both arms:

| | burden | SKAT |
| --- | --- | --- |
| same-direction, `--vc-MACthr 10` (default) | 0.92 | 0.84 |
| mixed-direction, `--vc-MACthr 10` | 0.24 | 0.28 |
| same-direction, `--vc-MACthr 1` | 0.92 | **0.64** |
| mixed-direction, `--vc-MACthr 1` | 0.24 | **0.52** |

With collapsing off, the textbook trade-off appears and reverses: burden wins
1.4x under aligned effects, SKAT wins 2.2x under opposing ones. With the
default, SKAT simply tracks burden.

REGENIE collapses variants with MAC <= 10 into a single burden pseudo-variant
before any variance-component test. At n = 2504 that is MAF <= 0.002, about 86%
of this panel's rare variants, and the collapse destroys the per-variant sign
information SKAT exists to use. For M1_lof, burden and all five VC tests
returned byte-identical statistics.

This is not a bug. Collapsing ultra-rare variants is a deliberate stability
choice, and in a 500,000-sample biobank MAC 10 is MAF 0.00001 so it is nearly
inert. It is a default whose effect depends entirely on cohort size.

### 4. Poor genotype quality inflates the LoF mask specifically

Median LoF share of rare coding variants is 6.27% in the seven genes with
median VQSLOD < 0, against 3.25% elsewhere. ZNF717 carries 32 LoF variants,
16.1% of its coding variants, more than ten times MAPT or FBN1 despite being
smaller than both, and has the panel's worst Mendelian error rate at 17.7%.

Under M1_lof, ZNF717 is the most powered gene in the panel and its power is
mismapped reads. Under M4_all it is unremarkable.

HBB is the control that makes this trustworthy: highest LoF fraction in the
panel at 35.7% with perfect quality metrics. That is real biology, and a filter
keyed on LoF enrichment alone would have discarded it.

---

## Data and design

1000 Genomes 30x NYGC, GRCh38, raw GT callset (not the phased panel, which
drops the low frequency variants this project is about). 2504 unrelated samples
across 5 superpopulations. 56 AD/PD candidate genes, 4.98 Mb of regions across
15 chromosomes, 205,205 variants of which 169,249 are MAF < 0.01.

Masks: M1_lof, M2_lof_mis, M4_all (adds synonymous as an annotation-level
negative control). AAF bins 0.001, 0.005, 0.01, plus singleton and uncapped.
Annotation by `bcftools csq`, protein_coding transcripts only, most severe
consequence per gene.

Simulation arms: unstructured quantitative null (500 reps), binary at 1:9 and
1:99 with and without Firth (200 each), structured null with and without PCs
(300), and a power grid over effect size, direction mix and causal gene.

---

## Stated limitations

- **Simulated phenotypes.** No claim about any of these genes.
- **n = 2504 is small.** Power curves describe relative ordering of tests and
  masks, not absolute detectability at biobank scale.
- **M1_lof is not really testable on a 56-gene panel.** Only 2 genes carry 10 or
  more LoF variants at MAF < 0.01. M1 results are reported with variant counts
  attached and nothing is concluded from M1 alone.
- **Exome-wide alpha 2.5e-6 is not estimable**, with about 28,000 tests per cell.
  Reported alphas stop at 1e-4.
- **PCs come from 673 LD-pruned SNPs inside the panel regions**, not a
  genome-wide scan, excluding 17q21.31 and the low-quality genes.
- **The panel is GWAS derived**, so it would be a poor choice for a discovery
  study. It carries ample rare variation (median 51 rare coding variants per
  gene) because that follows from sample size, not from gene selection.
- **The mixed-direction arm is sign-balanced**, which is the idealized case and
  the best case for SKAT. Randomly drawn signs would land between the arms.
- **Genotyping error here is non-differential** by construction. Differential
  error, the failure mode that actually breaks real studies, is not covered.
- **SKATO and ACAT were not emitted for M2_lof_mis** in any run, only SKAT.
  Unexplained rather than understood.

---

## Corrections made during the work

Recorded because the reasoning matters more than the final numbers.

- Heterozygote excess at MAF >= 0.05 was the first QC instrument. It was
  abandoned: that set is disjoint from the variants under test, and the test has
  no power below roughly 1/sqrt(N) regardless of threshold. Replaced with
  mapping quality, Mendelian error and VQSLOD, which are frequency independent.
- Five 17q21.31 genes were initially labeled "controls". They are the author's
  own candidate genes, and under a simulated null every gene is null anyway, so
  no biological control was ever needed.
- Type I error was going to be pooled across genes with a Wilson interval. The
  panel has 14 genes in 16p11.2 and 13 in 17q21.31, about 31 independent units
  out of 56, so rates are computed per gene over independent replicates and
  reported as a distribution.
- Predicted that PC correction would fail worst in the rarest bins. The opposite
  is true; see finding 2.
- Justified indexing power by variant count on the grounds that gene length
  predicts it. Measured correlation is r = 0.118. The conclusion stands for a
  different reason.

## Layout

```
config/    panel, coordinates, mappability, REGENIE annotation/set-list/masks
scripts/   numbered, run in order
docs/      design notes and per-stage results
nextflow/  DSL2 pipeline, parallel over chromosomes and phenotype chunks
results/   calibration and power tables
```
