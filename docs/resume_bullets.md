# Resume bullets, scoped to what was actually done

## What is claimable

- Gene-based rare variant association implemented end to end on individual-level
  WGS: region slicing, consequence annotation, mask construction, REGENIE step 2
  with burden, SKAT, SKAT-O and ACAT.
- Type I error calibration measured across masks, MAF thresholds, case-control
  ratios and population structure, on 500 to 300 null replicates per arm.
- Power benchmarked across four genes, three effect sizes and two effect
  direction models.
- Genotype quality QC on rare variants using mapping quality, Mendelian error
  rate and VQSLOD, and its effect on mask composition.
- Nextflow DSL2 pipeline, executed and validated on the association path with
  confirmed parallel task execution.

## What is NOT claimable

- **Not deployed on AWS Batch.** The profile is written and untested.
- **Containers unverified.** The validated run used `-profile local` with conda
  binaries; no image was ever pulled.
- **No discovery.** Simulated phenotypes throughout. No gene is implicated.
- **Not biobank scale.** n = 2504.

---

## Bullets

### Statistical genetics emphasis (recommended)

> Built and benchmarked a gene-based rare variant association pipeline (REGENIE
> burden, SKAT, SKAT-O, ACAT) on 1000 Genomes 30x WGS, characterizing type I
> error across mask definitions, MAF thresholds and case-control ratios;
> demonstrated that genomic inflation factor is uninformative for sparse binary
> traits, where an uncorrected score test at 1:99 case:control showed lambda 0.06
> alongside a 26-fold excess of false positives at alpha 0.001, corrected by
> Firth.

> Showed that population structure confounding in rare variant burden tests
> scales with the AAF ceiling (lambda 33 uncapped versus 0.92 for singletons),
> and that PC correction was sufficient for all MAF < 0.01 masks but not for
> uncapped ones.

### Methods depth emphasis

> Demonstrated that REGENIE's default `--vc-MACthr 10` collapses ~86% of rare
> variants into a single burden term at n = 2504, converting SKAT into a burden
> test; disabling it recovered SKAT's advantage under mixed-direction effects for
> variant-rich genes (power 0.60 versus 0.36) but left SKAT underpowered for
> genes dominated by ultra-rare variants, characterizing the default as a
> bias-variance trade-off rather than a defect.

### Workflow engineering emphasis

> Implemented the pipeline as a Nextflow DSL2 workflow parallelized over
> chromosomes and phenotype chunks, validated end to end on the association and
> analysis stages; phenotype chunking added after the serial implementation hit
> a per-process open-file limit at 246 concurrent REGENIE outputs.

### One-line version

> Benchmarked type I error calibration and power of gene-based rare variant
> tests (REGENIE burden/SKAT-O/ACAT) on 1000 Genomes 30x WGS across mask
> definitions and MAF thresholds, implemented as a Nextflow pipeline.

---

## Interview material, in order of strength

1. **Lambda is not a sufficient diagnostic.** Measured lambda 0.056 with a
   26.5x false positive excess at alpha 0.001, simultaneously, under a null with
   randomly assigned labels. Explain the bimodality: degenerate statistic in the
   bulk, extreme tail from the few genes with a case carrier.

2. **The AAF ceiling is a bias control, not only a power knob.** Uncorrected
   inflation rose monotonically from lambda 0.92 at singletons to 33.1 uncapped.

3. **`--vc-MACthr` and sample size.** MAC 10 is MAF 0.002 at n = 2504 and
   MAF 0.00001 at n = 500,000. The same default is nearly inert in a biobank and
   dominant in a small cohort.

4. **Mask choice changes which genes have power, and why.** ZNF717 carries 32
   LoF variants, 16% of its coding variants, and has the panel's worst Mendelian
   error rate. Under a LoF mask it is the most powered gene in the panel, from
   mismapped reads.

5. **A prediction that failed.** Predicted PC correction would be weakest in the
   rarest bins because rare variants carry sharper structure. The opposite held,
   and the mechanism (confounding needs precisely measured, differentiated
   burden scores; singleton burdens are too noisy) is more interesting than the
   prediction would have been.

6. **A result that did not replicate.** The pilot found a clean burden/SKAT
   reversal on one gene. Three of four genes in the full grid disagreed, and one
   contradicted the implied recommendation. Worth telling: it is the difference
   between a pilot that locates a grid and a pilot mistaken for a result.
