# Mask construction: what the panel can and cannot support

Masks built from `bcftools csq`, protein_coding transcripts only, most severe
consequence per gene. 56 genes in the set list.

## Testable genes per mask

Number of genes with at least N qualifying variants:

| AAF bin | N | M1_lof | M2_lof_mis | M4_all |
| --- | --- | --- | --- | --- |
| 0.001 | 1 | 43 | 55 | 56 |
| 0.001 | 3 | 18 | 55 | 55 |
| 0.001 | 10 | 1 | 50 | 54 |
| 0.01 | 1 | 46 | 56 | 56 |
| 0.01 | 3 | 21 | 55 | 56 |
| 0.01 | 10 | 2 | 51 | 55 |

**M1_lof is not a testable mask on this panel.** At MAF < 0.01 only two genes
carry ten or more LoF variants and only 21 carry three. A gene-based test on one
or two variants is a single-variant test with extra steps, and its null
distribution is not the one SKAT-O assumes. M1 stays in the grid because the
comparison is informative, but every M1 result is reported with its variant
count attached and nothing is concluded from M1 alone.

M2 and M4 are well populated: 50 to 55 genes at a floor of ten.

This is the honest ceiling of a 56 gene panel in 2504 samples, and it is worth
stating in the write-up rather than presenting three masks as equals.

## LoF fraction flags artifact and real biology, and the QC columns tell them apart

Raw LoF counts per kb are confounded by intron fraction: HBB and STH top that
ranking purely because they are small and almost entirely coding, while DLGAP1
at 959 kb looks LoF-poor because it is mostly intron. The internally normalized
version is LoF as a percentage of all rare coding variants in the gene.

| Gene | coding vars | LoF | LoF % | VQSLOD | ME>0% | median MQ |
| --- | --- | --- | --- | --- | --- | --- |
| HBB | 28 | 10 | 35.7 | +16.21 | 2.0 | 60.0 |
| ARL17A | 17 | 3 | 17.6 | **-6.18** | 12.4 | **36.2** |
| ZNF717 | 199 | 32 | 16.1 | **-2.13** | **17.7** | 57.4 |
| GOSR2 | 45 | 6 | 13.3 | +15.74 | 0.5 | 60.0 |

Median LoF fraction is **6.27% in the seven VQSLOD < 0 genes against 3.25% in
the other 46**, so poor genotype quality does inflate the LoF mask specifically,
roughly two fold, even though it deflates total variant counts.

Both effects are real and they are not in tension: duplicated regions get fewer
calls overall, but a larger share of the calls that survive are called as
loss-of-function, which is what mismapped reads produce.

**HBB is the control that makes the statistic trustworthy.** It has the highest
LoF fraction in the panel at 35.7%, and perfect quality: VQSLOD +16.21, MQ 60,
2.0% Mendelian error. That is real biology. HBB carries beta-thalassemia and
sickle alleles segregating at appreciable frequency across several of these
superpopulations under balancing selection for malaria resistance. A pipeline
that flagged genes on LoF enrichment alone would throw HBB out. The quality
columns are what separate it from ARL17A and ZNF717.

## Among large genes the outlier is unambiguous

Restricting to genes over 100 kb, where intron fraction is roughly comparable:

| Gene | kb | coding vars | LoF % | VQSLOD |
| --- | --- | --- | --- | --- |
| ZNF717 | 107 | 199 | **16.1** | **-2.13** |
| KCNIP1 | 383 | 32 | 9.4 | +16.14 |
| DYRK1A | 161 | 48 | 8.3 | +16.06 |
| MAPT | 134 | 96 | 3.1 | +15.52 |
| FBN1 | 237 | 204 | 1.5 | +16.18 |
| DLGAP1 | 959 | 73 | 1.4 | +15.95 |
| INPP5D | 148 | 107 | 0.0 | +15.56 |

ZNF717 carries 32 of the panel's LoF variants, more than ten times MAPT or FBN1
despite being smaller than both, and it is the gene with the worst Mendelian
error rate in the panel. This is the concrete answer to "when does mask choice
change the result": under M1_lof, ZNF717 is the most powered gene in the panel
and its power comes from mismapped reads. Under M4_all it is unremarkable.

## Under-calling confirmed in the mask counts

ARL17B has 3 qualifying variants in M4_all across 87 kb. ARL17A has 17 across
80 kb. Ordinary genes of that size carry 50 to 200. This is the same
under-calling seen in the QC step, now visible in the analysis set itself.

## Decisions recorded

- Only `protein_coding` transcript consequences count. NMD and retained_intron
  biotypes are dropped; a stop_gained on an NMD target is not gene LoF.
- Most severe consequence across the gene's protein_coding transcripts, not
  canonical or MANE only. More inclusive, and it does raise LoF counts.
- M3_lof_damaging (CADD >= 20) is not built. Deferred, and dropped entirely if
  time runs short; three masks already demonstrate mask sensitivity.
