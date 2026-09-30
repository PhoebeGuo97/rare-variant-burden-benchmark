#!/usr/bin/env nextflow
/*
 * Rare variant burden pipeline: calibration and power benchmarking.
 *
 * Two entry points.
 *
 *   default      full pipeline: slice -> annotate -> pgen -> assoc -> analyze
 *   assoc_only   start from the existing pgen and run only the association and
 *                analysis stages
 *
 * assoc_only exists so the pipeline can be tested for real in minutes instead
 * of re-downloading 2.3 GB. It covers the stages that carry the parallelism and
 * the stages most likely to break (phenotype chunking, REGENIE invocation,
 * aggregation), so a passing assoc_only run is meaningful evidence rather than
 * a syntax check.
 *
 *   nextflow run main_v2.nf -entry assoc_only -profile local \
 *       --pheno config/pheno_smoke.tsv --outdir results_nf_smoke
 */

nextflow.enable.dsl = 2

params.root        = "$projectDir/.."
params.panel_bed   = "${params.root}/config/panel_regions.bed"
params.keep        = "${params.root}/config/keep_2504.txt"
params.gff         = "${params.root}/ref/panel.chr.gff3.gz"
params.ref         = "${params.root}/ref/GRCh38.panel_chroms.fa"
params.annot       = "${params.root}/config/regenie_annot.tsv"
params.setlist     = "${params.root}/config/regenie_setlist.tsv"
params.maskdef     = "${params.root}/config/regenie_masks.def"
params.covar       = "${params.root}/config/covariates.tsv"
params.pheno       = "${params.root}/config/pheno_qt_null_qt.tsv"
params.pgen_prefix = "${params.root}/data/plink/panel"
params.analyze_py  = "${params.root}/scripts/11_analyze.py"
params.outdir      = "results_nf"

params.vcf_base    = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"
params.vcf_prefix  = "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05"

params.aaf_bins    = "0.001,0.005,0.01"
params.vc_tests    = "skat,skato,skato-acat,acatv,acato-full"
params.vc_maxaaf   = 0.01
params.vc_macthr   = 1        // 10 is the REGENIE default; see docs/power_result.md
params.chunk_size  = 100
params.pcs         = "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"


process SLICE {
    tag "$chrom"
    container params.container_bcftools

    input:
    tuple val(chrom), path(bed), path(keep), path(ref), path(fai)

    output:
    tuple val(chrom), path("${chrom}.vcf.gz"), path("${chrom}.vcf.gz.tbi")

    script:
    """
    awk -v c=$chrom '\$1==c' $bed > ${chrom}.bed
    bcftools view -R ${chrom}.bed -S $keep --force-samples -f PASS -c 1 -Ou \\
        ${params.vcf_base}/${params.vcf_prefix}_${chrom}.recalibrated_variants.vcf.gz \\
      | bcftools norm -m -any -f $ref -Ou \\
      | bcftools norm -d exact -Oz -o ${chrom}.vcf.gz
    tabix -p vcf ${chrom}.vcf.gz
    """
}

process ANNOTATE {
    tag "$chrom"
    container params.container_bcftools

    input:
    tuple val(chrom), path(vcf), path(tbi), path(gff), path(ref), path(fai)

    output:
    path "${chrom}.csq.vcf.gz"
    path "${chrom}.csq.vcf.gz.tbi"

    script:
    """
    bcftools csq -f $ref -g $gff --local-csq --ncsq 30 -Ou $vcf \\
      | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${chrom}.csq.vcf.gz
    tabix -p vcf ${chrom}.csq.vcf.gz
    """
}

process CONCAT {
    container params.container_bcftools

    input:
    path vcfs
    path tbis

    output:
    path "merged.vcf.gz"

    script:
    """
    ls *.csq.vcf.gz | sort -V > list.txt
    bcftools concat -f list.txt -Oz -o merged.vcf.gz
    """
}

process PGEN {
    container params.container_plink2
    publishDir "${params.outdir}/plink", mode: 'copy'

    input:
    path vcf

    output:
    tuple path("panel.pgen"), path("panel.pvar"), path("panel.psam")

    script:
    """
    plink2 --vcf $vcf --double-id --vcf-half-call m \\
           --make-pgen --sort-vars --out panel
    """
}

/*
 * Chunking is load bearing, not tidiness. REGENIE opens one output file per
 * phenotype and holds them all open; the serial run died at trait 246 against
 * a 256 open-file limit. Container defaults vary, so bound it explicitly.
 */
process SPLIT_PHENO {
    container params.container_python

    input:
    path pheno

    output:
    path "chunk_*.txt"

    script:
    """
#!/usr/bin/env python3
with open("$pheno") as fh:
    cols = fh.readline().rstrip("\\n").split("\\t")[2:]
n = ${params.chunk_size}
for i in range(0, len(cols), n):
    with open("chunk_%04d.txt" % (i // n + 1), "w") as out:
        out.write(",".join(cols[i:i + n]))
    """
}

process ASSOC {
    tag "${chunk.baseName}"
    container params.container_regenie
    cpus 4
    publishDir "${params.outdir}/regenie", mode: 'copy'

    input:
    tuple path(chunk), path(pgen), path(pvar), path(psam),
          path(pheno), path(covar), path(annot), path(setlist), path(maskdef)

    output:
    path "*.regenie.gz", optional: true, emit: results
    path "*.log",                       emit: logs

    script:
    """
    ulimit -n 10240 || true
    regenie --step 2 \\
      --pgen ${pgen.baseName} \\
      --phenoFile $pheno --phenoColList \$(cat $chunk) \\
      --covarFile $covar --covarColList ${params.pcs} \\
      --anno-file $annot --set-list $setlist --mask-def $maskdef \\
      --aaf-bins ${params.aaf_bins} \\
      --vc-tests ${params.vc_tests} \\
      --vc-maxAAF ${params.vc_maxaaf} \\
      --vc-MACthr ${params.vc_macthr} \\
      --build-mask max --minMAC 1 --ignore-pred --bsize 200 --gz \\
      --threads ${task.cpus} \\
      --out ${chunk.baseName}
    """
}

process CALIBRATION {
    container params.container_python
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path results
    path script

    output:
    path "calibration.tsv"

    script:
    """
    mkdir -p regenie
    # 11_analyze.py discovers arms by filename prefix, so the chunk outputs
    # are renamed into the a1_qt arm. The Arm B grouping is skipped because
    # its phenotype file is not staged here, which is correct: this entry
    # point only runs the unstructured null.
    for f in *.regenie.gz; do cp "\$f" regenie/a1_qt_\$f; done
    python3 $script --resdir regenie --out calibration.tsv
    """
}


workflow prep {
    main:
    chroms = Channel.fromPath(params.panel_bed)
                    .splitCsv(sep: '\t').map { it[0] }.unique()
    ref = file(params.ref); fai = file("${params.ref}.fai")

    sliced = SLICE(chroms.map { c ->
        tuple(c, file(params.panel_bed), file(params.keep), ref, fai) })

    ann = ANNOTATE(sliced.map { c, v, t -> tuple(c, v, t, file(params.gff), ref, fai) })
    merged = CONCAT(ann[0].collect(), ann[1].collect())
    pg = PGEN(merged)

    emit:
    pgen = pg
}

workflow assoc {
    take:
    pgen

    main:
    chunks = SPLIT_PHENO(file(params.pheno)).flatten()
    inp = chunks.combine(pgen).map { ch, pg, pv, ps ->
        tuple(ch, pg, pv, ps, file(params.pheno), file(params.covar),
              file(params.annot), file(params.setlist), file(params.maskdef)) }
    out = ASSOC(inp)
    CALIBRATION(out.results.collect(), file(params.analyze_py))
}

/*
 * Single entry workflow, branching on params.stage.
 *
 * Nextflow 26's strict parser removed `-entry`, so a named workflow can no
 * longer be selected from the command line. The replacement it points to is a
 * parameter read inside the entry workflow:
 *
 *   --stage all          full pipeline (default)
 *   --stage assoc_only   start from an existing pgen, skipping the 2.3 GB
 *                        download; used to test the pipeline in minutes
 */
params.stage = 'all'

workflow {
    if( params.stage == 'assoc_only' ) {
        pgen = Channel.of(tuple(file("${params.pgen_prefix}.pgen"),
                                file("${params.pgen_prefix}.pvar"),
                                file("${params.pgen_prefix}.psam")))
        assoc(pgen)
    }
    else {
        prep()
        assoc(prep.out.pgen)
    }
}
