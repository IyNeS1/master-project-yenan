#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam       = null
params.reference = null
params.outdir    = "results"

// Step3：Read QC
process READ_QC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path bam

    output:
    path "qc_output/"

    script:
    """
    NanoPlot --bam ${bam} --outdir qc_output/ --threads 4
    """
}

// Step4：Sequence Alignment
process ALIGN {
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    path bam
    path reference

    output:
    path "aligned.bam"
    path "aligned.bam.bai"

    script:
    """
    samtools fastq ${bam} > reads.fastq
    minimap2 -ax map-ont ${reference} reads.fastq | samtools sort -o aligned.bam
    samtools index aligned.bam
    """
}

// Step5+6：Modification Detection + bedMethyl Output
process MODKIT_PILEUP {
    publishDir "${params.outdir}/modkit", mode: 'copy'

    input:
    path bam
    path reference

    output:
    path "modifications.bed"

    script:
    """
    samtools index ${bam}
    modkit pileup --reference ${reference} ${bam} all_modifications.bed
    grep -E \$'\\ta\\t' all_modifications.bed | awk -F'\\t' '\$11 > 80' > modifications.bed
    """
}

workflow {
    bam_ch = Channel.fromPath(params.bam)
    ref_ch = Channel.fromPath(params.reference)

    READ_QC(bam_ch)
    ALIGN(bam_ch, ref_ch)
    MODKIT_PILEUP(bam_ch, ref_ch)
}