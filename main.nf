#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam       = null
params.reference = null
params.outdir    = "results"
params.ground_truth = null

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

// Calibration sweep
process CALIBRATION_SWEEP {
    publishDir "${params.outdir}/calibration", mode: 'copy'

    input:
    path bed
    path ground_truth

    output:
    path "calibration_sweep_combined.png"
    path "calibration_sweep_separate.png"

    script:
    """
    python ${projectDir}/calibration_sweep.py \
      --m6a ${bed} \
      --ground_truth ${ground_truth}
    """
}

workflow {
    bam_ch          = Channel.fromPath(params.bam)
    ref_ch          = Channel.fromPath(params.reference)
    ground_truth_ch = Channel.fromPath(params.ground_truth)

    READ_QC(bam_ch)
    ALIGN(bam_ch, ref_ch)
    modkit_out = MODKIT_PILEUP(bam_ch, ref_ch)
    CALIBRATION_SWEEP(modkit_out, ground_truth_ch)
}