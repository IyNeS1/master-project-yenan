#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Parameters ──────────────────────────────────────────
// Modified sample (修饰样本)
params.pod5             = null
params.bam              = null

// Control sample (对照样本)
params.control_pod5     = null
params.control_bam      = null

// References
params.reference            = null
params.ground_truth         = null      // m6A ground truth BED
params.control_ground_truth = null      // canonical A ground truth BED
params.outdir               = "results"
params.dorado_model         = "rna004_130bps_sup@v5.1.0"

// ── Step 1: Basecalling – modified sample (requires GPU) ─
process BASECALL {
    publishDir "${params.outdir}/basecall", mode: 'copy'
    conda null

    input:
    path pod5

    output:
    path "calls.bam"

    script:
    """
    dorado basecaller ${params.dorado_model} ${pod5} --modified-bases m6A_DRACH --batchsize 64 > calls.bam
    """
}

// ── Step 1b: Basecalling – control sample (requires GPU) ─
process BASECALL_CONTROL {
    publishDir "${params.outdir}/control/basecall", mode: 'copy'
    conda null

    input:
    path pod5

    output:
    path "calls.bam"

    script:
    """
    dorado basecaller ${params.dorado_model} ${pod5} --modified-bases m6A_DRACH --batchsize 64 > calls.bam
    """
}

// ── Step 2: Read QC ─────────────────────────────────────
process READ_QC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path bam

    output:
    path "qc_output/"

    script:
    """
    export PYTHONPATH=""
    /home/yenanwang/.conda/envs/nanopore_env/bin/NanoPlot --ubam ${bam} --outdir qc_output/ --threads 4
    """
}

// ── Step 3: Sequence Alignment – modified sample ────────
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
    /home/yenanwang/.conda/envs/nanopore_env/bin/samtools fastq -T MM,ML ${bam} > reads.fastq
    /home/yenanwang/.conda/envs/nanopore_env/bin/minimap2 -ax map-ont -y ${reference} reads.fastq \
        | /home/yenanwang/.conda/envs/nanopore_env/bin/samtools sort -o aligned.bam
    /home/yenanwang/.conda/envs/nanopore_env/bin/samtools index aligned.bam
    """
}

// ── Step 3b: Sequence Alignment – control sample ────────
process ALIGN_CONTROL {
    publishDir "${params.outdir}/control/aligned", mode: 'copy'

    input:
    path bam
    path reference

    output:
    path "control_aligned.bam"
    path "control_aligned.bam.bai"

    script:
    """
    /home/yenanwang/.conda/envs/nanopore_env/bin/samtools fastq -T MM,ML ${bam} > reads.fastq
    /home/yenanwang/.conda/envs/nanopore_env/bin/minimap2 -ax map-ont -y ${reference} reads.fastq \
        | /home/yenanwang/.conda/envs/nanopore_env/bin/samtools sort -o control_aligned.bam
    /home/yenanwang/.conda/envs/nanopore_env/bin/samtools index control_aligned.bam
    """
}

// ── Step 4: Modification Detection (bedMethyl) ──────────
process MODKIT_PILEUP {
    publishDir "${params.outdir}/modkit", mode: 'copy'

    input:
    path bam
    path reference

    output:
    path "modifications.bed"

    script:
    """
    /home/yenanwang/.conda/envs/nanopore_env/bin/samtools index ${bam}
    modkit pileup --reference ${reference} ${bam} all_modifications.bed
    awk -F'\\t' '\$4 == "a" && \$11 > 80' all_modifications.bed > modifications.bed
    """
}

// ── Step 5a: Calibration Sweep ──────────────────────────
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
    export PYTHONPATH="" && /home/yenanwang/.conda/envs/nanopore_env/bin/python ${projectDir}/bin/calibration_sweep.py \
        --m6a ${bed} \
        --ground_truth ${ground_truth}
    """
}

// ── Step 5b: Modkit Validate (ONT benchmark method) ─────
process MODKIT_VALIDATE {
    publishDir "${params.outdir}/validate", mode: 'copy'

    input:
    path m6a_bam
    path m6a_ground_truth
    path control_bam
    path control_ground_truth

    output:
    path "validate_output.txt"

    script:
    """
    /home/yenanwang/.conda/envs/nanopore_env/bin/modkit validate \
        --bam-and-bed ${m6a_bam} ${m6a_ground_truth} \
        --bam-and-bed ${control_bam} ${control_ground_truth} \
        2>&1 | tee validate_output.txt
    """
}

// ── Workflow ────────────────────────────────────────────
workflow {

    ref_ch          = Channel.fromPath(params.reference)
    ground_truth_ch = Channel.fromPath(params.ground_truth)

    // ── Modified sample (m6A): basecall or use pre-made BAM ──
    if (params.pod5) {
        pod5_ch = Channel.fromPath(params.pod5)
        m6a_bam = BASECALL(pod5_ch)
    } else {
        m6a_bam = Channel.fromPath(params.bam)
    }

    READ_QC(m6a_bam)
    ALIGN(m6a_bam, ref_ch)

    // ── Step 4 + 5a: pileup-based evaluation (m6A only) ──
    modkit_out = MODKIT_PILEUP(ALIGN.out[0], ref_ch)
    CALIBRATION_SWEEP(modkit_out, ground_truth_ch)

    // ── Control sample: basecall + align ──
    if (params.control_pod5 || params.control_bam) {

        if (params.control_pod5) {
            ctrl_pod5_ch = Channel.fromPath(params.control_pod5)
            ctrl_bam     = BASECALL_CONTROL(ctrl_pod5_ch)
        } else {
            ctrl_bam     = Channel.fromPath(params.control_bam)
        }

        ALIGN_CONTROL(ctrl_bam, ref_ch)

        // ── Step 5b: modkit validate (m6A + control) ──
        ctrl_gt_ch = Channel.fromPath(params.control_ground_truth)

        MODKIT_VALIDATE(
            ALIGN.out[0],
            ground_truth_ch,
            ALIGN_CONTROL.out[0],
            ctrl_gt_ch
        )
    }
}
