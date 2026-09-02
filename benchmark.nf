#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Parameters ──────────────────────────────────────────
params.pod5             = null
params.reference        = null
params.ground_truth     = null
params.outdir           = "results_benchmark"
params.models_dir       = "/data/projects/punim2974/software/dorado_models"

// 4 models to benchmark (comma-separated)
params.models = "rna004_130bps_sup@v5.1.0,rna004_130bps_sup@v5.2.0,rna004_130bps_hac@v5.1.0,rna004_130bps_hac@v5.2.0"

// 3 minimap2 k-mer sizes to benchmark (comma-separated)
params.kvalues = "12,14,16"

// ── Helper: clean model name for directory ──────────────
def modelLabel(model) {
    // rna004_130bps_sup@v5.1.0 -> sup_v5.1.0
    def m = model.replaceAll(/rna004_130bps_/, '').replaceAll(/@/, '_')
    return m
}

// ── Step 1: Basecalling (GPU, one per model) ────────────
process BASECALL {
    tag "${model_label}"
    publishDir "${params.outdir}/${model_label}/basecall", mode: 'copy'

    input:
    tuple val(model_name), val(model_label), path(pod5)

    output:
    tuple val(model_name), val(model_label), path("calls.bam")

    script:
    """
    dorado basecaller ${model_name} ${pod5} \
        --modified-bases m6A_DRACH \
        --batchsize 64 \
        --models-directory ${params.models_dir} \
        > calls.bam
    """
}

// ── Step 2: Alignment (one per model × k-value) ────────
process ALIGN {
    tag "${model_label}_k${kvalue}"
    publishDir "${params.outdir}/${model_label}_k${kvalue}/aligned", mode: 'copy'

    input:
    tuple val(model_name), val(model_label), val(kvalue), path(bam), path(reference)

    output:
    tuple val(model_name), val(model_label), val(kvalue), path("aligned.bam"), path("aligned.bam.bai")

    script:
    """
    samtools fastq -T MM,ML ${bam} > reads.fastq
    minimap2 -ax splice -uf -k${kvalue} -y -t ${task.cpus} ${reference} reads.fastq \
        | samtools sort -@ ${task.cpus} -o aligned.bam
    samtools index aligned.bam
    """
}

// ── Step 3: Modkit Pileup ───────────────────────────────
process MODKIT_PILEUP {
    tag "${model_label}_k${kvalue}"
    publishDir "${params.outdir}/${model_label}_k${kvalue}/modkit", mode: 'copy'

    input:
    tuple val(model_name), val(model_label), val(kvalue), path(bam), path(bai), path(reference)

    output:
    tuple val(model_name), val(model_label), val(kvalue), path("all_modifications.bed")

    script:
    """
    modkit pileup --reference ${reference} ${bam} all_modifications.bed
    """
}

// ── Step 4: Calibration Sweep ───────────────────────────
process CALIBRATION_SWEEP {
    tag "${model_label}_k${kvalue}"
    publishDir "${params.outdir}/${model_label}_k${kvalue}/calibration", mode: 'copy'

    input:
    tuple val(model_name), val(model_label), val(kvalue), path(bed), path(ground_truth), path(reference)

    output:
    tuple val(model_name), val(model_label), val(kvalue),
          path("calibration_sweep_results.csv"),
          path("calibration_sweep_combined.png"),
          path("calibration_sweep_separate.png")

    script:
    """
    python ${projectDir}/bin/calibration_sweep.py \
        --m6a ${bed} \
        --ground_truth ${ground_truth} \
        --reference ${reference}
    """
}

// ── Step 5: Summary table across all combinations ───────
process SUMMARISE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path "benchmark_summary.csv"

    script:
    script:
    """
    first=\$(ls *.csv | head -1)
    head -1 \$first | sed 's/^/model,kvalue,/' > benchmark_summary.csv
    for f in *.csv; do
        base=\$(basename \$f .csv)
        model=\$(echo \$base | sed 's/_k[0-9]*//')
        kval=\$(echo \$base | grep -o 'k[0-9]*' | sed 's/k//')
        tail -n +2 \$f | while IFS= read -r line; do
            echo "\${model},\${kval},\${line}"
        done
    done >> benchmark_summary.csv
    """
}
workflow {

    // Value channels (reusable)
    ref_ch = Channel.fromPath(params.reference).first()
    gt_ch  = Channel.fromPath(params.ground_truth).first()
    pod5_ch = Channel.fromPath(params.pod5, type: 'any').first()

    // Create model channel: [model_name, model_label]
    model_ch = Channel
        .of(params.models.tokenize(','))
        .flatten()
        .map { m -> tuple(m, modelLabel(m)) }

    // Create k-value channel
    kvalue_ch = Channel
        .of(params.kvalues.tokenize(','))
        .flatten()

    // ── Basecall: 4 runs (one per model) ──
    basecall_input = model_ch.combine(pod5_ch)
    basecalled = BASECALL(basecall_input)

    // ── Align: 12 runs (4 models × 3 k-values) ──
    align_input = basecalled
        .combine(kvalue_ch)
        .map { model_name, model_label, bam, kval ->
            tuple(model_name, model_label, kval, bam)
        }
        .combine(ref_ch)

    aligned = ALIGN(align_input)

    // ── Modkit: 12 runs ──
    modkit_input = aligned.combine(ref_ch)

    modkit_out = MODKIT_PILEUP(modkit_input)

    // ── Calibration: 12 runs ──
    cal_input = modkit_out
        .combine(gt_ch)
        .combine(ref_ch)

    cal_out = CALIBRATION_SWEEP(cal_input)

    // ── Summary: collect all calibration CSVs ──
    all_csvs = cal_out
        .map { model_name, model_label, kval, csv, png1, png2 ->
            file(csv).copyTo("${csv.parent}/${model_label}_k${kval}.csv")
        }
        .collect()

    SUMMARISE(all_csvs)
}
