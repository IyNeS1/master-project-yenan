#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam_dir      = null
params.reference    = null
params.ground_truth = null
params.outdir       = "results_benchmark"
params.kvalues      = "14,16"
params.min_coverage = 50
params.mod_threshold = 80

def modelLabel(bam_name) {
    return bam_name.replaceAll(/_calls\.bam$/, '')
}

// ── Reuse shared modules ────────────────────────────────
include { ALIGN             } from './modules/align'
include { MODKIT_PILEUP     } from './modules/modkit_pileup'
include { CALIBRATION_SWEEP } from './modules/calibration_sweep'
include { SUMMARISE         } from './modules/summarise'

// ── Workflow ────────────────────────────────────────────
workflow {
    ref_ch = Channel.fromPath(params.reference).first()
    gt_ch  = Channel.fromPath(params.ground_truth).first()

    bam_ch = Channel
        .fromPath("${params.bam_dir}/*_calls.bam")
        .map { bam -> tuple(modelLabel(bam.name), bam) }

    kvalue_ch = Channel
        .of(params.kvalues.toString().split(','))
        .flatten()

    // Cartesian: models × k-values
    align_input = bam_ch
        .combine(kvalue_ch)
        .map { model, bam, kval ->
            tuple("${model}_k${kval}", kval, bam)
        }
        .combine(ref_ch)

    aligned = ALIGN(align_input)

    modkit_input = aligned.combine(ref_ch)
    modkit_out = MODKIT_PILEUP(modkit_input)

    cal_input = modkit_out
        .map { sid, mod_bed, all_bed -> tuple(sid, all_bed) }
        .combine(gt_ch)
        .combine(ref_ch)
    cal_out = CALIBRATION_SWEEP(cal_input)

    // Collect CSVs for summary
    all_csvs = cal_out
        .map { sid, csv, png1, png2 ->
            file(csv).copyTo("${csv.parent}/${sid}_calibration_sweep_results.csv")
        }
        .collect()

    SUMMARISE(all_csvs)
}
