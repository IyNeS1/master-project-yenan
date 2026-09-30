#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Parameters ──────────────────────────────────────────
// Modified sample (WT)
params.pod5                 = null
params.bam                  = null

// Control sample (IVT)
params.control_pod5         = null
params.control_bam          = null

// References and ground truth
params.reference            = null
params.ground_truth         = null
params.control_ground_truth = null

// Pipeline settings
params.outdir               = "results"
params.dorado_model         = "rna004_130bps_sup@v5.1.0"
params.kvalue               = 14
params.min_coverage         = 50
params.mod_threshold        = 15

// ── Include modules ─────────────────────────────────────
include { BASECALL as BASECALL_WT      } from './modules/basecall'
include { BASECALL as BASECALL_CONTROL } from './modules/basecall'
include { READ_QC                      } from './modules/read_qc'
include { ALIGN as ALIGN_WT            } from './modules/align'
include { ALIGN as ALIGN_CONTROL       } from './modules/align'
include { MODKIT_PILEUP as PILEUP_WT      } from './modules/modkit_pileup'
include { MODKIT_PILEUP as PILEUP_CONTROL } from './modules/modkit_pileup'
include { CALIBRATION_SWEEP            } from './modules/calibration_sweep'
include { CALIBRATION_SWEEP_IVT        } from './modules/calibration_sweep_ivt'
include { MODKIT_VALIDATE              } from './modules/modkit_validate'
include { ACCURACY_COMPARISON          } from './modules/accuracy_comparison'

// ── Workflow ────────────────────────────────────────────
workflow {

    ref_ch = Channel.fromPath(params.reference)

    // ── WT sample: basecall or use pre-made BAM ──
    // Dorado models are required only when basecalling from POD5
    if ((params.pod5 || params.control_pod5) && !params.dorado_models_dir) {
        error "Please provide --dorado_models_dir when starting from POD5 input"
    }
    dorado_models_ch = params.dorado_models_dir ? Channel.value(file(params.dorado_models_dir, checkIfExists: true)) : null

    if (params.pod5) {
        pod5_ch = Channel.fromPath(params.pod5, type: 'any')
        wt_bam  = BASECALL_WT(Channel.of('wt').combine(pod5_ch), dorado_models_ch)
    } else {
        wt_bam  = Channel.of('wt').combine(Channel.fromPath(params.bam))
    }

    // Read QC on WT
    READ_QC(wt_bam)

    // Align WT
    wt_align_in = wt_bam
        .map { sid, bam -> tuple(sid, params.kvalue, bam) }
        .combine(ref_ch)
    wt_aligned  = ALIGN_WT(wt_align_in)
    // wt_aligned: [sample_id, aligned.bam, aligned.bam.bai]

    // Modkit pileup on WT
    wt_pileup_in = wt_aligned.combine(ref_ch)
    wt_modkit    = PILEUP_WT(wt_pileup_in)
    // wt_modkit: [sample_id, modifications.bed, all_modifications.bed]

    // ── Evaluation: only if ground truth provided ──
    if (params.ground_truth) {
        gt_ch = Channel.fromPath(params.ground_truth)

        // Calibration sweep (uses all_modifications.bed)
        cal_in = wt_modkit
            .map { sid, mod_bed, all_bed -> tuple(sid, all_bed) }
            .combine(gt_ch)
            .combine(ref_ch)
        CALIBRATION_SWEEP(cal_in)

        // GLORI correlation
        // GLORI correlation only applies to real-world data with a GLORI CSV ground truth;
        // synthetic BED ground truths have no per-site modification levels to correlate
        if (params.ground_truth.toString().endsWith('.csv')) {
            acc_in = wt_modkit
                .map { sid, mod_bed, all_bed -> tuple(sid, all_bed) }
                .combine(gt_ch)
            ACCURACY_COMPARISON(acc_in)
        }
    }

    // ── Control sample: basecall + align + pileup ──
    if (params.control_pod5 || params.control_bam) {

        if (params.control_pod5) {
            ctrl_pod5_ch = Channel.fromPath(params.control_pod5, type: 'any')
            ctrl_bam     = BASECALL_CONTROL(Channel.of('control').combine(ctrl_pod5_ch), dorado_models_ch)
        } else {
            ctrl_bam     = Channel.of('control').combine(Channel.fromPath(params.control_bam))
        }

        // Align control
        ctrl_align_in = ctrl_bam
            .map { sid, bam -> tuple(sid, params.kvalue, bam) }
            .combine(ref_ch)
        ctrl_aligned  = ALIGN_CONTROL(ctrl_align_in)

        // Modkit pileup on control
        ctrl_pileup_in = ctrl_aligned.combine(ref_ch)
        ctrl_modkit    = PILEUP_CONTROL(ctrl_pileup_in)

        // IVT calibration sweep (only if control ground truth provided)
        if (params.control_ground_truth) {
            ctrl_gt_ch = Channel.fromPath(params.control_ground_truth)

            ivt_cal_in = ctrl_modkit
                .map { sid, mod_bed, all_bed -> tuple(sid, all_bed) }
                .combine(ctrl_gt_ch)
            CALIBRATION_SWEEP_IVT(ivt_cal_in)
        }

        // Modkit validate (needs both WT and control ground truth)
        if (params.ground_truth && params.control_ground_truth) {
            gt_ch2      = Channel.fromPath(params.ground_truth)
            ctrl_gt_ch2 = Channel.fromPath(params.control_ground_truth)

            MODKIT_VALIDATE(
                wt_aligned.map { sid, bam, bai -> bam },
                gt_ch2,
                ctrl_aligned.map { sid, bam, bai -> bam },
                ctrl_gt_ch2
            )
        }
    }
}
