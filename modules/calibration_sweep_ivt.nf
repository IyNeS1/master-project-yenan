process CALIBRATION_SWEEP_IVT {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/calibration", mode: 'copy'

    input:
    tuple val(sample_id), path(bed), path(ground_truth)

    output:
    tuple val(sample_id),
          path("calibration_sweep_ivt_combined.png"),
          path("calibration_sweep_ivt_fp_trend.png"),
          path("calibration_sweep_ivt_results.csv")

    script:
    """
    python ${projectDir}/bin/calibration_sweep_ivt.py \
        --m6a ${bed} \
        --ground_truth ${ground_truth} \
        --min_coverage ${params.min_coverage}
    """
}
