process CALIBRATION_SWEEP {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/calibration", mode: 'copy'

    input:
    tuple val(sample_id), path(bed), path(ground_truth), path(reference)

    output:
    tuple val(sample_id),
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
