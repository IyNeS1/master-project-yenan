process ACCURACY_COMPARISON {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/accuracy", mode: 'copy'

    input:
    tuple val(sample_id), path(bed), path(ground_truth)

    output:
    tuple val(sample_id),
          path("glori_correlation.png"),
          path("glori_accuracy_results.csv")

    script:
    """
    python ${projectDir}/bin/glori_correlation.py \
        --modkit ${bed} \
        --glori ${ground_truth} \
        --min_coverage ${params.min_coverage}
    """
}
