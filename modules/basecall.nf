process BASECALL {
    tag "${sample_id}"
    label 'gpu'
    publishDir "${params.outdir}/${sample_id}/basecall", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)
    path models_dir

    output:
    tuple val(sample_id), path("calls.bam")

    script:
    """
    # Use pre-downloaded models so compute nodes do not need internet access
    dorado basecaller ${params.dorado_model} ${pod5} \\
        --models-directory ${models_dir} \\
        --modified-bases m6A_DRACH \\
        --batchsize 64 > calls.bam
    """
}
