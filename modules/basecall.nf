process BASECALL {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/basecall", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5)

    output:
    tuple val(sample_id), path("calls.bam")

    script:
    """
    dorado basecaller ${params.dorado_model} ${pod5} \
        --modified-bases m6A_DRACH --batchsize 64 > calls.bam
    """
}
