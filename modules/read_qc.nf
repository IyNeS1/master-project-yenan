process READ_QC {
    tag "${sample_id}"
    stageInMode 'copy'
    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("qc_output/")

    script:
    """
    NanoPlot --ubam ${bam} --outdir qc_output/ --threads ${task.cpus}
    """
}
