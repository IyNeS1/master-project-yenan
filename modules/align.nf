process ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), val(kvalue), path(bam), path(reference)

    output:
    tuple val(sample_id), path("aligned.bam"), path("aligned.bam.bai")

    script:
    """
    samtools fastq -T MM,ML ${bam} > reads.fastq
    minimap2 -ax splice -uf -k${kvalue} -y -t ${task.cpus} \
        ${reference} reads.fastq \
        | samtools sort -@ ${task.cpus} -o aligned.bam
    samtools index aligned.bam
    """
}
