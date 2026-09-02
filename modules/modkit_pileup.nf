process MODKIT_PILEUP {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/modkit", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference)

    output:
    tuple val(sample_id), path("modifications.bed"), path("all_modifications.bed")

    script:
    """
    modkit pileup --reference ${reference} ${bam} all_modifications.bed
    awk -F'\\t' '\$4 == "a" && \$5 >= ${params.min_coverage} && \$11 > ${params.mod_threshold}' \
        all_modifications.bed > modifications.bed
    """
}
