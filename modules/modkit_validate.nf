process MODKIT_VALIDATE {
    publishDir "${params.outdir}/validate", mode: 'copy'

    input:
    path m6a_bam
    path m6a_ground_truth
    path control_bam
    path control_ground_truth

    output:
    path "validate_output.txt"

    script:
    """
    modkit validate \
        --bam-and-bed ${m6a_bam} ${m6a_ground_truth} \
        --bam-and-bed ${control_bam} ${control_ground_truth} \
        2>&1 | tee validate_output.txt
    """
}
