process MODKIT_VALIDATE {
    publishDir "${params.outdir}/validate", mode: 'copy'

    input:
    // Stage into separate subdirectories: both BAMs come from the same ALIGN module
    // and share the file name aligned.bam, which would otherwise collide
    path m6a_bam, stageAs: 'wt/*'
    path m6a_ground_truth
    path control_bam, stageAs: 'control/*'
    path control_ground_truth

    output:
    path "validate_output.txt"

    script:
    """
    modkit validate \\
        --bam-and-bed ${m6a_bam} ${m6a_ground_truth} \\
        --bam-and-bed ${control_bam} ${control_ground_truth} \\
        2>&1 | tee validate_output.txt
    """
}
