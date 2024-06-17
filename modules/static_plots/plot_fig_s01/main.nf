process PLOT_FIG_S01 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    input:
    path 'Fig_S01.png'

    output:
    path "Fig_S01.png"

    script:
    // Empty process to transfer the static files to the output directory
    """
    """
}