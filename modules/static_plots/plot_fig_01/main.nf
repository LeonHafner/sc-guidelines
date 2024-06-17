process PLOT_FIG_01 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_01", mode: 'copy'

    input:
    path 'Fig_01.png'
    path 'Fig_01.drawio'

    output:
    path 'Fig_01.png', emit: 'png'
    path 'Fig_01.drawio', emit: 'drawio'

    script:
    // Empty process to transfer the static files to the output directory
    """
    """
}