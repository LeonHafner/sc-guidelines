process PLOT_FIG_07 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_07", mode: 'copy'

    input:
    path 'Fig_07.png'
    path 'Fig_07.drawio'

    output:
    path "Fig_07.png", emit: "png"
    path "Fig_07.drawio", emit: "drawio"

    script:
    // Empty process to transfer the static files to the output directory
    """
    """
}