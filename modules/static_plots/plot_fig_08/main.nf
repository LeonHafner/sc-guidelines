process PLOT_FIG_08 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_08", mode: 'copy'

    input:
    path 'Fig_08.png'
    path 'Fig_08.drawio'

    output:
    path "Fig_08.png", emit: "png"
    path "Fig_08.drawio", emit: "drawio"

    script:
    // Empty process to transfer the static files to the output directory
    """
    """
}