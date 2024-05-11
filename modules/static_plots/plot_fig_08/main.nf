process PLOT_FIG_08 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_08", mode: 'copy'

    output:
    path "Fig_08.png", emit: "png"
    path "Fig_08.drawio", emit: "drawio"

    script:
    """
    wget -O Fig_08.png https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_09/Fig_9.png
    wget -O Fig_08.drawio https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_09/Fig_9.drawio
    """

    stub:
    """
    touch Fig_08.png
    touch Fig_08.drawio
    """
}