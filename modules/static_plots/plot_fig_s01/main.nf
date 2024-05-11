process PLOT_FIG_S01 {
    container 'leonhafner/plotting'

    publishDir "${params.output}", mode: 'copy'

    output:
    path "Fig_S01.png"

    script:
    """
    wget -O Fig_S01.png https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_10/Fig_10.png
    """

    stub:
    """
    touch Fig_S01.png
    """
}