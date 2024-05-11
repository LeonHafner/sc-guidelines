process PLOT_FIG_01 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_01", mode: 'copy'

    output:
    path 'Fig_01.png', emit: 'png'
    path 'Fig_01.drawio', emit: 'drawio'
    
    script:
    """
    wget -O Fig_01.png https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_01/Fig_1.png
    wget -O Fig_01.drawio https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_01/Fig_1.drawio
    """

    stub:
    """
    touch Fig_01.png
    touch Fig_01.drawio
    """
}