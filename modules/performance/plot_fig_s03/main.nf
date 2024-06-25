process PLOT_FIG_S03 {
    tag "${meta.scenario}_${meta.run}"
    
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S03", mode: 'copy'

    input:
    tuple val(meta), path(input_anndata)

    output:
    path "Fig_S03_run${meta.run}.png"

    script:
    template 'plot_fig_s03.R'

    stub:
    """
    touch Fig_S03_run${meta.run}.png
    """
}