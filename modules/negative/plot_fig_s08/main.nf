process PLOT_FIG_S08 {
    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S08", mode: 'copy'

    input:
    tuple val(meta), path(plotting_data_kang)
    tuple val(meta2), path(plotting_data_atlas_negative)

    output:
    path 'Fig_S08.png'

    script:
    template 'plot_fig_s08.R'

    stub:
    """
    touch Fig_S08.png
    """
}