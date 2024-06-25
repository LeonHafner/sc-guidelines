process PLOT_FIG_S07 {
    tag "run_${meta.run}"

    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S07", mode: 'copy'

    input:
    tuple val(meta), path(prc_files)

    output:
    path "Fig_S07_run${meta.run}.png"

    script:
    path_string = prc_files.join(";")
    template 'plot_fig_s07.R'

    stub:
    """
    touch Fig_S07_run${meta.run}.png
    """
}