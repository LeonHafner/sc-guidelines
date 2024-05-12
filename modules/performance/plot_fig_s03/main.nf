process PLOT_FIG_S03 {
    tag "run_${meta.run}"

    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S03", mode: 'copy'

    input:
    tuple val(meta), val(path_list)

    output:
    path "Fig_S03_run${meta.run}.png"

    script:
    path_string = path_list.join(";")
    template 'plot_fig_s03.R'

    stub:
    """
    touch Fig_S03_run${meta.run}.png
    """
}