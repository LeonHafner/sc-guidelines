process PLOT_FIG_S03 {
    tag "run_${meta.run}"

    container 'leonhafner/plotting'

    publishDir "${params.output}/Fig_S03", mode: 'copy'

    input:
    tuple val(meta), path(prc_files)

    output:
    path "Fig_S03_run${meta.run}.png"

    script:
    path_string = prc_files.join(";")
    template 'plot_fig_s03.R'

    stub:
    """
    touch Fig_S03_run${meta.run}.png
    """
}