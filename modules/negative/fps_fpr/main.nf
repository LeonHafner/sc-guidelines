process FPS_FPR {
    tag "${meta.scenario}_${meta.run}"
    
    container 'leonhafner/python'

    input:
    tuple val(meta), path(pvalues), path(scvi_values)

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}_fpr.tsv"), emit: 'fpr'
    tuple val(meta), path("${meta.scenario}_${meta.run}_fps.tsv"), emit: 'fps'

    script:
    template 'fps_fpr.py'

    stub:
    """
    touch ${meta.scenario}_${meta.run}_fpr.tsv
    touch ${meta.scenario}_${meta.run}_fps.tsv
    """
}