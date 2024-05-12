process SCVI_PROCESSING {
    tag "${meta.scenario}_${meta.run}"
    
    container 'leonhafner/python'
    
    input:
    tuple(val(meta), path(scvi_outputs, stageAs: "result_files/*"))

    output:
    tuple val(meta), path('scvi_values.tsv')

    script:
    template 'scvi_processing.py'

    stub:
    """
    touch scvi_values.tsv
    """
}