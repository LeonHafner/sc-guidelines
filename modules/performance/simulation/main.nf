process SIMULATION {
    tag "${meta.scenario}_${meta.run}"
    
    container 'leonhafner/splatter'

    input:
    val meta
    val n_genes

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    simulation.R \
        --scenario ${meta.scenario} \
        --n_genes ${n_genes} \
        --output ${meta.scenario}_${meta.run}.h5ad \
        --seed ${meta.run}
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}