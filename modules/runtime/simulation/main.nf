process SIMULATION {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/splatter'

    input:
    val meta

    output:
    tuple val(meta), path("simulation.h5ad"), path("time.txt")

    script:
    def start = System.currentTimeMillis()
    """
    START=\$( date +%s )
    runtime-simulation.R \
        simulation.h5ad \
        ${meta.n_genes} \
        ${meta.n_cells} \
        ${meta.run}
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch simulation.h5ad
    echo 42 > time.txt
    """
}