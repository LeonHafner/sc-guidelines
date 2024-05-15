process DISTINCT {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/distinct'

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path('time.txt')

    script:
    """
    START=\$( date +%s )
    distinct.R \
        --input ${anndata} \
        --output results.tsv \
        --scenario dataset
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch results.tsv
    echo 46 > time.txt
    """
}