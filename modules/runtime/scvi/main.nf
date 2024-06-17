process SCVI {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/scvi'

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path('time.txt')

    script:
    """
    START=\$( date +%s )
    scvi-de.py \
        --input ${anndata} \
        --output results.tsv \
        --scenario dataset
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch results.tsv
    echo 50 > time.txt
    """
}