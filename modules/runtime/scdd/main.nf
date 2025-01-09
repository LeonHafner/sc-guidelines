process SCDD {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/scdd'

    input:
    tuple val(meta), path(anndata)

    output:
    tuple val(meta), path('time.txt')

    script:
    """
    START=\$( date +%s )
    scdd.R \
        --input ${anndata} \
        --output results.tsv
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch results.tsv
    echo 53 > time.txt
    """
}