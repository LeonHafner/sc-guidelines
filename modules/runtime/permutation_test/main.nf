process PERMUTATION_TEST {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/python'

    input:
    tuple val(meta), path(sc_anndata), path(pb_anndata)

    output:
    tuple val(meta), path(sc_anndata), path(pb_anndata), path('time.txt')

    script:
    """
    START=\$( date +%s )
    permutation_test.py \
        --input ${pb_anndata} \
        --output results.tsv \
        --scenario dataset
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch results.tsv
    echo 48 > time.txt
    """
}