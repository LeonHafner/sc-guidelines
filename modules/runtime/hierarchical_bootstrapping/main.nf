process HIERARCHICAL_BOOTSTRAPPING {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/python'

    input:
    tuple val(meta), path(sc_anndata), path(pb_anndata)

    output:
    tuple val(meta), path(sc_anndata), path(pb_anndata), path('time.txt')

    script:
    """
    START=\$( date +%s )
    hierarchical_bootstrapping.py \
        --input ${sc_anndata} \
        --output1 results_pre.tsv \
        --output2 results.tsv \
        --scenario dataset
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch results_pre.tsv
    touch results.tsv
    echo 49 > time.txt
    """
}