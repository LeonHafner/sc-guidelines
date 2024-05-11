process DESEQ2 {
    tag "${meta.scenario}_${meta.n_cells}_${meta.n_genes}_${meta.run}"
    
    container 'leonhafner/deseq2'

    input:
    tuple val(meta), path(sc_anndata), path(pb_anndata)

    output:
    tuple val(meta), path(sc_anndata), path(pb_anndata), path('time.txt')

    script:
    """
    START=\$( date +%s )
    deseq2.R \
        --input ${pb_anndata} \
        --output results.tsv \
        --scenario dataset
    END=\$( date +%s )
    echo \$((\$END-\$START)) > time.txt
    """

    stub:
    """
    touch results.tsv
    echo 47 > time.txt
    """
}