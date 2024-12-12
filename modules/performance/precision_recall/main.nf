process PRECISION_RECALL {
    tag "${meta.scenario}_${meta.run}"
    
    container 'leonhafner/python'

    input:
    tuple val(meta), path(pvalues, stageAs: 'pvalues.tsv'), path(ground_truth, stageAs: 'ground_truth.tsv')

    output:
    tuple val(meta), path("prc_${meta.scenario}_${meta.run}_*.tsv"), emit: 'prc'
    tuple val(meta), path("auc_${meta.scenario}_${meta.run}.tsv"), emit: 'auc'

    script:
    template 'precision_recall.py'

    stub:
    """
    touch prc_${meta.scenario}_${meta.run}_deseq2.tsv
    touch prc_${meta.scenario}_${meta.run}_mast.tsv
    touch prc_${meta.scenario}_${meta.run}_scvi.tsv
    touch prc_${meta.scenario}_${meta.run}_dream.tsv
    touch prc_${meta.scenario}_${meta.run}_distinct.tsv
    touch prc_${meta.scenario}_${meta.run}_hierarchical-bootstrapping.tsv
    touch prc_${meta.scenario}_${meta.run}_permutation-test.tsv
    touch prc_${meta.scenario}_${meta.run}_ttest.tsv

    touch auc_${meta.scenario}_${meta.run}.tsv
    """
}