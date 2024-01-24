params.output = "/nfs/home/students/l.hafner/fopra/nextflow_test/runner/output"
params.scenarios = [
    "atlas", "dataset", "atlas-ub-conditions", "dataset-ub-cells",
    "atlas_hvg", "dataset_hvg", "atlas-ub-conditions_hvg", "dataset-ub-cells_hvg",
    "atlas-less-de", "dataset-less-de", "atlas-ub-conditions-less-de", "dataset-ub-cells-less-de"
]
params.runs = 20
params.n_genes = 5000
params.hvg_ratio = 0.1
params.preprocessing_threshold = 0.1

process SIMULATION {
    tag "${meta.scenario}_${meta.run}"
    
    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    simulation.R \
        --scenario ${meta.scenario} \
        --n_genes ${params.n_genes} \
        --output ${meta.scenario}_${meta.run}.h5ad
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}

process UNBALANCE_ATLAS {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra_py'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    unbalance_atlas.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_${meta.run}.h5ad
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}

process FILTER_HIGHLY_VARIABLE_GENES {
    tag "${meta.scenario}_${meta.run}"
    
    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra_py'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    #!/usr/bin/env python3

    import scanpy as sc

    ad = sc.read_h5ad("${input_anndata}")

    sc.pp.highly_variable_genes(
        adata=ad,
        n_top_genes=${params.n_genes} * ${params.hvg_ratio},
        flavor='seurat_v3',
        subset=True,
    )

    ad.write_h5ad("${meta.scenario}_${meta.run}.h5ad")
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}

process PREPROCESSING {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra_py'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')
    val threshold

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    preprocessing.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_${meta.run}.h5ad \
        --threshold ${threshold}
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}

process PSEUDBULKING {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/py311'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata')

    output:
    tuple val(meta), path("${meta.scenario}_${meta.run}.h5ad")

    script:
    """
    pseudobulk.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_${meta.run}.h5ad \
        --scenario ${meta.scenario}
    """

    stub:
    """
    touch ${meta.scenario}_${meta.run}.h5ad
    """
}

process GET_GROUND_TRUTH {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra_py'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("de_${meta.scenario}_${meta.run}.tsv")

    script:
    """
    ground_truth.py \
        --input ${input_anndata} \
        --output de_${meta.scenario}_${meta.run}.tsv
    """

    stub:
    """
    touch de_${meta.scenario}_${meta.run}.tsv
    """
}

process RUN_MAST {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_mast_${meta.run}.tsv")

    script:
    meta = meta + [method: 'mast']
    """
    mast.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_mast_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'mast']
    """
    touch ${meta.scenario}_mast_${meta.run}.tsv
    """
}

process RUN_DISTINCT {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'
    
    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_distinct_${meta.run}.tsv")

    script:
    meta = meta + [method: 'distinct']
    """
    distinct.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_distinct_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'distinct']
    """
    touch ${meta.scenario}_distinct_${meta.run}.tsv
    """
}

process RUN_DESEQ2 {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_deseq2_${meta.run}.tsv")

    script:
    meta = meta + [method: 'deseq2']
    """
    deseq2.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_deseq2_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'deseq2']
    """
    touch ${meta.scenario}_deseq2_${meta.run}.tsv
    """
}

process RUN_PERMUTATION_TEST {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra_py'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_permutation-test_${meta.run}.tsv")

    script:
    meta = meta + [method: 'permutation-test']
    """
    permutation.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_permutation-test_${meta.run}.tsv \
        -n 10000 \
        -n_max 100000 \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'permutation-test']
    """
    touch ${meta.scenario}_permutation-test_${meta.run}.tsv
    """
}

process RUN_HIERARCHICAL_BOOTSTRAPPING {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra_py'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_hierarchical-bootstrapping_${meta.run}.tsv")

    script:
    meta = meta + [method: 'hierarchical-bootstrapping']
    """
    hierarchical_bootstrapping.py \
        --input ${input_anndata} \
        --output1 ${meta.scenario}_hierarchical-bootstrapping_${meta.run}_pre.tsv \
        --output2 ${meta.scenario}_hierarchical-bootstrapping_${meta.run}.tsv \
        --scenario ${meta.scenario} 
    """

    stub:
    meta = meta + [method: 'hierarchical-bootstrapping']
    """
    touch ${meta.scenario}_hierarchical-bootstrapping_${meta.run}.tsv
    """
}

process RUN_SCVI {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/scvi_latest'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_scvi_${meta.run}.tsv")

    script:
    meta = meta + [method: 'scvi']
    """
    scvi-de.py \
        --input ${input_anndata} \
        --output ${meta.scenario}_scvi_${meta.run}.tsv \
        --scenario ${meta.scenario}
    """

    stub:
    meta = meta + [method: 'scvi']
    """
    touch ${meta.scenario}_scvi_${meta.run}.tsv
    """
}

process RUN_DREAM {
    tag "${meta.scenario}_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    input:
    tuple val(meta), path(input_anndata, stageAs: 'input_anndata.h5ad')

    output:
    tuple val(meta), path("${meta.scenario}_dream_${meta.run}.tsv")

    script:
    meta = meta + [method: 'dream']
    """
    dream.R \
        --input ${input_anndata} \
        --output ${meta.scenario}_dream_${meta.run}.tsv \
        --scenario ${meta.scenario} \
        --threads 1
    """

    stub:
    meta = meta + [method: 'dream']
    """
    touch ${meta.scenario}_dream_${meta.run}.tsv
    """
}

process GET_PVALUES {
    tag "${tuple_list[0][0].scenario}_${tuple_list[0][0].run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/py311'

    input:
    val tuple_list

    output:
    tuple val(meta), path("pval_${meta.scenario}_${meta.run}.tsv")

    script:
    // Process channel input which is a list of tuples consisting of a meta map and a path object
    def meta_list = []
    def path_list = []
    for (element in tuple_list) {
        meta_list.add(element[0])
        path_list.add(element[1])
    }
    def meta_string = meta_list.join(';')
    def path_string = path_list.join(';')
    meta = [scenario: tuple_list[0][0].scenario, run: tuple_list[0][0].run]
    """
    pvalues.py \
        --meta_string "${meta_string}" \
        --path_string "${path_string}" \
        --output pval_${meta.scenario}_${meta.run}.tsv
    """

    stub:
    meta = [scenario: tuple_list[0][0].scenario, run: tuple_list[0][0].run]
    """
    touch pval_${meta.scenario}_${meta.run}.tsv
    """
}

process GET_PRC {
    tag "${meta.scenario}_${meta.run}"
    
    conda '/nfs/home/students/l.hafner/miniconda3/envs/py311'

    input:
    tuple val(meta), path(pvalues, stageAs: 'pvalues.tsv'), path(ground_truth, stageAs: 'ground_truth.tsv')

    output:
    tuple val(meta), path("prc_${meta.scenario}_${meta.run}_*.tsv", arity: 7), emit: 'prc'
    tuple val(meta), path("auc_${meta.scenario}_${meta.run}.tsv"), emit: 'auc'

    script:
    """
    precision-recall.py \
        --pvalues ${pvalues} \
        --ground_truth ${ground_truth} \
        --output ${meta.scenario}_${meta.run}
    """

    stub:
    """
    touch prc_${meta.scenario}_${meta.run}_deseq2.tsv
    touch prc_${meta.scenario}_${meta.run}_mast.tsv
    touch prc_${meta.scenario}_${meta.run}_scvi.tsv
    touch prc_${meta.scenario}_${meta.run}_dream.tsv
    touch prc_${meta.scenario}_${meta.run}_distinct.tsv
    touch prc_${meta.scenario}_${meta.run}_hierarchical-bootstrapping.tsv
    touch prc_${meta.scenario}_${meta.run}_permutation-test.tsv

    touch auc_${meta.scenario}_${meta.run}.tsv
    """

}

process PREPARE_COUNTS_FIG_03 {
    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra_py'

    input:
    tuple val(meta), val(path_list)

    output:
    tuple val(meta), path('cell_counts.txt')

    script:
    """
    Fig_03_prepare_counts.py \
        --input "${path_list}" \
        --output "cell_counts.txt"
    """

    stub:
    """
    touch cell_counts.txt
    """
}

process PLOT_FIG_03 {
    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(meta), path(cell_counts)

    output:
    path 'Fig_03.png'

    script:
    """
    Fig_03.R \
        --input ${cell_counts} \
        --output Fig_03.png
    """

    stub:
    """
    touch Fig_03.png
    """
}

process PLOT_FIG_07 {
    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    publishDir "${params.output}/Fig_07", mode: 'copy'

    input:
    val tuple_list

    output:
    path 'Fig_07.png'

    script:
    // Process channel input which is a list of tuples consisting of a meta map and a path object
    def meta_list = []
    def path_list = []
    for (element in tuple_list) {
        meta_list.add(element[0])
        path_list.add(element[1])
    }
    def meta_string = meta_list.join(';')
    def path_string = path_list.join(';')
    """
    Fig_07.R \
        --meta_string "${meta_string}" \
        --path_string "${path_string}" \
        --output Fig_07.png
    """

    stub:
    """
    touch Fig_07.png
    """
}

process PLOT_FIG_07_HVG {
    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    publishDir "${params.output}/Fig_07", mode: 'copy'

    input:
    val tuple_list

    output:
    path 'Fig_07_HVG.png'

    script:
    // Process channel input which is a list of tuples consisting of a meta map and a path object
    def meta_list = []
    def path_list = []
    for (element in tuple_list) {
        meta_list.add(element[0])
        path_list.add(element[1])
    }
    def meta_string = meta_list.join(';')
    def path_string = path_list.join(';')
    """
    Fig_07_HVG.R \
        --meta_string "${meta_string}" \
        --path_string "${path_string}" \
        --output Fig_07_HVG.png
    """

    stub:
    """
    touch Fig_07_HVG.png
    """
}

process PLOT_FIG_11 {
    tag "run_${meta.run}"

    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    publishDir "${params.output}/Fig_11", mode: 'copy'

    input:
    tuple val(meta), val(path_list)

    output:
    path "Fig_11_run${meta.run}.png"

    script:
    def path_string = path_list.join(";")
    """
    Fig_11.R \
        --path_string "${path_string}" \
        --output Fig_11_run${meta.run}.png
    """

    stub:
    def path_string = path_list.join(";")
    """
    touch Fig_11_run${meta.run}.png
    """
}

process PLOT_FIG_14 {
    conda '/nfs/home/students/l.hafner/miniconda3/envs/fopra'

    publishDir "${params.output}", mode: 'copy'

    input:
    val tuple_list

    output:
    path 'Fig_14.png'

    script:
    def meta_list = []
    def path_list = []
    for (element in tuple_list) {
        meta_list.add(element[0])
        path_list.add(element[1])
    }
    def meta_string = meta_list.join(';')
    def path_string = path_list.join(';')
    """
    Fig_14.R \
        --meta_string "${meta_string}" \
        --path_string "${path_string}" \
        --output Fig_14.png
    """

    stub:
    """
    touch Fig_14.png
    """
}


workflow {
    ch_scenarios = Channel.from(params.scenarios)
    ch_runs = Channel.from(1..params.runs)
    
    ch_meta = ch_scenarios.combine(ch_runs).map{scenario, run -> [scenario: scenario, run: run]}

    SIMULATION(ch_meta)
    
    // Split channel as UNBALANCE_ATLAS need only to be executed on the 'atlas-ub-conditions' scenarios
    ch_atlas_ub = SIMULATION.out.filter{item -> ['atlas-ub-conditions', 'atlas-ub-conditions_hvg', 'atlas-ub-conditions-less-de'].contains(item[0].scenario)}
    ch_non_atlas_ub = SIMULATION.out.filter{item -> !['atlas-ub-conditions', 'atlas-ub-conditions_hvg', 'atlas-ub-conditions-less-de'].contains(item[0].scenario)}
    
    // Execute unbalancing only for 'atlas-ub-conditions' scenario
    UNBALANCE_ATLAS(ch_atlas_ub)

    // Re-combine channels
    ch_preprocessing_input = ch_non_atlas_ub.mix(UNBALANCE_ATLAS.out)
    
    GET_GROUND_TRUTH(ch_preprocessing_input)
    
    PREPROCESSING(ch_preprocessing_input, params.preprocessing_threshold)
    
    ch_hvg = PREPROCESSING.out.filter{item -> ["atlas_hvg", "dataset_hvg", "atlas-ub-conditions_hvg", "dataset-ub-cells_hvg"].contains(item[0].scenario)}
    ch_all_genes = PREPROCESSING.out.filter{item -> !["atlas_hvg", "dataset_hvg", "atlas-ub-conditions_hvg", "dataset-ub-cells_hvg"].contains(item[0].scenario)}
    
    FILTER_HIGHLY_VARIABLE_GENES(ch_hvg)

    ch_methods_input = ch_all_genes.mix(FILTER_HIGHLY_VARIABLE_GENES.out)
    
    PSEUDBULKING(ch_methods_input)
    
    RUN_MAST(ch_methods_input)
    RUN_DISTINCT(ch_methods_input)
    RUN_DESEQ2(PSEUDBULKING.out)
    RUN_PERMUTATION_TEST(PSEUDBULKING.out)
    RUN_HIERARCHICAL_BOOTSTRAPPING(ch_methods_input)
    RUN_SCVI(ch_methods_input)
    RUN_DREAM(PSEUDBULKING.out)
    
    // Mix all channels, introduce grouping criterion 'scenario_run', group them by 'scenario_run' and remove the grouping criterion at the end
    ch_all_methods = RUN_MAST.out.mix(RUN_DISTINCT.out, RUN_DESEQ2.out, RUN_PERMUTATION_TEST.out, RUN_HIERARCHICAL_BOOTSTRAPPING.out, RUN_SCVI.out, RUN_DREAM.out)
        .map{
            item -> [item[0].scenario + '_' + item[0].run, [item[0], item[1]]]}
        .groupTuple()
        .map{
            item -> item[1]
        }
    
    GET_PVALUES(ch_all_methods)

    // Cartesian product of the channels filtered by matching meta object and mapped to linear tuple [meta, path_pvalues, path_ground_truth]
    ch_meta_pval_degenes = GET_PVALUES.out.combine(GET_GROUND_TRUTH.out).filter{it[0] == it[2]}.map{item -> [item[0], item[1], item[3]]}
    
    GET_PRC(ch_meta_pval_degenes)
    
    // Filter for scenario 'dataset-ub-cells', slice of meta object and enclose into tuple, collect and add single metadata back
    ch_fig_03 = PREPROCESSING.out
        .filter{item -> item[0].scenario == 'dataset-ub-cells'}
        .map{item -> [item[1]]}
        .collect()
        .map{item -> [[scenario: 'dataset-ub-cells'], item]}
    
    PREPARE_COUNTS_FIG_03(ch_fig_03)

    PLOT_FIG_03(PREPARE_COUNTS_FIG_03.out)
    
    // Filter for the scenarios needed for Fig_07, enclose the items in single value tuples and collect them
    ch_fig_07 = GET_PRC.out.auc
        .filter{item -> ["atlas", "dataset", "atlas-ub-conditions", "dataset-ub-cells"].contains(item[0].scenario)}
        .map{item -> [item]}
        .collect()
    
    PLOT_FIG_07(ch_fig_07)

    // Filter for the scenarios needed for scenarios filtered for highly variable genes, enclose the items in single value tuples and collect them
    ch_fig_07_hvg = GET_PRC.out.auc
        .filter{item -> ["atlas_hvg", "dataset_hvg", "atlas-ub-conditions_hvg", "dataset-ub-cells_hvg"].contains(item[0].scenario)}
        .map{item -> [item]}
        .collect()
    
    PLOT_FIG_07_HVG(ch_fig_07_hvg)

    
    // Filter prc files by scenario, map "run" to first argument to group them by run and map them back attaching a new meta object and flattening the list of paths
    ch_fig_11 = GET_PRC.out.prc
        .filter{item -> ["atlas", "dataset", "atlas-ub-conditions", "dataset-ub-cells"].contains(item[0].scenario)}
        .map{item -> [item[0].run, [item[1]]]}
        .groupTuple()
        .map{item -> [[run: item[0]], item[1].flatten()]}

    PLOT_FIG_11(ch_fig_11)

    ch_fig_14 = GET_PRC.out.auc
        .filter{item -> [
            "atlas", "dataset", "atlas-ub-conditions", "dataset-ub-cells",
            "atlas-less-de", "dataset-less-de", "atlas-ub-conditions-less-de", "dataset-ub-cells-less-de"
            ].contains(item[0].scenario)}
        .map{item -> [item]}
        .collect()

    PLOT_FIG_14(ch_fig_14)

}