include { SIMULATION } from '../modules/performance/simulation'
include { UNBALANCE_ATLAS } from '../modules/performance/unbalance_atlas'
include { GROUND_TRUTH } from '../modules/performance/ground_truth'
include { PREPROCESSING } from '../modules/performance/preprocessing'
include { FILTER_HVG } from '../modules/performance/filter_hvg'
include { PSEUDOBULKING } from '../modules/performance/pseudobulking'
include { MAST } from '../modules/performance/mast'
include { DISTINCT } from '../modules/performance/distinct'
include { DESEQ2 } from '../modules/performance/deseq2'
include { PERMUTATION_TEST } from '../modules/performance/permutation_test'
include { HIERARCHICAL_BOOTSTRAPPING } from '../modules/performance/hierarchical_bootstrapping'
include { SCVI } from '../modules/performance/scvi'
include { DREAM } from '../modules/performance/dream'
include { PVALUES } from '../modules/performance/pvalues'
include { PRECISION_RECALL } from '../modules/performance/precision_recall'
include { PREPARE_FIG_02 } from '../modules/performance/prepare_fig_02'
include { PLOT_FIG_02 } from '../modules/performance/plot_fig_02'
include { PLOT_FIG_03 } from '../modules/performance/plot_fig_03'
include { PLOT_FIG_04 } from '../modules/performance/plot_fig_04'
include { PLOT_FIG_S02 } from '../modules/performance/plot_fig_s02'
include { PLOT_FIG_S03 } from '../modules/performance/plot_fig_s03'
include { PLOT_FIG_S04 } from '../modules/performance/plot_fig_s04'
include { PLOT_FIG_S05 } from '../modules/performance/plot_fig_s05'
include { PLOT_FIG_S07 } from '../modules/performance/plot_fig_s07'
include { PLOT_FIG_S09 } from '../modules/performance/plot_fig_s09'

workflow PERFORMANCE {
    take:
        kang
        reactome
        n_runs
        n_genes
        hvg_ratio
        preprocessing_threshold

    main:
        scenarios_regular = ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells']
        scenarios_hvg = ['atlas_hvg', 'dataset_hvg', 'atlas-ub-conditions_hvg', 'dataset-ub-cells_hvg']
        scenarios_less_de = ['atlas-less-de', 'dataset-less-de', 'atlas-ub-conditions-less-de', 'dataset-ub-cells-less-de']

        ch_scenarios = Channel.from(scenarios_regular + scenarios_hvg + scenarios_less_de)
        ch_runs = Channel.from(1..n_runs)
        
        ch_meta = ch_scenarios.combine(ch_runs).map{scenario, run -> [scenario: scenario, run: run]}
        ch_kang = Channel.from([[scenario: 'kang2018', run:'1']]).combine(Channel.from(kang))
        
        SIMULATION(ch_meta, n_genes)
        
        // Split channel as UNBALANCE_ATLAS need only to be executed on the 'atlas-ub-conditions' scenarios
        ch_atlas_ub = SIMULATION.out.filter{item -> ['atlas-ub-conditions', 'atlas-ub-conditions_hvg', 'atlas-ub-conditions-less-de'].contains(item[0].scenario)}
        ch_non_atlas_ub = SIMULATION.out.filter{item -> !['atlas-ub-conditions', 'atlas-ub-conditions_hvg', 'atlas-ub-conditions-less-de'].contains(item[0].scenario)}

        UNBALANCE_ATLAS(ch_atlas_ub)
        
        // Re-combine channels
        ch_preprocessing_input = ch_non_atlas_ub.mix(UNBALANCE_ATLAS.out).mix(ch_kang)

        GROUND_TRUTH(ch_preprocessing_input, reactome)
        
        // Create ground truth objects for fixed effect comparison
        ch_ground_truth_fixed_effect_comparison = GROUND_TRUTH.out
            .filter{item -> item[0].scenario == 'dataset-ub-cells'}
            .map{item -> [[scenario: 'dataset-ub-cells_pb-fixed-effect', run: item[0].run], item[1]]}

        PREPROCESSING(ch_preprocessing_input, preprocessing_threshold)

        ch_hvg = PREPROCESSING.out.filter{item -> ['atlas_hvg', 'dataset_hvg', 'atlas-ub-conditions_hvg', 'dataset-ub-cells_hvg'].contains(item[0].scenario)}
        ch_all_genes = PREPROCESSING.out.filter{item -> !['atlas_hvg', 'dataset_hvg', 'atlas-ub-conditions_hvg', 'dataset-ub-cells_hvg'].contains(item[0].scenario)}

        FILTER_HVG(ch_hvg, hvg_ratio)
        
        ch_methods_input = ch_all_genes.mix(FILTER_HVG.out)

        PSEUDOBULKING(ch_methods_input)

        ch_pseudobulking_fixed_effect_comparison = PSEUDOBULKING.out
            .filter{item -> item[0].scenario == 'dataset-ub-cells'}
            .map{item -> [[scenario: 'dataset-ub-cells_pb-fixed-effect', run: item[0].run], item[1]]}

        MAST(ch_methods_input)
        DISTINCT(ch_methods_input)
        HIERARCHICAL_BOOTSTRAPPING(ch_methods_input)
        SCVI(ch_methods_input)
        DESEQ2(PSEUDOBULKING.out.mix(ch_pseudobulking_fixed_effect_comparison))
        PERMUTATION_TEST(PSEUDOBULKING.out)
        DREAM(PSEUDOBULKING.out.mix(ch_pseudobulking_fixed_effect_comparison))

        // Mix all channels, introduce grouping criterion 'scenario_run', group them by 'scenario_run' and remove the grouping criterion at the end
        ch_all_methods = MAST.out
            .mix(
                DISTINCT.out,
                DESEQ2.out,
                PERMUTATION_TEST.out,
                HIERARCHICAL_BOOTSTRAPPING.out,
                SCVI.out,
                DREAM.out)
            .map{meta, path -> [meta.scenario + '_' + meta.run, meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}

        PVALUES(ch_all_methods)

        // Cartesian product of the channels filtered by matching meta object and mapped to linear tuple [meta, path_pvalues, path_ground_truth]
        // Fixed effect comparison object ground truth is merged into this channel
        ch_meta_pval_groundtruth = PVALUES.out.combine(GROUND_TRUTH.out.mix(ch_ground_truth_fixed_effect_comparison)).filter{it[0] == it[2]}.map{item -> [item[0], item[1], item[3]]}

        PRECISION_RECALL(ch_meta_pval_groundtruth)

        // Filter for scenario 'dataset-ub-cells', remove meta object, collect and add new metadata back
        ch_fig_02 = PREPROCESSING.out
            .filter{meta, path -> meta.scenario == 'dataset-ub-cells'}
            .map{meta, path -> [path]}
            .collect()
            .map{path -> [[scenario: 'dataset-ub-cells'], path]}

        ch_fig_02_png = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/80e1bd085596a121fb490ef924da34fa18f250af/plotting/Fig_02/Fig_2.png')
        ch_fig_02_drawio = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/80e1bd085596a121fb490ef924da34fa18f250af/plotting/Fig_02/Fig_2.drawio')

        PREPARE_FIG_02(ch_fig_02)

        PLOT_FIG_02(PREPARE_FIG_02.out, ch_fig_02_png, ch_fig_02_drawio)

        // Filter for the scenarios needed for Fig_03, introduce pseudokey, group to get list of metas and list of paths, remove pseudokey
        ch_fig_03 = PRECISION_RECALL.out.auc
            .filter{item -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(item[0].scenario)}
            .map{meta, path -> ["key", meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}
        
        PLOT_FIG_03(ch_fig_03)

        ch_fig_04 = PRECISION_RECALL.out.prc
            .filter{meta, path -> meta.scenario == 'kang2018'}
    
        PLOT_FIG_04(ch_fig_04)

        ch_fig_s02 = PRECISION_RECALL.out.auc
            .filter{meta, path -> ['dataset-ub-cells', 'dataset-ub-cells_pb-fixed-effect'].contains(meta.scenario)}
            .map{meta, path -> ["key", meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}
        
        PLOT_FIG_S02(ch_fig_s02)

        // Filter prc files by scenario, map 'run' to first argument to group them by run and map them back attaching a new meta object and flattening the list of paths
        ch_fig_s03 = PRECISION_RECALL.out.prc
            .filter{meta, path -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(meta.scenario)}
            .map{meta, path -> [meta.run, [path]]}
            .groupTuple()
            .map{run, paths -> [[run: run], paths.flatten()]}

        PLOT_FIG_S03(ch_fig_s03)

        // Filter for the scenarios needed for Fig_S06, introduce pseudokey, group to get list of metas and list of paths, remove pseudokey
        ch_fig_s04 = PRECISION_RECALL.out.auc
            .filter{meta, path -> [
                'atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells',
                'atlas-less-de', 'dataset-less-de', 'atlas-ub-conditions-less-de', 'dataset-ub-cells-less-de'
                ].contains(meta.scenario)}
            .map{meta, path -> ["key", meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}

        PLOT_FIG_S04(ch_fig_s04)

        // Filter for the scenarios needed for Fig_S02, introduce pseudokey, group to get list of metas and list of paths, remove pseudokey
        ch_fig_s05 = PRECISION_RECALL.out.auc
            .filter{item -> ['atlas_hvg', 'dataset_hvg', 'atlas-ub-conditions_hvg', 'dataset-ub-cells_hvg'].contains(item[0].scenario)}
            .map{meta, path -> ["key", meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}
    
        PLOT_FIG_S05(ch_fig_s05)

        ch_fig_s07_prc = PRECISION_RECALL.out.prc
            .filter{meta, path -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(meta.scenario)}
            .map{meta, path -> [meta.run, [path]]}
            .groupTuple()
            .map{run, paths -> [[run: run], paths.flatten()]}
        
        // Filter for scenario, introduce pseudokey, group to get list of meta and list of paths, remove pseudokey, collect to get value channel
        ch_fig_s07_auc = PRECISION_RECALL.out.auc
            .filter{meta, path -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(meta.scenario)}
            .map{meta, path -> ["key", meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}
            .collect()

        PLOT_FIG_S07(ch_fig_s07_prc, ch_fig_s07_auc)

        ch_fig_s09 = PSEUDOBULKING.out.filter{meta, file -> meta.scenario == "atlas-ub-conditions"}
        
        PLOT_FIG_S09(ch_fig_s09)
}