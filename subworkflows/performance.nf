include { SIMULATION } from '../modules/performance/simulation/main'
include { UNBALANCE_ATLAS } from '../modules/performance/unbalance_atlas/main'
include { GROUND_TRUTH } from '../modules/performance/ground_truth/main'
include { PREPROCESSING } from '../modules/performance/preprocessing/main'
include { FILTER_HVG } from '../modules/performance/filter_hvg/main'
include { PSEUDOBULKING } from '../modules/performance/pseudobulking/main'
include { MAST } from '../modules/performance/mast/main'
include { DISTINCT } from '../modules/performance/distinct/main'
include { DESEQ2 } from '../modules/performance/deseq2/main'
include { PERMUTATION_TEST } from '../modules/performance/permutation_test/main'
include { HIERARCHICAL_BOOTSTRAPPING } from '../modules/performance/hierarchical_bootstrapping/main'
include { SCVI } from '../modules/performance/scvi/main'
include { DREAM } from '../modules/performance/dream/main'
include { PVALUES } from '../modules/performance/pvalues/main'
include { PRECISION_RECALL } from '../modules/performance/precision_recall/main'
include { PREPARE_FIG_02 } from '../modules/performance/prepare_fig_02/main'
include { PLOT_FIG_02 } from '../modules/performance/plot_fig_02/main'
include { PLOT_FIG_05 } from '../modules/performance/plot_fig_05/main'
include { PLOT_FIG_06 } from '../modules/performance/plot_fig_06/main'
include { PLOT_FIG_S02 } from '../modules/performance/plot_fig_s02/main'
include { PLOT_FIG_S03 } from '../modules/performance/plot_fig_s03/main'
include { PLOT_FIG_S06 } from '../modules/performance/plot_fig_s06/main'
include { PLOT_FIG_S07 } from '../modules/performance/plot_fig_s07/main'
include { PLOT_FIG_S08 } from '../modules/performance/plot_fig_s08/main'

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
            .map{item -> [item[0].scenario + '_' + item[0].run, [item[0], item[1]]]}
            .groupTuple()
            .map{item -> item[1]}
        
        PVALUES(ch_all_methods)

        // Cartesian product of the channels filtered by matching meta object and mapped to linear tuple [meta, path_pvalues, path_ground_truth]
        // Fixed effect comparison object ground truth is merged into this channel
        ch_meta_pval_groundtruth = PVALUES.out.combine(GROUND_TRUTH.out.mix(ch_ground_truth_fixed_effect_comparison)).filter{it[0] == it[2]}.map{item -> [item[0], item[1], item[3]]}

        PRECISION_RECALL(ch_meta_pval_groundtruth)

        // Filter for scenario 'dataset-ub-cells', slice of meta object and enclose into tuple, collect and add single metadata back
        ch_fig_02 = PREPROCESSING.out
            .filter{item -> item[0].scenario == 'dataset-ub-cells'}
            .map{item -> [item[1]]}
            .collect()
            .map{item -> [[scenario: 'dataset-ub-cells'], item]}

        PREPARE_FIG_02(ch_fig_02)

        PLOT_FIG_02(PREPARE_FIG_02.out)

        // Filter for the scenarios needed for Fig_07, enclose the items in single value tuples and collect them
        ch_fig_05 = PRECISION_RECALL.out.auc
            .filter{item -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(item[0].scenario)}
            .map{item -> [item]}
            .collect()
        
        PLOT_FIG_05(ch_fig_05)

        ch_fig_06 = PRECISION_RECALL.out.prc
            .filter{item -> item[0].scenario == 'kang2018'}
    
        PLOT_FIG_06(ch_fig_06)

        ch_fig_s02 = PRECISION_RECALL.out.auc
            .filter{item -> ['atlas_hvg', 'dataset_hvg', 'atlas-ub-conditions_hvg', 'dataset-ub-cells_hvg'].contains(item[0].scenario)}
            .map{item -> [item]}
            .collect()
    
        PLOT_FIG_S02(ch_fig_s02)

        // Filter prc files by scenario, map 'run' to first argument to group them by run and map them back attaching a new meta object and flattening the list of paths
        ch_fig_s03 = PRECISION_RECALL.out.prc
            .filter{item -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(item[0].scenario)}
            .map{item -> [item[0].run, [item[1]]]}
            .groupTuple()
            .map{item -> [[run: item[0]], item[1].flatten()]}

        PLOT_FIG_S03(ch_fig_s03)

        ch_fig_s06 = PRECISION_RECALL.out.auc
            .filter{item -> [
                'atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells',
                'atlas-less-de', 'dataset-less-de', 'atlas-ub-conditions-less-de', 'dataset-ub-cells-less-de'
                ].contains(item[0].scenario)}
            .map{item -> [item]}
            .collect()

        PLOT_FIG_S06(ch_fig_s06)

        ch_fig_s07_prc = PRECISION_RECALL.out.prc
            .filter{item -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(item[0].scenario)}
            .map{item -> [item[0].run, [item[1]]]}
            .groupTuple()
            .map{item -> [[run: item[0]], item[1].flatten()]}
        
        ch_fig_s07_auc = PRECISION_RECALL.out.auc
            .filter{item -> ['atlas', 'dataset', 'atlas-ub-conditions', 'dataset-ub-cells'].contains(item[0].scenario)}
            .map{item -> [item]}
            .collect()

        PLOT_FIG_S07(ch_fig_s07_prc, ch_fig_s07_auc)

        ch_fig_s08 = PRECISION_RECALL.out.auc
            .filter{item -> ['dataset-ub-cells', 'dataset-ub-cells_pb-fixed-effect'].contains(item[0].scenario)}
            .map{item -> [item]}
            .collect()
        
        PLOT_FIG_S08(ch_fig_s08)
    }