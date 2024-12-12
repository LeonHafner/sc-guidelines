include { addRuntimeToMeta } from '../subworkflows/utils'
include { SIMULATION } from '../modules/runtime/simulation'
include { PREPROCESSING } from '../modules/runtime/preprocessing'
include { PSEUDOBULKING } from '../modules/runtime/pseudobulking'
include { MAST } from '../modules/runtime/mast'
include { DISTINCT } from '../modules/runtime/distinct'
include { DESEQ2 } from '../modules/runtime/deseq2'
include { PERMUTATION_TEST } from '../modules/runtime/permutation_test'
include { HIERARCHICAL_BOOTSTRAPPING } from '../modules/runtime/hierarchical_bootstrapping'
include { SCVI } from '../modules/runtime/scvi'
include { DREAM } from '../modules/runtime/dream'
include { TTEST } from '../modules/runtime/ttest'
include { PREPARE_FIG_06 } from '../modules/runtime/prepare_fig_06'
include { PLOT_FIG_06 } from '../modules/runtime/plot_fig_06'


workflow RUNTIME {
    take:
        n_runs
        n_fixed_cells
        n_fixed_genes
        preprocessing_threshold

    main:
        // Channel for fixed cells simulations
        ch_fixed_cells = Channel.from([scenario: 'fixed_cells', n_cells: n_fixed_cells])
                .combine(1..n_runs)
                .combine([(100..900).step(100), (1000..10000).step(500)].flatten())
                .map{meta, run, n_genes -> meta + [n_genes: n_genes, run: run]}

        // Channel for fixed genes simulations
        ch_fixed_genes = Channel.from([scenario: 'fixed_genes', n_genes: n_fixed_genes])
                .combine(1..n_runs)
                .combine([(100..900).step(100), (1000..10000).step(500)].flatten())
                .map{meta, run, n_cells -> meta + [n_cells: n_cells, run: run]}

        meta = ch_fixed_cells.mix(ch_fixed_genes)

        SIMULATION(meta, preprocessing_threshold)

        ch_out_simulation = SIMULATION.out
            .map{
                meta, anndata, time_file -> 
                addRuntimeToMeta(meta, anndata, time_file, 'simulation')
            }
                
        PREPROCESSING(ch_out_simulation, preprocessing_threshold)

        ch_out_preprocessing = PREPROCESSING.out
            .map{
                meta, anndata, time_file ->
                addRuntimeToMeta(meta, anndata, time_file, 'preprocessing')
            }
                
        PSEUDOBULKING(ch_out_preprocessing)

        ch_out_pseudobulking = PSEUDOBULKING.out
            .map{
                meta, anndata, time_file ->
                addRuntimeToMeta(meta, anndata, time_file, 'pseudobulking')
            }

        MAST(ch_out_preprocessing)
        DISTINCT(ch_out_preprocessing)
        DESEQ2(ch_out_pseudobulking)
        PERMUTATION_TEST(ch_out_pseudobulking)
        HIERARCHICAL_BOOTSTRAPPING(ch_out_preprocessing)
        SCVI(ch_out_preprocessing)
        DREAM(ch_out_pseudobulking)
        TTEST(ch_out_preprocessing)

        ch_out_mast = MAST.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'mast')}
        ch_out_distinct = DISTINCT.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'distinct')}
        ch_out_deseq2 = DESEQ2.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'deseq2')}
        ch_out_permutation_test = PERMUTATION_TEST.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'permutation_test')}
        ch_out_hierarchical_bootstrapping = HIERARCHICAL_BOOTSTRAPPING.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'hierarchical_bootstrapping')}
        ch_out_scvi = SCVI.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'scvi')}
        ch_out_dream = DREAM.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'dream')}
        ch_out_ttest = TTEST.out.map{meta, time_file -> addRuntimeToMeta(meta, time_file, 'ttest')}
        
        // Combine all channels, introduce grouping key, group by key and merge maps together
        ch_combined = ch_out_mast.mix(ch_out_distinct, ch_out_deseq2, ch_out_permutation_test, ch_out_hierarchical_bootstrapping, ch_out_scvi, ch_out_dream, ch_out_ttest)
            .map{meta -> [(meta.scenario + meta.n_cells + meta.n_genes + meta.run).join('_'), meta]}
            .groupTuple()
            .map{key, meta -> meta.flatten().sum()}

        PREPARE_FIG_06(ch_combined.collect())

        PLOT_FIG_06(PREPARE_FIG_06.out.fixed_cells, PREPARE_FIG_06.out.fixed_genes)
}