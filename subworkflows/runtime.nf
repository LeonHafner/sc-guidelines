include { addRuntimeToMeta } from '../subworkflows/utils'
include { SIMULATION } from '../modules/runtime/simulation/main'
include { PREPROCESSING } from '../modules/runtime/preprocessing/main'
include { PSEUDOBULKING } from '../modules/runtime/pseudobulking/main'
include { MAST } from '../modules/runtime/mast/main'
include { DISTINCT } from '../modules/runtime/distinct/main'
include { DESEQ2 } from '../modules/runtime/deseq2/main'
include { PERMUTATION_TEST } from '../modules/runtime/permutation_test/main'
include { HIERARCHICAL_BOOTSTRAPPING } from '../modules/runtime/hierarchical_bootstrapping/main'
include { SCVI } from '../modules/runtime/scvi/main'
include { DREAM } from '../modules/runtime/dream/main'
include { PREPARE_FIG_04 } from '../modules/runtime/prepare_fig_04/main'
include { PLOT_FIG_04 } from '../modules/runtime/plot_fig_04/main'


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

        SIMULATION(meta)

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
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'pseudobulking')
            }
                
        MAST(ch_out_pseudobulking)

        ch_out_mast = MAST.out
            .map{
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'mast')
            }
        
        DISTINCT(ch_out_mast)
    
        ch_out_distinct = DISTINCT.out
            .map{
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'distinct')
            }

        DESEQ2(ch_out_distinct)

        ch_out_deseq2 = DESEQ2.out
            .map{
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'deseq2')
            }
        
        PERMUTATION_TEST(ch_out_deseq2)

        ch_out_permutation_test = PERMUTATION_TEST.out
            .map{
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'permutation_test')
            }
        
        HIERARCHICAL_BOOTSTRAPPING(ch_out_permutation_test)

        ch_out_hierarchical_bootstrapping = HIERARCHICAL_BOOTSTRAPPING.out
            .map{
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'hierarchical_bootstrapping')
            }

        SCVI(ch_out_hierarchical_bootstrapping)

        ch_out_scvi = SCVI.out
            .map{
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'scvi')
            }

        DREAM(ch_out_scvi)

        ch_out_dream = DREAM.out
            .map{
                meta, sc_anndata, pb_anndata, time_file ->
                addRuntimeToMeta(meta, sc_anndata, pb_anndata, time_file, 'dream')
            }
            .map{meta, sc_anndata, pb_anndata -> meta}
                
        PREPARE_FIG_04(ch_out_dream.collect())

        PLOT_FIG_04(PREPARE_FIG_04.out.fixed_cells, PREPARE_FIG_04.out.fixed_genes)
}