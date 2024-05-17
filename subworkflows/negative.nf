include { SIMULATION } from '../modules/performance/simulation/main'
include { PREPROCESSING } from '../modules/performance/preprocessing/main'
include { PSEUDOBULKING } from '../modules/performance/pseudobulking/main'
include { MAST } from '../modules/performance/mast/main'
include { DISTINCT } from '../modules/performance/distinct/main'
include { DESEQ2 } from '../modules/performance/deseq2/main'
include { PERMUTATION_TEST } from '../modules/performance/permutation_test/main'
include { HIERARCHICAL_BOOTSTRAPPING } from '../modules/performance/hierarchical_bootstrapping/main'
include { SCVI } from '../modules/negative/scvi/main'
include { DREAM } from '../modules/performance/dream/main'
include { PVALUES } from '../modules/performance/pvalues/main'
include { SCVI_PROCESSING } from '../modules/negative/scvi_processing/main'
include { FPS_FPR } from '../modules/negative/fps_fpr/main'
include { PLOT_FIG_07 } from '../modules/negative/plot_fig_07/main'
include { PLOT_FIG_S04 } from '../modules/negative/plot_fig_s04/main'
include { PLOT_FIG_S05 } from '../modules/negative/plot_fig_s05/main'


workflow NEGATIVE {
    take:
        kang
        n_genes
        preprocessing_threshold
        

    main:
        ch_kang = Channel.from([[scenario: "kang2018", run:"1"]]).combine(Channel.from(kang))
        ch_meta = Channel.from([[scenario: "atlas-negative", run:"1"]])
        
        SIMULATION(ch_meta, n_genes)

        ch_input_preprocessing = SIMULATION.out.mix(ch_kang)

        PREPROCESSING(ch_input_preprocessing, preprocessing_threshold)

        PSEUDOBULKING(PREPROCESSING.out)

        MAST(PREPROCESSING.out)
        DISTINCT(PREPROCESSING.out)
        DESEQ2(PSEUDOBULKING.out)
        PERMUTATION_TEST(PSEUDOBULKING.out)
        HIERARCHICAL_BOOTSTRAPPING(PREPROCESSING.out)
        SCVI(PREPROCESSING.out)
        DREAM(PSEUDOBULKING.out)

        // Mix all channels (except scvi), introduce grouping criterion 'scenario_run', group them by 'scenario_run' and remove the grouping criterion at the end
        ch_all_methods = MAST.out
            .mix(
                DISTINCT.out,
                DESEQ2.out,
                PERMUTATION_TEST.out,
                HIERARCHICAL_BOOTSTRAPPING.out,
                DREAM.out)
            .map{meta, path -> [meta.scenario + '_' + meta.run, meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}
        
        PVALUES(ch_all_methods)

        SCVI_PROCESSING(SCVI.out)
        
        // Combine all pvalues with scvi pvalues and group them by scenario and run
        ch_input_fpr_fps = PVALUES.out.concat(SCVI_PROCESSING.out)
            .map{meta, path -> [meta.scenario + '_' + meta.run, meta,  path]}
            .groupTuple()
            .map{key, metas, paths -> [metas[0], paths[0], paths[1]]}

        FPS_FPR(ch_input_fpr_fps)

        PLOT_FIG_07(FPS_FPR.out.fpr.filter{meta, path -> meta.scenario == 'atlas-negative'})

        // Channel for kang and simulated scenario - use 'first' for stability in case more runs of these scenarios are added in later versions
        ch_fig_s04_kang = FPS_FPR.out.fpr.filter{meta, path -> meta.scenario == 'kang2018'}.first()
        ch_fig_s04_atlas_negative = FPS_FPR.out.fpr.filter{meta, path -> meta.scenario == 'atlas-negative'}.first()

        PLOT_FIG_S04(
            ch_fig_s04_kang,
            ch_fig_s04_atlas_negative,
            )
        
        // Channel for kang and simulated scenario - use 'first' for stability in case more runs of these scenarios are added in later versions
        ch_fig_s05_kang = FPS_FPR.out.fps.filter{meta, path -> meta.scenario == 'kang2018'}.first()
        ch_fig_s05_atlas_negative = FPS_FPR.out.fps.filter{meta, path -> meta.scenario == 'atlas-negative'}.first()

        PLOT_FIG_S05(
            ch_fig_s05_kang,
            ch_fig_s05_atlas_negative,
        )
}
