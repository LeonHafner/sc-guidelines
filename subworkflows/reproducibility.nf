include { PREPARE_FILES } from '../modules/reproducibility/prepare_files'
include { PREPROCESSING } from '../modules/performance/preprocessing'
include { PSEUDOBULKING } from '../modules/performance/pseudobulking'
include { MAST } from '../modules/performance/mast'
include { DISTINCT } from '../modules/performance/distinct'
include { DESEQ2 } from '../modules/performance/deseq2'
include { PERMUTATION_TEST } from '../modules/performance/permutation_test'
include { HIERARCHICAL_BOOTSTRAPPING } from '../modules/performance/hierarchical_bootstrapping'
include { SCVI } from '../modules/performance/scvi'
include { DREAM } from '../modules/performance/dream'
include { SCDD } from '../modules/performance/scdd'
include { TTEST } from '../modules/performance/ttest'
include { PVALUES } from '../modules/performance/pvalues'
include { PLOT_FIG_S10 } from '../modules/reproducibility/plot_fig_s10'
include { PLOT_FIG_S11 } from '../modules/reproducibility/plot_fig_s11'
include { PLOT_FIG_S12 } from '../modules/reproducibility/plot_fig_s12'

workflow REPRODUCIBILITY {
    take:
        lung_cancer_atlas
        preprocessing_threshold
    
    main:
        PREPARE_FILES(lung_cancer_atlas)

        ch_preprocessing_input = PREPARE_FILES.out
            .flatten()
            .map{file -> [[run: file.name.split('_')[0], scenario: "luca"], file]}

        PREPROCESSING(ch_preprocessing_input, preprocessing_threshold)
        
        PSEUDOBULKING(PREPROCESSING.out)
        
        MAST(PREPROCESSING.out)
        DISTINCT(PREPROCESSING.out)
        DESEQ2(PSEUDOBULKING.out)
        PERMUTATION_TEST(PSEUDOBULKING.out)
        HIERARCHICAL_BOOTSTRAPPING(PREPROCESSING.out)
        SCVI(PREPROCESSING.out)
        DREAM(PSEUDOBULKING.out)
        SCDD(PREPROCESSING.out)
        TTEST(PREPROCESSING.out)

        ch_all_methods = MAST.out
            .mix(
                DISTINCT.out,
                DESEQ2.out,
                PERMUTATION_TEST.out,
                HIERARCHICAL_BOOTSTRAPPING.out,
                SCVI.out,
                DREAM.out,
                SCDD.out,
                TTEST.out)
            .map{meta, path -> [meta.scenario + '_' + meta.run, meta, path]}
            .groupTuple()
            .map{key, meta, path -> [meta, path]}
        
        PVALUES(ch_all_methods)

        ch_figs = PVALUES.out
            .map{meta, file -> file}
            .collect()

        PLOT_FIG_S10(ch_figs)

        PLOT_FIG_S11(ch_figs)

        PLOT_FIG_S12(ch_figs)
}