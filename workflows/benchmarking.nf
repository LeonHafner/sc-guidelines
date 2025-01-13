include { RUNTIME } from '../subworkflows/runtime'
include { STATIC_PLOTS } from '../subworkflows/static_plots'
include { CORRELATION } from '../subworkflows/correlation'
include { NEGATIVE } from '../subworkflows/negative'
include { PERFORMANCE } from '../subworkflows/performance'
include { REPRODUCIBILITY } from '../subworkflows/reproducibility'


workflow BENCHMARKING {
    take:
        kang
        reactome
        
        runtime_enabled
        runtime_n_runs
        runtime_n_fixed_cells
        runtime_n_fixed_genes
        runtime_preprocessing_threshold

        static_plots_enabled

        correlation_enabled
        correlation_n_runs
        correlation_n_genes

        negative_enabled
        negative_n_genes
        negative_preprocessing_threshold

        performance_enabled
        performance_n_runs
        performance_n_genes
        performance_hvg_ratio
        performance_preprocessing_threshold

        reproducibility_enabled
        reproducibility_lung_cancer_atlas
        reproducibility_preprocessing_threshold
        
    main:
        if (runtime_enabled) {
            RUNTIME(
                runtime_n_runs,
                runtime_n_fixed_cells,
                runtime_n_fixed_genes,
                runtime_preprocessing_threshold,
            )
        }
        
        if (static_plots_enabled) {
            STATIC_PLOTS()
        }

        if (correlation_enabled) {
            CORRELATION(
                correlation_n_runs,
                correlation_n_genes,
            )
        }

        if (negative_enabled) {
            NEGATIVE(
                kang,
                negative_n_genes,
                negative_preprocessing_threshold,
            )
        }

        if (performance_enabled) {
            PERFORMANCE(
                kang,
                reactome,
                performance_n_runs,
                performance_n_genes,
                performance_hvg_ratio,
                performance_preprocessing_threshold,
            )
        }

        if (reproducibility_enabled) {
            REPRODUCIBILITY(
                reproducibility_lung_cancer_atlas,
                reproducibility_preprocessing_threshold
            )
        }
}