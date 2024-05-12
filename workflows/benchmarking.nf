include { RUNTIME } from '../subworkflows/runtime'
include { STATIC_PLOTS } from '../subworkflows/static_plots'
include { CORRELATION } from '../subworkflows/correlation'
include { NEGATIVE } from '../subworkflows/negative'


workflow BENCHMARKING {
    take:
        kang
        
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
        
}