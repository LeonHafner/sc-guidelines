include { RUNTIME } from '../subworkflows/runtime'
include { STATIC_PLOTS } from '../subworkflows/static_plots'


workflow BENCHMARKING {
    take:
        runtime_enabled
        runtime_n_runs
        runtime_n_fixed_cells
        runtime_n_fixed_genes
        runtime_preprocessing_threshold

        static_plots_enabled

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
        
}