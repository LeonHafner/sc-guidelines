include { RUNTIME } from '../subworkflows/runtime'


workflow BENCHMARKING {
    take:
        runtime_enabled
        runtime_n_runs
        runtime_n_fixed_cells
        runtime_n_fixed_genes
        runtime_preprocessing_threshold

    main:
        if (runtime_enabled) {
            RUNTIME(
                runtime_n_runs,
                runtime_n_fixed_cells,
                runtime_n_fixed_genes,
                runtime_preprocessing_threshold,
            )
        }
        
}