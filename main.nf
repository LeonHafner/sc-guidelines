include { BENCHMARKING } from './workflows/benchmarking'


workflow {
    BENCHMARKING(
        params.runtime.enabled,
        params.runtime.n_runs,
        params.runtime.n_fixed_cells,
        params.runtime.n_fixed_genes,
        params.runtime.preprocessing_threshold,

        params.static_plots.enabled,
    )
}