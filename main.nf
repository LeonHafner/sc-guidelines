include { BENCHMARKING } from './workflows/benchmarking'


workflow {
    BENCHMARKING(
        params.kang,
        
        params.runtime.enabled,
        params.runtime.n_runs,
        params.runtime.n_fixed_cells,
        params.runtime.n_fixed_genes,
        params.runtime.preprocessing_threshold,

        params.static_plots.enabled,

        params.correlation.enabled,
        params.correlation.n_runs,
        params.correlation.n_genes,

        params.negative.enabled,
        params.negative.n_genes,
        params.negative.preprocessing_threshold,
    )
}