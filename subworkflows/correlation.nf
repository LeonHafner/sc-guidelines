include { SIMULATION } from '../modules/performance/simulation/main'
include { UNCORRELATED_GENES } from '../modules/correlation/uncorrelated_genes/main'
include { CORRELATIONS } from '../modules/correlation/correlations/main'
include { PLOT_FIG_03 } from '../modules/correlation/plot_fig_03/main'

workflow CORRELATION {
    take:
        n_runs
        n_genes
        
    main:
        ch_meta = Channel.from('dataset').combine(Channel.from(1..n_runs))
            .map{scenario, run -> [scenario: scenario, run: run]}
        
        SIMULATION(ch_meta, n_genes)

        UNCORRELATED_GENES(SIMULATION.out)

        CORRELATIONS(UNCORRELATED_GENES.out)

        PLOT_FIG_03(CORRELATIONS.out.map{meta, tsv -> tsv}.collect())
}