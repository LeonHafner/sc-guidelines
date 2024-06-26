include { SIMULATION } from '../modules/performance/simulation'
include { UNCORRELATED_GENES } from '../modules/correlation/uncorrelated_genes'
include { CORRELATIONS } from '../modules/correlation/correlations'
include { PLOT_FIG_S09 } from '../modules/correlation/plot_fig_s09'

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

        PLOT_FIG_S09(CORRELATIONS.out.map{meta, tsv -> tsv}.collect())
}