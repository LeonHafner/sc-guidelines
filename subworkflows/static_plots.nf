include { PLOT_FIG_01 } from '../modules/static_plots/plot_fig_01/main'
include { PLOT_FIG_08 } from '../modules/static_plots/plot_fig_08/main'
include { PLOT_FIG_S01 } from '../modules/static_plots/plot_fig_s01/main'

workflow STATIC_PLOTS {
    main:
        PLOT_FIG_01()
        PLOT_FIG_08()
        PLOT_FIG_S01()
}