include { PLOT_FIG_01 } from '../modules/static_plots/plot_fig_01'
include { PLOT_FIG_07 } from '../modules/static_plots/plot_fig_07'
include { PLOT_FIG_S01 } from '../modules/static_plots/plot_fig_s01'

workflow STATIC_PLOTS {
    main:
        ch_plot_fig_01_png = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/main/assets/Fig_01/Fig_01.png')
        ch_plot_fig_01_drawio = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/main/assets/Fig_01/Fig_01.drawio')
        
        PLOT_FIG_01(ch_plot_fig_01_png, ch_plot_fig_01_drawio)


        ch_plot_fig_07_png = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/main/assets/Fig_07/Fig_07.png')
        ch_plot_fig_07_drawio = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/main/assets/Fig_07/Fig_07.drawio')
        
        PLOT_FIG_07(ch_plot_fig_07_png, ch_plot_fig_07_drawio)


        ch_plot_fig_s01 = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/main/assets/Fig_S01/Fig_S01.png')
        
        PLOT_FIG_S01(ch_plot_fig_s01)
}