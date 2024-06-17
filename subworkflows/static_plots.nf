include { PLOT_FIG_01 } from '../modules/static_plots/plot_fig_01/main'
include { PLOT_FIG_08 } from '../modules/static_plots/plot_fig_08/main'
include { PLOT_FIG_S01 } from '../modules/static_plots/plot_fig_s01/main'

workflow STATIC_PLOTS {
    main:
        ch_plot_fig_01_png = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_01/Fig_1.png')
        ch_plot_fig_01_drawio = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_01/Fig_1.drawio')
        
        PLOT_FIG_01(ch_plot_fig_01_png, ch_plot_fig_01_drawio)


        ch_plot_fig_08_png = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_09/Fig_9.png')
        ch_plot_fig_08_drawio = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_09/Fig_9.drawio')
        
        PLOT_FIG_08(ch_plot_fig_08_png, ch_plot_fig_08_drawio)


        ch_plot_fig_s01 = Channel.from('https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_10/Fig_10.png')
        
        PLOT_FIG_S01(ch_plot_fig_s01)
}