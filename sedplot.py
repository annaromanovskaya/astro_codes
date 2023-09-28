#!/usr/bin/python3
# pyright: reportUnusedImport=false


from ast import arg
import numpy as np
# from astropy.io import ascii
# import sys
import argparse, textwrap
import pandas as pd
import matplotlib.gridspec as gridspec
import spectrophotometry2erg as s2e
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from astropy.io import ascii
import warnings
import os

def main():
    parser = argparse.ArgumentParser()
    # parser.add_argument('-i', '--input', nargs='+', type=str, required=True)
    parser.add_argument('-star', '--star', type=str, required=True, help= textwrap.dedent('''\
        Write the star number from the HD catalogue.'''))
    parser.add_argument('-s', '--showsed', action='store_true', required=False, help= textwrap.dedent('''\
        Mode for interactive figure: 
            zoom in, 
            zoom out,
            save as PNG, 
            etc. 
        Exit from the server: ctrl+C (macOS)'''))
    parser.add_argument('-o', '--output', type=str, required=False, help= textwrap.dedent('''\
        Mode for not interactive figure: \
            Write output file name for .eps format. '''))
    args = parser.parse_args()

    # HD of star
    hd = args.star
    print('Plottind SED for HD{}...'.format(hd))
    # PATHS
    path_Adelman = r'/home/astroram/INASAN/SED/shulyak/databases/observations/Adelman/data.dat'
    path_orig_obs = r'/home/astroram/INASAN/SED/sed-package/HD{}_sed/Adelman/fit-HD{}-obs-original_9071g36.obs'.format(*(hd,hd))
    path_model_1 = r'/home/astroram/INASAN/SED/sed-package/HD{}_sed/Adelman/fit-HD{}-model-fit_9071g36.obs'.format(*(hd,hd))
    path_model_1conv = r'/home/astroram/INASAN/SED/sed-package/HD{}_sed/Adelman/fit-HD{}-model-convolved-fit_9071g36.obs'.format(*(hd,hd))

    # DATAFRAMES OF PHOTOMETRY DATA
    wavelength = []
    flux = []

    catalog_names = ['TD1', 'IUE', 'Adelman', 'Photometry', 'Breger', 'Alekseeva', 'Burnashev']
    for name in catalog_names:
        data = s2e.find_catalog(path_orig_obs, wavelength, flux, name)
        d = {'Wavelength': data[0], 'Flux': data[1]}
        globals()[name] = pd.DataFrame(data=d, index=data[0])
        # print(globals()[name].head())
        if not globals()[name].empty:
            print(name, 'photometry found!')
            continue
        else:
            print(name, 'photometry not found!')
            
    x = input('Continue? [Y/n]:')
    if x.lower() == 'y':
        # find flux value from 5000A in Adelman df
        # it need for correct plot models
        if not Adelman.empty:
            flux_5000 = Adelman.loc[float(5000)][1]
            last_digits = str(flux_5000)[-2:]
        # find flux value from 12390A in 2MASS df
        if not Photometry.empty:
            flux_12390 = Photometry.loc[float(12350)][1]#12390
            last_digits = str(flux_12390)[-2:]
        # find name of model:
        model = open(path_model_1, 'r')
        model_name = model.readline().strip('\n').split(' ')[1]    
            
        """
        # find of m5556 flux FOR ADELMAN
        wl_Adelman = data_Adelman[0]
        flux_Adelman_init = data_Adelman[1]
        for i in range(0, len(wl_Adelman)):
            if wl_Adelman[i] == float(5556):
                m5556 = flux_Adelman_init[i]
        # calculate flux in correct units
        data_Adelman2erg = s2e.adelman2erg(wl_Adelman, flux_Adelman_init, m5556, 1.92)
        flux_Adelman = data_Adelman2erg[1]
        """



        print('Reading model file...')
        # dictionary for change flux units
        units_dict = {'11': 10**(-8), '12': 10**(-9), '13': 10**(-10), '14': 10**(-11)}
        # units_dict = {'11': 10**(-9), '12': 10**(-10), '13': 10**(-11), '14': 10**(-12)}
        # models
        data_m1 = ascii.read(path_model_1, delimiter=' ')
        d_m1 = {'Wavelength': data_m1['col1'], 'Flux_init': data_m1['col2']}
        df_m1 = pd.DataFrame(data=d_m1)
        for key in units_dict:
            if last_digits == key:#check for last digits in observation for correctly plot models
                df_m1['Flux'] = df_m1['Flux_init'].apply(lambda x: x*units_dict[key])
        print('Model file found!')
        
        # check if file not empty
        if os.path.getsize(path_model_1) > 14:
            data_model_1conv = ascii.read(path_model_1conv, delimiter=' ')
            d_conv = {'Wavelength': data_model_1conv['col1'], 'Flux_init': data_model_1conv['col2']}
            df_conv = pd.DataFrame(data=d_conv)
            for key in units_dict:
                if last_digits == key:#check for last digits in observation for correctly plot models
                    df_conv['Flux'] = df_conv['Flux_init'].apply(lambda x: x*units_dict[key])
            print('Convolved model file added!')
        else:
            pass

        # split photometry dataframe into GAIA, 2MASS, etc
        filter_name = ['gaia','johnson','hipparcos','panstarrs', 'twomass', 'wise']
        for i in filter_name:
            path_filters = r'/home/astroram/IDL/photometry/{}.filters'.format(i)
            data_filter = open(path_filters, 'r')
            Lines = data_filter.read().split('\n')
            globals()[i] = []
            for line in Lines:
                if 'w0 = ' not in line:
                    continue
                else:
                    wavelength = float(line.replace('w0 = ',''))
                    globals()[i].append(wavelength)
            globals()[i+'_df'] = Photometry[Photometry['Wavelength'].isin(globals()[i])]
            
        
        dfs = {
        "Model {}".format(model_name) : {'df': df_m1        , 'color': 'black'  , 'kind': 'line'   , 'marker': False    , 'size': 40, 'fill': 'black', 'linestyle': False   , 'linewidth': 1.0}, 
        "Convolved model"  : {'df': df_conv      , 'color': 'blue'   , 'kind': 'scatter', 'marker': 'd', 'size': 50, 'fill': None, 'linestyle': False   , 'linewidth': 1.0}, 
        "IUE"            : {'df': IUE          , 'color': 'm' , 'kind': 'line'   , 'marker': False    , 'size': 30, 'fill': None, 'linestyle': 'dashed', 'linewidth': 1.0},
        "Adelman"        : {'df': Adelman      , 'color': 'black'  , 'kind': 'scatter', 'marker': 'o'      , 'size': 30, 'fill': None, 'linestyle': False   , 'linewidth': False}, 
        # "Alekseeva"      : {'df': Alekseeva    , 'color': 'green'  , 'kind': 'scatter', 'marker': 'x'      , 'size': 30, 'fill': 'none', 'linestyle': False   , 'linewidth': False}, 
        # "Breger"         : {'df': Breger       , 'color': 'grey'   , 'kind': 'scatter', 'marker': 'o'      , 'size': 30, 'fill': 'none', 'linestyle': False   , 'linewidth': False}, 
        "TD1"            : {'df': TD1          , 'color': 'm' , 'kind': 'scatter', 'marker': 'o'      , 'size': 30, 'fill': 'magenta', 'linestyle': False   , 'linewidth': False}, 
        # "Burnashev"        : {'df': Burnashev    , 'color': 'green'  , 'kind': 'scatter', 'marker': 'o'      , 'size': 30, 'fill': 'none', 'linestyle': False   , 'linewidth': False}, 
        # "GAIA DR3"         : {'df': gaia_df      , 'color': 'green'  , 'kind': 'scatter', 'marker': 's'      , 'size': 30, 'fill': 'green', 'linestyle': False   , 'linewidth': False},
        # "Johnson"          : {'df': johnson_df   , 'color': 'blue'   , 'kind': 'scatter', 'marker': 'o'      , 'size': 30, 'fill': 'blue', 'linestyle': False   , 'linewidth': False}, 
        # "Hipparcos"        : {'df': hipparcos_df , 'color': 'orange' , 'kind': 'scatter', 'marker': 's'      , 'size': 30, 'fill': 'orange', 'linestyle': False   , 'linewidth': False}, 
        # "PAN-STARRS"       : {'df': panstarrs_df , 'color': 'orange' , 'kind': 'scatter', 'marker': 'o'      , 'size': 30, 'fill': 'orange', 'linestyle': False   , 'linewidth': False}, 
        "2MASS"            : {'df': twomass_df   , 'color': 'red'    , 'kind': 'scatter', 'marker': 'o'      , 'size': 30, 'fill': None, 'linestyle': False   , 'linewidth': False},
        # "WISE "            : {'df': wise_df      , 'color': 'black'    , 'kind': 'scatter', 'marker': 'o'    , 'size': 30, 'fill': 'black', 'linestyle': False   , 'linewidth': False}
            }

        

            
        # if show for interactive mode
        if args.showsed == True:
            # create interactive plot
            import hvplot
            import hvplot.pandas
            import holoviews as hv

            # setting bokeh as backend
            hv.extension('bokeh')

            # going to use show() to open plot in browser
            from bokeh.plotting import show
            print('Drawing for interactive mode')
            # pd.set_option('plotting.backend', 'hvplot')
            # dict for the dataframes and their names

            plots = []
            for key, value in dfs.items():
                df = value['df']
                color = value['color']
                kind = value['kind']
                linewidth = value['linewidth']
                marker = value['marker']
                size = value['size']
                label=key
                fill = value['fill']
                # Создание графика для текущего датафрейма
                if kind == 'line':
                    plot = df.hvplot(x='Wavelength', y='Flux', kind=kind, color=color, width=1600, height=820, line_width=linewidth, legend=False, label=label).opts(show_legend=True)
                else:
                    plot = df.hvplot(x='Wavelength', y='Flux', kind=kind, line_color=color, width=1600, height=820, marker=marker, size=size, fill_color=fill, legend=False, label=label).opts(show_legend=True)
                plots.append(plot)

            # Объединение графиков в один
            combined_plot = hv.Overlay(plots)
                
            hvplot.show(combined_plot)
            print('Done!')
        

        # if show for drawing
        else:
            print('Drawing for saving mode')
            # Create 2x2 sub plots
            fig = fig = plt.figure(figsize=(15,8))
            gs = gridspec.GridSpec(2, 2, figure=fig)
            ax1 = fig.add_subplot(gs[0, 1])#UV
            ax2 = fig.add_subplot(gs[ : , 0])#Visual
            ax3 = fig.add_subplot(gs[1, 1])#IR
            # plt.figure()

            # lists for UV, Optical and IR ranges
            uv_list = ["Model {}".format(model_name), 'Convolved model', 'Adelman', 'IUE', 'TD1']
            opt_list = ["Model {}".format(model_name), 'Convolved model', 'Adelman', 'GAIA', 'PANSTARRS', 'Hipparcos', 'Alekseeva', 'Breger', 'Johnson', 'Burnashev']
            ir_list = ["Model {}".format(model_name), 'Convolved model', '2MASS', 'WISE']
            # upcont UV fluxes
            IUE['Flux_up'] = IUE['Flux'].apply(lambda x: x*1.05)
            TD1['Flux_up'] = TD1['Flux'].apply(lambda x: x*1.05)

            # cut IUE
            IUE = IUE[IUE['Wavelength'].between(1250, 3300)]
            
            # VISUAL
            # format Y-axis
            # def y_fmt(x, y):
            #     # return '${:.2e}'.format(x).replace('e', '') + '}$'
            #     return '{:.2}'.format(x).replace(x[len(x)-3:], '')
            # label = "Model {}".format(model_name)
            for key, value in dfs.items():
                df = value['df']
                color = value['color']
                kind = value['kind']
                linewidth = value['linewidth']
                marker = value['marker']
                size = value['size']
                label = key
                fill = value['fill']
                # ULTRAVIOLET
                if key in uv_list:
                    ax = ax1
                    ax.set_title('Ultraviolet', weight='bold')
                    ax1.set(xlim=(1000, 3600), ylim=(0.0*1e-08, 0.10*1e-08))
                    ax.set(xlim=(1000, 3600))
                    ax.set_xlabel(None)
                    ax.set_ylabel(r'F$_{\lambda}$, $erg/s/cm^2/\AA$', fontsize=12)
                    if kind == 'line':
                        df.plot(x="Wavelength", y="Flux", ax = ax, label=label, color=color, linewidth=linewidth, linestyle='-', legend=True)
                    else:
                        if fill == None:
                            df.plot(x='Wavelength', y='Flux', ax = ax, label=label, color=color, marker=marker, markersize=5, mfc='none', linewidth=0, legend=True)
                        else:
                            df.plot(x='Wavelength', y='Flux', ax = ax, label=label, color=color, marker=marker, markersize=5, linewidth=0, legend=True)
                    ax.legend(['t9071g36', 'Convolved model', 'Adelman', 'IUE up to 1.05', 'TD1 up to 1.05'], loc='lower right', fontsize = 12)
                    # ax.legend(loc='upper right', fontsize = 12)
                # VISUAL
                if key in opt_list:
                    ax = ax2
                    ax.set_title('Visual', weight='bold')
                    # ax2.set(xlim=(3000, 8000), ylim=(0, 10*1e-10))
                    ax.set(xlim=(3000, 8000))
                    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True, useOffset=False)
                    # ax.yaxis.set_major_formatter(tick.FuncFormatter(y_fmt))
                    ax.set_xlabel(r'Wavelength, $\AA$', fontsize=12)
                    ax.set_ylabel(r'F$_{\lambda}$, $erg/s/cm^2/\AA$', fontsize=12)
                    if kind == 'line':
                        df.plot(x="Wavelength", y="Flux", ax = ax, label=label, color=color, linewidth=linewidth, legend=True)
                    else:
                        if fill == None:
                            df.plot(x='Wavelength', y='Flux', ax = ax, label=label, color=color, marker=marker, markersize=5, linewidth=0, mfc='none', legend=True)
                        else:
                            df.plot(x='Wavelength', y='Flux', ax = ax, label=label, color=color, marker=marker, markersize=5, linewidth=0, legend=True)
                    ax.legend(loc='upper right', fontsize = 12)
                # INFRARED
                if key in ir_list:
                    ax = ax3
                    ax.set_title('Infrared', weight='bold')
                    ax3.set(xlim=(10000, 22000), ylim=(0, 1.5*1e-10))
                    ax.set(xlim=(10000, 22000))
                    ax.set_xlabel(r'Wavelength, $\AA$', fontsize=12)
                    if kind == 'line':
                        df.plot(x="Wavelength", y="Flux", ax = ax, label=label, color=color, linewidth=linewidth, legend=True)
                    else:
                        if fill == None:
                            df.plot(x='Wavelength', y='Flux', ax = ax, label=label, color=color, marker=marker, markersize=5, linewidth=0, mfc='none', legend=True)
                        else:
                            df.plot(x='Wavelength', y='Flux', ax = ax, label=label, color=color, marker=marker, markersize=5, linewidth=0, legend=True)
                    ax.legend(loc='upper right', fontsize = 12)
            # plt.suptitle('Omi Peg SED')
            plt.tight_layout()
            plt.show()
            ax.figure.savefig(args.output+'', format='eps', dpi=300)
            warnings.filterwarnings( "ignore")
            print('Done!')
    if x.lower() == 'n':
        print('End program')
        pass

if __name__ == '__main__':
    main()
    

# type: ignore