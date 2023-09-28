#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from matplotlib import rc
import pandas as pd
import math
import argparse
from ast import arg

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-p', '--print', type=str, required=False, help='Show abundance table. ions - for abundance in each species, aver - for averaged abundance.')
    args = parser.parse_args()
        
    # path_init = r'/Users/astroram/INASAN/STARS_A/HD145788/'
    # plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams.update({'font.size': 12})

    order = ['He',
    'C',
    'N',
    'O',
    'Na',
    'Mg',
    'Al',
    'Si',
    'S',
    'K',
    'Ca',
    'Sc',
    'Ti',
    'V',
    'Cr',
    'Mn',
    'Fe',
    'Co',
    'Ni',
    'Zn',
    'Sr',
    'Y',
    'Zr',
    'Ba',
    'La',
    'Ce',
    'Pr',
    'Nd',
    #  'Sm',
    #  'Eu',
    #  'Gd',
    #  'Tb',
    #  'Dy',
    #  'Er',
    #  'Yb',
    #  'Lu'
    ]

    order_dict = dict(zip(order, range(len(order))))

    order_nlte = ['He',
    'C',
    'N',
    'O',
    'Na',
    'Mg',
    'Al',
    'Si',
    'S',
    'K',
    'Ca',
    'Ti',
    'Fe',
    'Zn',
    'Sr',
    'Y',
    'Zr',
    'Ba']


    order_nlte_dict = dict(zip(order_nlte, range(len(order_nlte))))

    # SOLAR abundance
    path_sun = r'/home/astroram/INASAN/Reference_Abunds/Info_Comparison_Stars/Solar_abundance_2021.dat'
    data_sun = ascii.read(path_sun, delimiter=' ')
    df_dict_sun = {'El': data_sun['El'], 'init_abund': data_sun['Lodd_21']}
    df_sun_init = pd.DataFrame(data = df_dict_sun)
    df_sun = df_sun_init[df_sun_init['El'].apply(lambda x: x in order)]
    df_sun['Sun'] = df_sun['init_abund']-12.04

    # STAR abundance
    # path = r'/hdd/home/astroram/INASAN/STARS_A/HD145788/Line_abund_HD145788.dat'
    path = args.input
    data = ascii.read(path, delimiter=',')
    df_dict = {'Ion': data['col1'], 'Wavelength': data['col2'], 'Ei': data['col3'], 'loggf': data['col4'], 'loggf_ref': data['col10'], 'hfs_ref': data['col11'], 'iso_ref': data['col12']}
    df_init = pd.DataFrame(data = df_dict)

    # replacements for loggs hfs and iso references
    replacement = {
        'gf:': '',
        'hfs:': '',
        '\t\t\t': '\t\t',
        '\t\t': '\t'
        # "HELLO": "CORP"
    }

    df1 = df_init.replace(replacement, regex=True)
    df2 = df1['iso_ref'].str.split('\t', expand=True).rename(columns={0:'ISO_ref', 1:'LTE', 2:'NLTE'})
    df3 = df1['Ion'].str.strip("'").str.split(' ', expand=True).rename(columns={0:'Element', 1:'ion'})
    pd_List = [df3,df1,df2]
    df4 = pd.concat(pd_List, axis=1)
    df = df4.drop(['Ion','iso_ref'], axis=1)

    # create df for abundances
    df_a1 = df[['Element','LTE','NLTE']]
    df_a1['Element'] = df1['Ion'].str.strip("'")
    df_a2 = df_a1[['LTE','NLTE']].apply(pd.to_numeric, errors='coerce')
    df_a = pd.concat([df_a1['Element'],df_a2], axis=1)

    # create dataframe for IONS
    df_abund = df_a.groupby(by=['Element']).mean().round(2)
    df_abund_err = df_a.groupby(by=['Element']).std(ddof=0).round(2).rename(columns={'LTE':'LTE_err','NLTE':'NLTE_err'})
    df_ions1 = pd.concat([df_abund,df_abund_err], axis=1).reset_index()
    df_ions1[['El','Ion']] = df_ions1['Element'].str.split(' ',expand=True)
    df_ions = pd.merge(df_ions1, df_sun, on="El")
    df_ions['X/H lte'] = df_ions['LTE'] - df_ions['Sun'] # pure lte
    df_ions['X/H nlte'] = df_ions['NLTE'] - df_ions['Sun']
    df_ions['X/H lte+nlte'] = df_ions.apply(lambda row: row['X/H lte'] if np.isnan(row['X/H nlte']) else row['X/H nlte'], axis = 1) # for nlte+lte
    df_ions['ERR lte+nlte'] = df_ions.apply(lambda row: row['LTE_err'] if np.isnan(row['NLTE_err']) else row['NLTE_err'], axis = 1)
    elements = ['Co 2', 'La 2', 'N 1', 'Si 1']
    columns = ['ERR lte+nlte', 'LTE_err']
    df_ions.loc[df_ions['Element'].isin(elements), columns] = 0.20
    df_neutral = df_ions[df_ions.apply(lambda row: row['Ion'] == '1' and row['El'] in order or row['El'] in ['Ba'], axis = 1)]
    df_1st = df_ions[df_ions.apply(lambda row: row['Ion'] in ['2','3'] and row['El'] in order and row['El'] not in ['Ba'], axis = 1)]

    # # create datagrame for AVERAGE
    df_a['El'] = df3['Element']
    # df_abund = df_a.groupby(by=['El']).mean().round(2)
    df_abund = df_a.drop(columns=['Element']).groupby(by=['El']).mean().round(2)
    df_abund_err = df_a.drop(columns=['Element']).groupby(by=['El']).std(ddof=0).round(2).rename(columns={'LTE':'LTE_err','NLTE':'NLTE_err'})
    df_average1 = pd.concat([df_abund,df_abund_err], axis=1).reset_index()
    df_average = pd.merge(df_average1, df_sun, on="El")
    df_average['X/H lte'] = df_average['LTE'] - df_average['Sun'] # pure lte
    df_average['X/H nlte'] = df_average['NLTE'] - df_average['Sun']
    df_average['X/H lte+nlte'] = df_average.apply(lambda row: row['X/H lte'] if np.isnan(row['X/H nlte']) else row['X/H nlte'], axis = 1) # for nlte+lte
    df_average['ERR lte+nlte'] = df_average.apply(lambda row: row['LTE_err'] if np.isnan(row['NLTE_err']) else row['NLTE_err'], axis = 1)
    elements = ['Co', 'La', 'N']
    columns = ['ERR lte+nlte', 'LTE_err']
    df_average.loc[df_average['El'].isin(elements), columns] = 0.20
    
    if args.print == 'ions':
        df = df_ions.drop(['El','Ion', 'init_abund'], axis=1).round(2)
        print(df)
        df.to_csv('Abundances_ions.dat', sep='\t', index=False)
    if args.print == 'aver':
        df = df_average.drop(['init_abund'], axis=1).round(2)
        print(df)
        df.to_csv('Abundances_average.dat', sep='\t', index=False)

    # # replace the values of S to S1 (except S average)!!!
    # colname = ['GG-sun_lte', 'GG-sun_nlte', 'OP-sun_lte', 'OP-sun_nlte', 'TV-sun_lte', 'TV-sun_nlte']
    # df_final_mean[colname].loc[['S']] = df_final[colname].loc[df_final['Element'] == 'S 1']
    # df_final_mean[['NC-sun_lte', 'NC-sun_nlte']].loc[['S']] = df_final[['NC-sun_lte', 'NC-sun_nlte']].loc[df_final['Element'] == 'S 2']

    # DRAWING

    # IONS ALL
    # PURE LTE abundance
    X_ions = df_ions['El']
    Y_ions_lte = df_ions['X/H lte']
    dY_ions_lte = df_ions['LTE_err']
    # LTE and final NLTE abundance
    Y_ions_nlte = df_ions['X/H lte+nlte']
    dY_ions_nlte = df_ions['ERR lte+nlte']
    # IONS neutral+1st
    # PURE LTE abundance
    X_neutral = df_neutral['El']
    Y_neutral_lte = df_neutral['X/H lte']
    dY_neutral_lte = df_neutral['LTE_err']
    # LTE and final NLTE abundance
    Y_neutral_nlte = df_neutral['X/H lte+nlte']
    dY_neutral_nlte = df_neutral['ERR lte+nlte']
    # PURE LTE abundance
    X_1st = df_1st['El']
    Y_1st_lte = df_1st['X/H lte']
    dY_1st_lte = df_1st['LTE_err']
    # LTE and final NLTE abundance
    Y_1st_nlte = df_1st['X/H lte+nlte']
    dY_1st_nlte = df_1st['ERR lte+nlte']

    # AVERAGE
    # PURE LTE abundance IONS
    X_average = df_average['El']
    Y_aver_lte = df_average['X/H lte']
    dY_aver_lte = df_average['LTE_err']

    # LTE and final NLTE abundance
    Y_aver_nlte = df_average['X/H lte+nlte']
    dY_aver_nlte = df_average['ERR lte+nlte']

    fig = plt.figure(figsize=(7,4))
    plt.ylabel('[X/H]', fontsize=12)

    if args.print == 'ions':
        # IONS
        # plt.errorbar([order.index(i) for i in X_ions], Y_ions_nlte, dY_ions_nlte, color='r', fmt='o', capsize=4, linewidth=1)
        plt.errorbar([order.index(i) for i in X_neutral], Y_neutral_lte, dY_neutral_lte, color='k', fmt='o', capsize=4, linewidth=1)
        plt.errorbar([order.index(i) for i in X_1st], Y_1st_lte, dY_1st_lte, color='r', fmt='o', capsize=4, linewidth=1, markerfacecolor='none', markeredgewidth=0.7)
    if args.print == 'aver':
        # AVERAGE
        plt.errorbar([order.index(i) for i in X_average], Y_aver_nlte, dY_aver_nlte, color='k', fmt='o', capsize=4, linewidth=1)
    
    plt.xticks([i for i in range(len(order))], order)
    plt.yticks(np.arange(-3, 8, step=0.5))

    # plt.legend([r'HD 145788 abundance'], loc='upper left', fontsize=10)
    plt.legend([r'Neutral + REE 1st', r'1st ions'], loc='upper left', fontsize=10)

    plt.xlim(-1, 28)
    plt.ylim(-1.2, 2)
    plt.plot([-2, 26],[0,0], color='k', linestyle='-', linewidth=0.7)

    # if order_nlte == True:
        # plt.vlines(x = [order.index(i) for i in order_nlte], ymin = -1.2, ymax = 2.0, linewidth = 7.0,
        #         #    colors = 'gainsboro',
        #             colors = 'k',
        #         alpha = 0.1
        #         )
    

    plt.title(r'$\alpha$Equ A')
    fig.tight_layout()
    fig.savefig(args.output, format='pdf', dpi=300)
    plt.show()