#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from astropy.io import ascii
import pandas as pd
import argparse, textwrap

def main():
    parser = argparse.ArgumentParser(description='some information',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', type=str, required=True, help= textwrap.dedent('''\
        Write the path of input file '''))
    parser.add_argument('-s', '--spectrograph', type=str, required=True, help= textwrap.dedent('''\
        Choose the spectrograph:
            harps - HARPS R = 80000
            sophie - SOPHIE R = 75000
            elodie - ELODIE R = 48000
            feros - FEROS R =
            uves - UVES R = 
            espadons - ESPADONS R = '''))
    parser.add_argument('-o', '--output', type=str, required=False, help= textwrap.dedent('''\
        Write the path of output file. '''))
    parser.add_argument('-p', '--plot', action='store_true', required=False, help= textwrap.dedent('''\
        Plot figure normalised flux vs wavelength. '''))
    parser.add_argument('-wl1', '--wavelength1', type=float, required=False, help= textwrap.dedent('''\
        Write first value of wavelength window for extract to file. '''))
    parser.add_argument('-wl2', '--wavelength2', type=float, required=False, help= textwrap.dedent('''\
        Write second value of wavelength window for extract to file. '''))

    args = parser.parse_args()

    if args.spectrograph == 'harps':
        file = fits.open(args.input)
        header = file[1].columns
        print(header)
        data = file[1].data
        wavelength = data[0][0]
        flux = data[0][1]
        err = data[0][2]
    
    if args.spectrograph == 'uves':
        file = fits.open(args.input)
        header = file[1].columns
        print(header)
        data = file[1].data
        wavelength = data[0][0]
        flux = data[0][3]
        err = data[0][4]

    if args.spectrograph == 'narval':
        file = open(args.input, 'r')
        Lines = file.readlines()
        wl = []
        flux = []
        for line in Lines[2:]:
            line = line.rstrip('\n').replace('  ',' ').split(' ')
            wl.append(float(line[1]))
            flux.append(float(line[2]))
        wavelength = [i*10 for i in wl]
        data_df = {'Wavelength': [round(elem, 4) for elem in wavelength], 'Flux': [round(elem, 4) for elem in flux]}
        df = pd.DataFrame(data=data_df).astype(np.float32)

        order = 0
        path = f"spc_{order}.dat"
        for i in range(0, len(df)):
            delta = df['Wavelength'].iloc[i] - df['Wavelength'].iloc[i-1]
            if delta < 0:
                path = f"spc_{order}.dat"
                order += 1
            else:
                with open(path, 'a') as file:
                    file.write(f"{df['Wavelength'].iloc[i]} {df['Flux'].iloc[i]}\n")

    if args.spectrograph == 'espadons':
        file = fits.open(args.input)
        header = file[0].header
        data = file[0].data

        npix = float(header['NAXIS1'])
        naxis = float(header['NAXIS'])
        ncols = float(header['NAXIS2'])

        header = file[0].header
        data = file[0].data
        wavelength = [float(i)*10 for i in data[0]]
        flux = data[1]
        err = data[2]

        data_df = {'Wavelength': [round(elem, 4) for elem in wavelength], 'Flux': [round(elem, 4) for elem in flux]}
        df = pd.DataFrame(data=data_df).astype(np.float32)

        order = 0
        for i in range(1, len(df)-1):
            delta = df['Wavelength'].iloc[i] - df['Wavelength'].iloc[i-1]
            if delta < 0:
                path = f"spc_{order}.dat"
                order += 1
            else:
                with open(path, 'a') as file:
                    file.write(f"{df['Wavelength'].iloc[i]} {df['Flux'].iloc[i]}\n")

        
    if args.spectrograph == 'sophie' or args.spectrograph == 'elodie':
        file = fits.open(args.input)
        header = file[0].header
        data = file[0].data
        crval = float(header['CRVAL1'])
        crpix = int(header['CRPIX1'])
        cdelt = float(header['CDELT1'])
        npix = float(header['NAXIS1'])
        wavelength = (np.arange(1.0,npix+1) - crpix)*cdelt + crval
        flux = data.astype(np. float64)

    if args.spectrograph == 'elodie':
        file = fits.open(args.input)
        header = file[0].header
        data = file[0].data
        crval = float(header['CRVAL1'])
        crpix = int(header['CRPIX1'])
        cdelt = float(header['CDELT1'])
        npix = float(header['NAXIS1'])
        wavelength = (np.arange(1.0,npix+1) - crpix)*cdelt + crval
        flx = data.astype(np.float64)
        flux = [x*10**(-8) for x in flx]
    
    if args.spectrograph == 'lamost':
        file = fits.open(args.input)
        header = file[0].header
        data = file[0].data
        wavelength = data[2]
        flux = data[1]

    # if args.spectrograph == 'feros':
    # if args.spectrograph == 'uves':

    data_df = {'Wavelength': [round(elem, 4) for elem in wavelength], 'Flux': ["{:e}".format(elem) for elem in flux]}
    df = pd.DataFrame(data=data_df)
    # print(df.head())
    if args.plot == True:
        plt.figure(figsize=(12, 6), dpi=80)
        # df.plot(x='Wavelength', y='Flux')
        plt.plot(wavelength,flux)
        plt.xlabel('Wavelength (Angstrom)')
        plt.ylabel('Normalized flux')
        plt.show()

    # filter for Hydrogen lines
    # df.apply(lambda x: x >= args.wavelength1 and x <= args.wavelength2)
    if args.wavelength1:
        df_wl = df[df['Wavelength'].between(args.wavelength1, args.wavelength2)]
        # df_wl = df.apply(lambda x: x >= args.wavelength1 and x <= args.wavelength2)
        df_wl.to_csv(args.output, sep=' ', header=False, index=False)
    else:
        df.to_csv(args.output, sep=' ', header=False, index=False)

if __name__ == '__main__':
    main()