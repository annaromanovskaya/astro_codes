#!/usr/local/bin/python3
# find header, they contains only words
def is_header(line):
    return len([w for w in line.split() if w.isalpha()]) > 0

def find_adelman(path, wavelength, flux, star):

    file = open(path, 'r')
    Lines = file.read().split('\n')
    full_data = []
    vector = ''
    header = ''
    vec = []

    for line in Lines:
        if is_header(line):
            if len(vector) > 0:
                full_data.append([header, vector])
                header = ''
                vector = ''
            header += line
        else:
            vector += line 
    
    if len(vector) > 0:
        full_data.append([header, vector])
    
    # star = 'Gam Gem'
    for header, vector in full_data:
        if star in header:
            vec = vector.split()

    # wavelength = []
    # flux = []
    for i in range(0, len(vec)):
        if i%2 == 0:
            wavelength.append(vec[i].rstrip('.'))
            # change type from str to float
            wavelength = list(map(float, wavelength))
        else:
            flux.append(vec[i])
            # change type from str to float
            flux = list(map(float, flux))
    return wavelength, flux, star


def adelman2erg(wavelength, flux, m5556, V):
    # V = 4.94 #magnitude of filter V in HD 118022
    fluxes = []
    for (x, y) in zip(wavelength,flux):
        flux_calc = 3.46e-9*((5556.0/x)**2)*10**(-0.4*(V-0.026+y-m5556))
        fluxes.append(flux_calc)
    return wavelength, fluxes

    # return wavelength, fluxes, m5556

def find_catalog(path, wavelength, flux, catalog):
    file = open(path, 'r')
    Lines = file.read().split('\n')
    full_data = []
    vector = ''
    header = ''
    vec = []
    for line in Lines:
        if "# " in line:
            if len(vector) > 0:
                full_data.append([header, vector])
                header = ''
                vector = ''
            header += line
        else:
            vector += line

    if len(vector) > 0:
        full_data.append([header, vector])

    # catalog = 'Alekseeva'
    for header, vector in full_data:
        if catalog in header:
            vec = vector.split()
        else:
            pass

    if catalog == 'IUE':
        vec = []
        for header, vector in full_data:
            if 'IUE' in header:
                vec1 = vector.split()
                vec.extend(vec1)
            
    wavelength = []
    flux = []

    # [start_index::spacing]
    wavelength = vec[0::3]
    wavelength = list(map(float, wavelength))
    if catalog == 'TD1':
        flux = vec[2::3]
    else:
        flux = vec[1::3]
    flux = list(map(float, flux))

    return wavelength, flux