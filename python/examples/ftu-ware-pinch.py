"""
Neoclassical fluxes for modelling FTU pulse 30582 r/a=0.6 t=0.3

densities are in [1e19 m^-3], temperature values in [keV]
"""
from math import pi
import numpy as np

import neoart


rln = 2.9061
eps = 0.185
q = 2.761
bn = 5.685  # [T]
rgeo = 0.97  # [m]
vloop = -1.3  # [V]

electron = {'density': 6.029,
            'temperature': 0.57,
            'density gradient': rln,
            'temperature gradient': 9.92}
deuterium = {'density': 3.23546285,
             'temperature': 0.47,
             'density gradient': rln,
             'temperature gradient': 6.9258}
lithium = {'density': 0.9313394214,  # from quasi-neutrality
           'temperature': 0.47,
           'density gradient': rln,
           'temperature gradient': 6.9258}

species = [('electron', electron),
           ('deuterium', deuterium),
           ('lithium', lithium)]

for label, profiles in species:  # add eps for each species
    profiles['eps'] = eps

geometry = neoart.CircularGeometry(minor_radius=rgeo * eps, major_radius=rgeo,
                                   safety_factor=q, magnetic_field=bn)
neo = neoart.Neoart(species, geometry)

# compute the fluxes with/without the parallel electric field
flux_eparr = neo.fluxes(eparr=vloop / (2 * pi))
flux = neo.fluxes()

density = np.array(neo._species.get_profile('density'))
ware_flux = flux_eparr - flux
ware_pinch = ware_flux / density / 1e19

# get rid of the singular dimensions
ware_flux = ware_flux.squeeze()
flux_eparr = flux_eparr.squeeze()
ware_pinch = ware_pinch.squeeze()

# print results
order = '[' + ', '.join(s[0] for s in species) + ']'
print 'Total neoclassical particle flux in [m^-2 s^-1] for', order
print flux_eparr[:, 0], '[m^-2 s^-1]'
print
print 'Ware contribution for', order
print ware_flux[:, 0], '[m^-2 s^-1]'
print
print 'Ware-pinch for ', order
print ware_pinch[:, 0], '[m/s]'
