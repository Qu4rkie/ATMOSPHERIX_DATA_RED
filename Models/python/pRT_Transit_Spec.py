#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:51:18 2023

@author: yarivv
"""

import numpy as np
import matplotlib.pyplot as plt
from petitRADTRANS.radtrans import Radtrans
from petitRADTRANS import physical_constants as cst
from petitRADTRANS.physics import temperature_profile_function_guillot_global
from Functions import *

pipeline_rep = "/Users/yarivv/ATMOSPHERIX_DATA_RED/"

#Choose planet and import parameters
planet_name = "WASP-43b"
planet_file = pipeline_rep + planet_name + "_params.py"

#emission or transit
emission = True

with open(planet_file) as file:
    exec(file.read())

#%%
#Setting up atmosphere object with CO, He and H 

atmosphere = Radtrans(
                      pressures = np.logspace(-6, 2, 100),
                      line_species = ['16O-1H'], #, 'CO_all_iso'
                      rayleigh_species = ['H2', 'He'],
                      gas_continuum_contributors = ['H2-H2', 'H2-He'],
                      wavelength_boundaries = [0.3, 2.5],
                      scattering_in_emission = True,
                      line_opacity_mode='lbl')

#Gravity
reference_gravity = G*(Mp*1e3)/(Rp*1e3)**2    #[m/s2]
reference_pressure = 0.1                     #[bar]

pressures = atmosphere.pressures*1e-6 # cgs to bar
infrared_mean_opacity = 0.01
gamma = 0.4
intrinsic_temperature = 100
equilibrium_temperature = 1000 #Nightside


#Isothermal profile
#temperatures = 1400 * np.ones_like(atmosphere.pressures)
#Guillot Profile
temperatures = temperature_profile_function_guillot_global(
    pressures=pressures,
    infrared_mean_opacity=infrared_mean_opacity,
    gamma=gamma,
    gravities=reference_gravity,
    intrinsic_temperature=intrinsic_temperature,
    equilibrium_temperature=equilibrium_temperature
)

#Set up mass fractions and mean molecular weight for pRT atmosphere
mass_fractions = {}
mass_fractions['H2'] = 0.7 * np.ones_like(temperatures)
mass_fractions['He'] = 0.2 * np.ones_like(temperatures)
mass_fractions['1H2-16O'] = 1e-3 * np.ones_like(temperatures)


MMW = 2.33 * np.ones_like(temperatures)

#Compute emission spectrum
if emission:
    frequencies, flux, _ = atmosphere.calculate_flux(
    temperatures=temperatures,
    mass_fractions=mass_fractions,
    mean_molar_masses = MMW,
    reference_gravity = reference_gravity,
    frequencies_to_wavelengths=False
    )
    wavelengths = cst.c / frequencies / 1e-8                             #[Angstrom]
    flux_planet = flux * cst.c / (wavelengths*1e-8)**2                   #[erg/s/cm2/cm]
    
    #Convert units
    wav_micron = wavelengths / 1e4
    fs = bb_to_flux(bb(wavelengths*10**-10, Ts))                    #[erg/s/cm2/cm]
    Rp_Rs =  - flux_planet/fs * (Rp/Rs)**2                      #Transit Depth (negative for emission)
    
else:
    wavelengths, transit_radii, _ = atmosphere.calculate_transit_radii(
        temperatures=temperatures,
        mass_fractions=mass_fractions,
        mean_molar_masses=MMW,
        reference_gravity=reference_gravity,
        planet_radius=Rp*1e5,                   #[cm]
        reference_pressure=reference_pressure
        )
    
    #Convert units
    wav_micron = wavelengths * 1e4
    Rp_Rs      = transit_radii * 1e-5 / Rs

#Plot spectrum
fig, ax = plt.subplots(figsize = (10, 6))
ax.plot(wav_micron, Rp_Rs)
ax.set_xlabel('Wavelength [microns]')
ax.set_ylabel(r'Transit depth [$R_p$/$R_*$]')
plt.show()

#%%
#Plot T-P profile and Opacities
fig = plt.figure(figsize = (8,8))
ax1 = plt.subplot(121)
ax1.plot(temperatures, pressures)
ax1.set_yscale('log')
ax1.set_ylim([1e0, 1e-5])
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Pressure (bar)')

ax2 = plt.subplot(122)
for i in mass_fractions:
    ax2.plot(mass_fractions[i], pressures, label=i)
ax2.legend()
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylim([1e0, 1e-5])
ax2.set_xlabel('Mass Fraction')
plt.show()
#%%
#Plot Emission spectrum
# fig, ax = plt.subplots(figsize=(12,5))
# ax.plot(nc.c/atmosphere.freq/1e-4, atmosphere.flux/1e-6, linewidth='0.25', c="black")
# ax.set_xlim([2.1,2.5])
# ax.set_ylim([15,30])
# ax.set_xlabel('Wavelength (microns)')
# ax.set_ylabel(r'Planet flux $F_\nu$ (10$^{-6}$ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)')
# ax.set_title("Model Planet Spectrum")
# plt.show()

#%%
#save spectrum

save_dir  = '/Users/yarivv/ATMOSPHERIX_DATA_RED/Data_Simulator/Model/Results/'
save_name = 'WASP_43b_H2O_Guillot_nightside'


np.savetxt(save_dir+save_name+'.dat', np.column_stack((wav_micron, Rp_Rs)))


