#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 14:51:18 2023

@author: yarivv
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions import *

pipeline_rep = "/user/home/yarivv/ATMOSPHERIX_DATA_RED/Planets/"
#Choose planet and import parameters
planet_name = "WASP-33b"
planet_file = pipeline_rep + planet_name + "_params.py"

#emission or transit
emission = True
prt_version = 3
if prt_version == 2:
    from petitRADTRANS import Radtrans
    from petitRADTRANS import nat_cst as nc
elif prt_version ==3:
    from petitRADTRANS.radtrans import Radtrans
    from petitRADTRANS import physical_constants as cst


with open(planet_file) as file:
    exec(file.read())

#%%

#Setting up two point T-P profile with strong temperature inversion 
#Parameters for WASP-33 chosen according to Finnerty et al. (2023) paper retrieved parameters
#H- opacity not included as ~constant in K band but careful over other bands H- may be important
T_bottom = 3500.
T_top = 4500.
pressures = np.logspace(-6, 2, 100)
meso = np.ones_like(np.flatnonzero(pressures<=10**-3))*T_top
tropo = np.ones_like(np.flatnonzero(pressures>=10**-1))*T_bottom
strato = np.linspace(T_top, T_bottom, len(pressures)-len(meso)-len(tropo))
temperature = np.concatenate((meso,strato,tropo), axis=0)

#Gravity for WASP-33b
gravity = G*Mp/(Rp*1e3)**2 #37.62m/s2 from Wong et al (but inconsistent with their masses)
print(gravity)



#Set up mass fractions and mean molecular weight for pRT atmosphere
mass_fractions = {}
mass_fractions['H2'] = 0.7 * np.ones_like(temperature)
mass_fractions['He'] = 0.28 * np.ones_like(temperature)
mass_fractions['13C-16O'] = 1e-3 * np.ones_like(temperature)
#mass_fractions['OH'] = 1*10**-2 * np.ones_like(temperature)
#mass_fractions['Ti'] = 1*10**-4 * np.ones_like(temperature)
#mass_fractions['Fe'] = 1*10**-4 * np.ones_like(temperature)

MMW = 2.33 * np.ones_like(temperature)

#Setting up atmosphere object with CO, He and H 
if prt_version == 2:
    atmosphere = Radtrans(line_species = ['13C-16O'], #, 'CO_all_iso', "Fe", "Ti"...
                      rayleigh_species = ['H2', 'He'],
                      continuum_opacities = ['H2-H2', 'H2-He'],
                      wlen_bords_micron = [1.5, 2.5],
                      do_scat_emis=True,
                      mode='lbl')

    atmosphere.setup_opa_structure(pressures)
    #Compute planet spectrum
    if emission:
        atmosphere.calc_flux(temperature, mass_fractions, gravity, MMW)
    else:
        atmosphere.calc_flux(temperature, mass_fractions, gravity, MMW)


    #Convert Units
    wavelengths = nc.c / atmosphere.freq / 1e-8                             #[Angstrom]
    flux_planet = atmosphere.flux * nc.c / (wavelengths*1e-8)**2            #[erg/s/cm2/cm]
elif prt_version == 3:
    atmosphere = Radtrans(
                      pressures = np.logspace(-6, 2, 100),
                      line_species = ['13C-16O'], #, 'CO_all_iso', "Fe" ...
                      rayleigh_species = ['H2', 'He'],
                      gas_continuum_contributors = ['H2-H2', 'H2-He'],
                      wavelength_boundaries = [0.3, 2.5],
                      scattering_in_emission = True,
                      line_opacity_mode='lbl')
    #Compute emission spectrum
    if emission:
        frequencies, flux, _ = atmosphere.calculate_flux(
        temperatures=temperature,
        mass_fractions=mass_fractions,
        mean_molar_masses = MMW,
        reference_gravity = gravity,
        frequencies_to_wavelengths=False
        )
        wavelengths = cst.c / frequencies / 1e-8                             #[Angstrom]
        flux_planet = flux * cst.c / (wavelengths*1e-8)**2                   #[erg/s/cm2/cm]
    
    #Convert units
    #    wav_micron = wavelengths / 1e4
    #    fs = bb_to_flux(bb(wavelengths*10**-10, Ts))                    #[erg/s/cm2/cm]
    #    Rp_Rs =  - flux_planet/fs * (Rp/Rs)**2                      #Transit Depth (negative for emission)
    
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



#%%
volume_fractions = {}
volume_fractions['H2'] = mass_fractions['H2'] * MMW / 2
volume_fractions['He'] = mass_fractions['He'] * MMW / 4
volume_fractions['13C-16O'] = mass_fractions['13C-16O'] * MMW / 29
#volume_fractions['OH'] = mass_fractions['OH'] * MMW / 17  #
#volume_fractions['Fe'] = mass_fractions['Fe'] * MMW / 56  # 
#volume_fractions['Ti'] = mass_fractions['Ti'] * MMW / 48  # 

#%%
#Plot Emission spectrum
fig, ax = plt.subplots(figsize=(12,5))
ax.plot(wavelengths/1e4, flux_planet, linewidth='0.25', c="black")
#ax.set_xlim([2.1,2.5])
#ax.set_ylim([15,30])
ax.set_xlabel('Wavelength (microns)')
ax.set_ylabel(r'Planet flux $F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ cm$^{-1}$)')
ax.set_title("Model Planet Spectrum")
plt.savefig("/user/home/yarivv/ATMOSPHERIX_DATA_RED/Figures/Model_13CO.png")
#plt.show()

#save spectrum



save_dir  = '/user/home/yarivv/ATMOSPHERIX_DATA_RED/Models/Results/'
save_name = 'WASP33_13CO'

np.savetxt(save_dir+'lambdas'+save_name+'.txt', wavelengths/10)  #[nm]
np.savetxt(save_dir+'Fp'+save_name+'.txt', flux_planet)          #[erg/s/cm2/cm]




