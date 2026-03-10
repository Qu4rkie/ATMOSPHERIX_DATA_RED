#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 11:56:56 2024

@author: yarivv
"""

#Constants
Rsun     = 6.9634e5         #Solar radius (km)
Rjup     = 6.9911e4       #Jupiter Radius [km]
Au       = 1.4960e8        #AU [km]
pc       = 3.08567e13      #pc [km]
Msun     = 1.99e30         #[kg]
Mjup     = 1.898e27        #[kg]


###############################################################################
### Parameters for TOI2109b
### Ephemerides (to compute orbital phase)
T0         =  2453885.6888            #Mid-Transit epoch Alvaro-Montes et al, 2025
#T0        =  2459378.459370          #Mid-Transit epoch Wong et al, 2021
Porb       =  0.67247455              #Orbital period [d] Alvaro-Montes et al, 2025
#Porb      =  0.67247414              #Orbital period [d] Wong et al, 2021
T_peri     =  2459378.459370          #Time of peri astron passage for an elliptical orbit
T_eclipse  =  T0 + Porb/2           #Mid-Eclipse time

T_star       = 6540                  #Star temperature [K]
Fp       = 0.0015                    #Planet/Star emission ratio (used if in emission) 

transiting  = True  ## Only used in emission, to calculate the window function during eclipse if need be

### Transit parameters -- Compute the transit window
### Using batman python package https://lweb.cfa.harvard.edu/~lkreidberg/batman/
### Get the limb-darkening coefficients in H band from Claret+2011: https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/529/A75
Rp       = 1.347 * Rjup        #Planet radius  [km] WOng et al. 2021
Rs       = 1.698 * Rsun        #Stellar radius [km]      ""
Ms       = 1.453 *Msun             #Stellar mass [kg]     ""
Mp       = 5.02  *Mjup             #Planet mass [kg]  ""
ip       = 70.74               #Transit incl.  [deg]  Wong et al, 2021
ap       = 2.268            #Semi-maj axis  [R_star]    ""
ep       = 0.                 #Eccentricity of Pl. orbit
wp       = 90.                #Arg of periaps [deg]
ld_mod   = "quadratic"         #Limb-darkening model ["nonlinear", "quadratic", "linear"]
ld_coef  = [0.28, -0.40]     #Limb-darkening coefficients (Harre et al. 2024)


### Stellar radial velocity info
Ks        =  0.860   #RV semi-amplitude of the star orbital motion due to planet [km/s] Wong et al, 2021
V0        =  -25.64  #Stellar systemic velocity [km/s]  ""
Vsini     = 81.2     #Stellar Rotation velocity [km/s] ""
sy_dist   = 262.041 * pc    #System distance [km]

#Derived Parameters
Vrot      = 2*np.pi*Rp / (Porb*24*3600)     #Synchronous rotation velocity [km/s]
###############################################################################