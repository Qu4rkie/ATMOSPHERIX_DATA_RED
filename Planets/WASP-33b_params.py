#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:14:37 2024

@author: yarivv
"""

#Useful Constants
Rsun     = 6.9634e5        #Solar radius (km)
Rjup     = 6.9911e4        #Jupiter Radius [km]
Au       = 1.4960e8        #AU [km]
Msun     = 1.99e30         #[kg]
Mjup       = 1.898e27        #[kg]


###############################################################################
### Parameters for WASP-33
### Ephemerides (to compute orbital phase)
T0       =  2454163.22367                #Mid-Transit epoch, Zhang et al.
Porb     =  1.21987089                      #Orbital period [d] Zhang et al.

Ts       = 7430                      #Star temperature [K]
Fp       = 0.001678                       #Planet/Star emission ratio (bb emission at 2.3 micron)
Te       = T0 + Porb/2                    #Mid-Eclipse time 

### Transit parameters -- Compute the transit window
### Using batman python package https://lweb.cfa.harvard.edu/~lkreidberg/batman/
### Get the limb-darkening coefficients in H band from Claret+2011: https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/529/A75
Rp       = 1.603 * Rjup        #Planet radius  [km]
Rs       = 1.509 * Rsun        #Stellar radius [km] 
Ms       = 1.561 *Msun             #Stellar mass [kg]
Mp       = 2.16  *Mjup             #Planet mass [kg]
ip       = 86.2               #Transit incl.  [deg]
ap       = 3.69                #Semi-maj axis  [R_star]
ep       = 0.0                 #Eccentricity of Pl. orbit
wp       = 0.0                 #Arg of periaps [deg]
ld_mod   = "quadratic"         #Limb-darkening model ["nonlinear", "quadratic", "linear"]
ld_coef  = [0.0951,0.1586]     #Limb-darkening coefficients 


### Stellar radial velocity info
Ks        = 0.304    #RV semi-amplitude of the star orbital motion due to planet [km/s] 
V0        = -2.7    #Stellar systemic velocity [km/s] (SIMBAD) ( -3.02 km/s in Nugroho et al 2017)
###############################################################################