#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 03/2022 12:54:11 2021

@author: Baptiste & Florian
"""
import numpy as np
import os
from astropy.io import fits
from astropy.time import Time
from barycorrpy import get_BC_vel, JDUTC_to_BJDTDB
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
from global_parameters import c0, G

import batman
import glob


# -----------------------------------------------------------
# This function reads a time series of DRS-provided SPIRou files
# It stores some of the relevant information into "Order" objects
# and returns time series relevant for the analysis
# The spectra to read must be stored in a dedicated repository
# For the time being, the function can only read t.fits extensions 
# -----------------------------------------------------------
def read_data_spirou(repp,list_ord,nord):

    """
    --> Inputs:     - repp:      Path to the directory containing all the '.fits' files to read
                                 NOTE: files must be ordered in the chronologic order
                    - list_ord:  List of Order object
                    - nord:      Number of orders -- 49 for SPIRou

    --> Outputs:    - Attributes of Order objects:
                      1. W_raw (Wavelengths vectors)
                      2. I_raw (Time series of spectra)
                      3. blaze (Time series of blaze functions)
                      4. A_raw (Time series of telluric spectra computed from the DRS)
                      5. SNR (list of order S/N values for all observations)
                    - list_ord: upgraded list of orders
                    - airmass: airmass value for each observation
                    - bjd: time vector
                    - snr_mat: 2D matrix containing the S/N value for each observation and order (N_observation,N_order)
    """
    print(repp)
    directory = os.listdir(repp)
    nam_t = []
    
    #Get only t.fits files (in case other files in the directory)
    for filename in directory:
        if "t.fits" in filename:
            nam_t.append(filename)
    
    nam_t     = sorted(nam_t)
    nobs      = len(nam_t)
    airmass   = np.zeros(nobs)
    bjd       = np.zeros(nobs)
    berv      = np.zeros(nobs)
    snr_mat   = np.zeros((nobs,nord))
    for nn in range(nobs):
        nmn          = repp + "/" + str(nam_t[nn])
        hdul_t       = fits.open(nmn)
        airmass[nn]  = float(hdul_t[0].header["AIRMASS"])   #airmass
        bjd[nn]      = float(hdul_t[1].header["BJD"])       #mid-exposure time [BJD-TDB]
        berv[nn]     = float(hdul_t[1].header["BERV"])      #Barycentric Earth Velocity [km/s]
        i            = np.array(hdul_t[1].data,dtype=float) # intensity spectrum
        w            = np.array(hdul_t[2].data,dtype=float) # wavelength vector [nm] - confirm unit
        bla          = np.array(hdul_t[3].data,dtype=float) # blaze vector
        atm          = np.array(hdul_t[4].data,dtype=float) # telluric spectrum
        ### Get S/N values
        for mm in range(nord):
            num = 79 - list_ord[mm].number
            if num < 10: key = "EXTSN00" + str(num)
            else: key = "EXTSN0" + str(num)
            sn  = float(hdul_t[1].header[key]) # S/N for each order
            list_ord[mm].SNR.append(sn)
            snr_mat[nn,mm] = sn
        hdul_t.close()
        ## Store Order's attributes
        for mm in range(nord):
            O = list_ord[mm]
            num = 79 - list_ord[mm].number
            O.W_raw.append(w[num])
            O.I_raw.append(i[num])
            O.blaze.append(bla[num])
            O.I_atm.append(atm[num])
    for mm in range(nord):
        O       = list_ord[mm]
        O.SNR   = np.array(O.SNR,dtype=float)
        O.W_raw = np.array(O.W_raw,dtype=float)
        O.I_raw = np.array(O.I_raw,dtype=float)
        O.blaze = np.array(O.blaze,dtype=float)
        O.I_atm = np.array(O.I_atm,dtype=float)
    return list_ord,airmass,bjd,berv,snr_mat
 
 
# -----------------------------------------------------------
# This function reads a time series of APERO-DRS NIRPS files
# It stores some of the relevant information into "Order" objects
# and returns time series relevant for the analysis
# The spectra to read must be stored in a dedicated repository
# For the time being, the function can only read t.fits extensions 
# -----------------------------------------------------------
def read_data_nirps(repp,list_ord,nord):

    """
    --> Inputs:     - repp:      Path to the directory containing all the '.fits' files to read
                                 NOTE: files must be ordered in the chronologic order
                    - list_ord:  List of Order object
                    - nord:      Number of orders -- xx for NIRPS

    --> Outputs:    - Attributes of Order objects:
                      1. W_raw (Wavelengths vectors)
                      2. I_raw (Time series of spectra)
                      3. blaze (Time series of blaze functions)
                      4. A_raw (Time series of telluric spectra computed from the DRS)
                      5. SNR (list of order S/N values for all observations)
                    - list_ord: upgraded list of orders
                    - airmass: airmass value for each observation
                    - bjd: time vector
                    - snr_mat: 2D matrix containing the S/N value for each observation and order (N_observation,N_order)
    """
    directory = os.listdir(repp)
    nam_t = []

    #Get only t.fits files (in case other files in the directory)
    for filename in directory:
        if "t.fits" in filename:
            nam_t.append(filename)
    
    nam_t     = sorted(nam_t)
    nobs      = len(nam_t)
    airmass   = np.zeros(nobs)
    bjd       = np.zeros(nobs)
    berv      = np.zeros(nobs)
    snr_mat   = np.zeros((nobs,nord))
    for nn in range(nobs):
        nmn          = repp + "/" + str(nam_t[nn])
        hdul_t       = fits.open(nmn)
        airmass[nn]  = (float(hdul_t[0].header["HIERARCH ESO TEL AIRM START"])
                        +float(hdul_t[0].header["HIERARCH ESO TEL AIRM END"]))/2 #airmass
        bjd[nn]      = float(hdul_t[1].header["BJD"])       #mid-exposure time [BJD-TDB]
        berv[nn]     = float(hdul_t[1].header["BERV"])      #Barycentric Earth Velocity [km/s]
        i            = np.array(hdul_t[1].data,dtype=float) # intensity spectrum
        w            = np.array(hdul_t[2].data,dtype=float) # wavelength vector [nm]
        bla          = np.array(hdul_t[3].data,dtype=float) # blaze vector
        atm          = np.array(hdul_t[4].data,dtype=float) # telluric spectrum
        ### Get S/N values
        for mm in range(nord):
            num = list_ord[mm].number
            if num < 10: key = "EXTSN00" + str(num)
            else: key = "EXTSN0" + str(num)
            sn  = float(hdul_t[1].header[key]) # S/N for each order
            list_ord[mm].SNR.append(sn)
            snr_mat[nn,mm] = sn
        hdul_t.close()
        ## Store Order's attributes
        for mm in range(nord):
            O = list_ord[mm]
            num = list_ord[mm].number
            O.W_raw.append(w[num])
            O.I_raw.append(i[num])
            O.blaze.append(bla[num])
            O.I_atm.append(atm[num])
    for mm in range(nord):
        O       = list_ord[mm]
        O.SNR   = np.array(O.SNR,dtype=float)
        O.W_raw = np.array(O.W_raw,dtype=float)
        O.I_raw = np.array(O.I_raw,dtype=float)
        O.blaze = np.array(O.blaze,dtype=float)
        O.I_atm = np.array(O.I_atm,dtype=float)
    
    return list_ord,airmass,bjd,berv,snr_mat


# -----------------------------------------------------------
# This function reads a time series of HARPS DRS (S2D) files
# It stores some of the relevant information into "Order" objects
# and returns time series relevant for the analysis
# The spectra to read must be stored in a dedicated repository
# -----------------------------------------------------------
def read_data_harps(repp,list_ord,nord):

    """
    --> Inputs:     - repp:      Path to the directory containing all the '.fits' files to read
                                 NOTE: files must be ordered in the chronologic order
                    - list_ord:  List of Order object
                    - nord:      Number of orders -- xx for NIRPS

    --> Outputs:    - Attributes of Order objects:
                      1. W_raw (Wavelengths vectors)
                      2. I_raw (Time series of spectra)
                      3. blaze (Time series of blaze functions)
                      4. A_raw (Time series of telluric spectra computed from the DRS)
                      5. SNR (list of order S/N values for all observations)
                    - list_ord: upgraded list of orders
                    - airmass: airmass value for each observation
                    - bjd: time vector
                    - snr_mat: 2D matrix containing the S/N value for each observation and order (N_observation,N_order)
    """
    directory = os.listdir(repp)
    nam_t = []
    nam_raw = []

    #Get only t.fits files (in case other files in the directory)
    for filename in directory:
        if "S2D_A.fits" in filename:
            nam_t.append(filename)
        if "S2D_BLAZE_A.fits" in filename:
            nam_raw.append(filename)
    
    nam_t     = sorted(nam_t)
    nam_raw     = sorted(nam_raw)
    nobs      = len(nam_t)
    airmass   = np.zeros(nobs)
    bjd       = np.zeros(nobs)
    berv      = np.zeros(nobs)
    snr_mat   = np.zeros((nobs,nord))
    for nn in range(nobs):
        nmn          = repp + "/" + str(nam_t[nn])
        hdul_t       = fits.open(nmn)
        nmn_raw      = repp + "/" + str(nam_raw[nn])
        hdul_raw     = fits.open(nmn_raw)
        airmass[nn]  = (float(hdul_t[0].header["HIERARCH ESO TEL AIRM START"])
                        +float(hdul_t[0].header["HIERARCH ESO TEL AIRM END"]))/2 #airmass
        bjd[nn]      = float(hdul_t[0].header["HIERARCH ESO QC BJD"])           # mid-exposure time [BJD-TDB] - check this is indeed in TDB
        berv[nn]     = float(hdul_t[0].header["HIERARCH ESO QC BERV"])          # Barycentric Earth Velocity [km/s]
        i            = np.array(hdul_raw[1].data,dtype=float)                   # intensity spectrum - raw
        w            = np.array(hdul_t[4].data,dtype=float) / 10                # wavelength vector [nm]
        bla          = np.array(hdul_raw[1].data/hdul_t[1].data,dtype=float) # blaze vector - computed from blaze corrected
        atm          = np.ones_like(i) # telluric spectrum - no telluric correction for the moment
        ### Get S/N values
        for mm in range(nord):
            num = list_ord[mm].number
            key = "HIERARCH ESO QC ORDER" + str(num+1) + " SNR"
            sn  = float(hdul_t[0].header[key]) # S/N for each order
            list_ord[mm].SNR.append(sn)
            snr_mat[nn,mm] = sn
        hdul_t.close()
        ## Store Order's attributes
        for mm in range(nord):
            O = list_ord[mm]
            num = list_ord[mm].number
            O.W_raw.append(w[num])
            O.I_raw.append(i[num])
            O.blaze.append(bla[num])
            O.I_atm.append(atm[num])
    for mm in range(nord):
        O       = list_ord[mm]
        O.SNR   = np.array(O.SNR,dtype=float)
        O.W_raw = np.array(O.W_raw,dtype=float)
        O.I_raw = np.array(O.I_raw,dtype=float)
        O.blaze = np.array(O.blaze,dtype=float)
        O.I_atm = np.array(O.I_atm,dtype=float)
    
    return list_ord,airmass,bjd,berv,snr_mat

def read_data_igrins(repp, list_ord, nord, fmt="PLP"):
    """
    --> Inputs:     - repp:      Path to the directory containing all the '.fits' files to read
                                 NOTE: files must be ordered in the chronologic order
                    - list_ord:  List of Order object
                    - nord:      Number of orders -- 54 for IGRINS
                    - fmt:       Format of the fits files -- allowed values:
                                 "PLP" for data reduced with IGRINS PLP (default)
                                 "GOA" for science quality data downloaded from Gemini Observatory Archive

    --> Outputs:    - Attributes of Order objects:
                      1. W_raw (Wavelengths vectors)            [microns]
                      2. I_raw (Time series of spectra)
                      3. blaze (Time series of blaze functions) - coming soon tm
                      4. A_raw (Time series of telluric spectra)  - coming soon tm
                      5. SNR (list of order S/N values for all observations)
                    - list_ord: upgraded list of orders
                    - airmass: airmass value for each observation
                    - bjd: time vector
                    - snr_mat: 2D matrix containing the S/N value for each observation and order (N_observation,N_order)
    """
    
    ### Get all data to read
    if fmt == 'PLP':
        #Indexes for Wavelength and spectrum
        datidx, wlidx = 0,1
        #Spectrum files
        specfilesH=sorted(glob.glob(repp+'*SDCH*spec.fits'))
        specfilesK=sorted(glob.glob(repp+'*SDCK*spec.fits'))

        #remove .sum.spec files for IGRINS-2 PLP
        specfilesH = [file for file in specfilesH if "sum" not in file]
        specfilesK = [file for file in specfilesK if "sum" not in file]

        #Variance files
        varfilesH=sorted(glob.glob(repp+'*SDCH*variance.fits'))
        varfilesK=sorted(glob.glob(repp+'*SDCK*variance.fits'))

        varfilesH = [file for file in varfilesH if "sum" not in file]
        varfilesK = [file for file in varfilesK if "sum" not in file]

        #SNR files
        snfilesH=sorted(glob.glob(repp+'*SDCH*sn.fits'))
        snfilesK=sorted(glob.glob(repp+'*SDCK*sn.fits'))

        snfilesH = [file for file in snfilesH if "sum" not in file]
        snfilesK = [file for file in snfilesK if "sum" not in file]

        assert len(specfilesH)==len(varfilesH)==len(specfilesK)==len(varfilesK), "Unequal no. of H and K files"
    elif fmt == "GOA":
        #Indexes for Wavelength and spectrum
        datidx, wlidx = 1, 3

        #GOA files contain all the info need
        specfilesH = sorted(glob.glob(repp+'*_H.spec.*'))
        specfilesK = sorted(glob.glob(repp+'*_K.spec.*'))

        assert len(specfilesH)==len(specfilesK), "Unequal no. of H and K files"
    else :
        print("Format '" + str(fmt) + "' unknown")

    ### Construct Wavelength grid from final spec file
    wlfile = fits.open(specfilesH[-1])  
    wlensH = wlfile[wlidx].data                  #wavelengths [microns]
    wlfile = fits.open(specfilesK[-1])
    wlensK = wlfile[wlidx].data                  #wavelengths [microns]
    wlens  = np.concatenate([wlensK,wlensH]) #descending order
    wlens  = wlens[::-1,:]                   #ascending order


    ### Initialisation
    ndet, npix = wlens.shape
    time_JD    = np.zeros(len(specfilesH))
    airmass    = np.zeros_like(time_JD)
    humidity   = np.zeros_like(time_JD)
    data_RAW   = np.zeros((ndet,len(specfilesH),npix))
    data_var   = np.zeros_like(data_RAW)
    data_snr   = np.zeros_like(data_RAW)


    ### open and read files -- arrange data into spectral matrix
    for ifile in range(len(specfilesH)):
        #H-band
        hdu_listH = fits.open(specfilesH[ifile])
        image_dataH = hdu_listH[datidx].data

        #K-band
        hdu_listK = fits.open(specfilesK[ifile])
        image_dataK = hdu_listK[datidx].data

        #SNR and variances
        if fmt == "PLP":
            #H-band
            hdu_list   = fits.open(varfilesH[ifile])
            var_H = hdu_list[0].data
            hdu_list   = fits.open(snfilesH[ifile])
            snr_H  = hdu_list[0].data

            #K-band
            hdu_list = fits.open(varfilesK[ifile])
            var_K = hdu_list[0].data
            hdu_list   = fits.open(snfilesK[ifile])
            snr_K  = hdu_list[0].data
        elif fmt == "GOA":
            #H-band
            var_H = hdu_listH[2].data
            snr_H = image_dataH / np.sqrt(var_H)
            #K-band
            var_K = hdu_listK[2].data
            snr_K = image_dataK / np.sqrt(var_K)

        #Date, airmass, humidity info
        hdr = hdu_listK[0].header
        instrument = hdr["INSTRUME"]
        if instrument=="IGRINS":
            time_JD[ifile] = 0.5 * (hdr["JD-OBS"] + hdr["JD-END"]) # mid-exposure time, JD UTC
            humidity[ifile] = hdr['HUMIDITY']                      # humidity
        else:
            frame_start = Time(hdr['UTSTART'], format="isot", scale="utc")
            frame_end = Time(hdr['UTEND'], format="isot", scale="utc")
            time_JD[ifile] = 0.5 * (frame_start.jd + frame_end.jd)

        airmass[ifile] = 0.5 * (hdr['AMSTART']+hdr['AMEND'])     # average airmass

        #concatenating K and H spectra
        data = np.concatenate([image_dataK,image_dataH])
        var  = np.concatenate([var_K,var_H])
        snr   = np.concatenate([snr_K,snr_H])
        data_RAW[:,ifile,:] = data # master matrix
        data_var[:,ifile,:] = var
        data_snr[:,ifile,:]  = snr

    #invert arrays - ascending order in wavelength
    data_RAW = data_RAW[::-1,:,:]
    data_var = data_var[::-1,:,:]
    data_snr  = data_snr[::-1,:,:]

    # compute barycentric correction
    print("\ncompute barycentric velocities")
    location = wlfile[0].header['TELESCOP']                                 #Observatory name (Gemini-N/S)
    if location == 'Gemini-North':
        location = 'gemini_north'
        RA, DEC = wlfile[0].header['RA'], wlfile[0].header['DEC']           #Object RA, DEC [degrees]
    elif location == 'Gemini-South':
        location='gemini_south'
        RA, DEC = wlfile[0].header['OBJRA'], wlfile[0].header['OBJDEC']      #Object RA, DEC [degrees]
    try:
        pmRA  = float(wlfile[0].header["PMRA"])*1e3 
        pmDEC = float(wlfile[0].header["PMDEC"])*1e3    #proper motion if any [mas/y]
    except:
        pmRA, pmDEC = 0, 0                              #if none set to 0

    time_JD = Time(time_JD, format="jd", scale="utc")
    bjd, warning, status = JDUTC_to_BJDTDB(JDUTC = time_JD, ra=RA, dec=DEC, pmra=pmRA, pmdec=pmDEC, obsname=location)
    if status==2:
        print(warning)
    berv, warning, status = get_BC_vel(JDUTC = time_JD, ra=RA, dec=DEC, pmra=pmRA, pmdec=pmDEC, obsname=location)
    if status==2:
        print(warning)
    bjd = np.array(bjd)                 #Observation time in BJD TDB
    berv = np.array(berv) / 1e3         #Barycentric Earth Velocity [km/s]
    
    #average SNR over each order
    snr_mat   = np.nanmean(data_snr,axis=2)
    #duplicate wavelength array to each integration sequence
    wl_raw    = np.zeros_like(data_RAW)
    for nn in range(len(bjd)):
        wl_raw[:,nn,:] = wlens
    
    
    for mm in range(nord):
        O       = list_ord[mm]
        blaze   = np.ones_like(data_RAW[mm])
        atm     = np.ones_like(data_RAW[mm])
        O.SNR   = np.array(snr_mat[mm],dtype=float)
        O.W_raw = np.array(wl_raw[mm],dtype=float)
        O.I_raw = np.array(data_RAW[mm],dtype=float)
        O.blaze = np.array(blaze,dtype=float)
        O.I_atm = np.array(atm,dtype=float)
    
    return list_ord,airmass,bjd,berv,snr_mat.T

# -----------------------------------------------------------
# Get transit window -- requires batman python module
# Uncomment lines below to use batman module to compute transit flux
# See information in https://lweb.cfa.harvard.edu/~lkreidberg/batman/
# -----------------------------------------------------------
def compute_transit(Rp,Rs,ip,T0,ap,Porb,ep,wp,limb_dark,uh,T_obs,ttype="primary",T_eclipse=0,fp=0):

    """
    --> Inputs:     - Rp:        Planetary radius
                    - Rs:        Stellar radius (same unit as Rp)
                    - ip:        Transit inclination [deg]
                    - T0:        Mid-transit time (same unit as T_obs -- here: bjd)
                    - ap:        Semi-major-axis [Stellar radius]
                    - Porb:      Planet orbital period (same unit as T_obs)
                    - ep:        Eccentricity of the planet orbit
                    - wp:        Argument of the periapsis for the planet orbit [deg] 
                    - limb_dark: Type of limb-darkening model: "linear", "quadratic", "nonlinear" see https://lweb.cfa.harvard.edu/~lkreidberg/batman/
                    - uh:        Limb-darkening coefficients matching the type of model and in the SPIRou band (Typically H or K)
                    - T_obs:     Time vector

    --> Outputs:    - flux:      Relative transit flux (1 outside transit) 
    """
#

    params           = batman.TransitParams()
    params.rp        = Rp/Rs                       
    params.inc       = ip
    params.t0        = T0
    params.a         = ap
    params.per       = Porb
    params.ecc       = ep
    params.w         = wp         
    params.limb_dark = limb_dark
    params.u         = uh
    params.fp        = fp
    params.t_secondary = T_eclipse
    if ttype=="secondary":
        bat              = batman.TransitModel(params,T_obs,transittype="secondary")
        flux             = bat.light_curve(params)
    else:
        bat              = batman.TransitModel(params,T_obs)
        flux             = bat.light_curve(params)
    return flux


# -----------------------------------------------------------
# Compute RV signature induced by the planet on its host star
# Assuming circular orbit for the planet
# -----------------------------------------------------------
def get_rvs(t,k,p,t0):

    """
    --> Inputs:     - t:   Time vector
                    - k:   Semi-amplitude of the planet-induced RV signal on the host star
                    - p:   Planet orbital period
                    - t0:  Planet mid-transit time

    --> Outputs:    - Planet-induced RV signature for the input time values 
    """

    return - k*np.sin(2.*np.pi/p * (t-t0))

# -----------------------------------------------------------
# Compute planet RV in the stellar rest frame
# -----------------------------------------------------------
def rvp(phase,k,v):
    """
    --> Inputs:     - phase: Planet orbital phase (T-T_obs)/Porb
                    - k:     Semi-amplitude of the planet RV
                    - v:     Planet RV at mid-transit

    --> Outputs:    - Planet RV for the input orbital phases
    """
    return k*np.sin(2.*np.pi*phase)+v
 
#### Main class -- Order
    
class Order:


    def __init__(self,numb):
        
        ### Generic information
        self.number       = numb    # Order number (in absolute unit -- 79: bluest; 31: reddest)
        self.W_mean       = 0.0     # Mean order wavelength
        self.SNR          = []      # DRS-computed S/N at the center of the order
        self.SNR_mes      = []      # Empirical estimate of the SNR before PCA
        self.SNR_mes_pca  = []      # Empirical estimate of the SNR after PCA

        ### Raw data information
        self.W_raw    = []      # Wavelength vectors for the raw observations - 2D matrix (time-wavelength)
        self.I_raw    = []      # Time series of observed spectra from the SPIRou DRS - 2D matrix (time-wavelength)                
        self.I_atm    = []      # Time series of Earth atmosphere spectra computed from the observations using Artigau+2014 method - DRS-provided 2D matrix (time-wavelength)
        self.blaze    = []      # Time series of blaze functions - 2D matrix (time-wavelength)

        ### Data reduction information
        self.W_fin    = []      # Final wavelength vector in the Geocentric frame
        self.W_bary   = []      # Final wavelength vector in the stellar rest frame
        self.I_fin    = []      # Reduced flux matrix before PCA cleaning
        self.I_pca    = []      # Reduced flux matrix after PCA cleaning
        self.proj     = []      # Reduced projected flux matrix after Gibson+21 cleaning



    # -----------------------------------------------------------
    # Pre-process of the DRS-provided spectra:
    # 1. Blaze normalization process
    # 2. Remove NaNs from each spectrum and convert sequences of
    #    spectra into np.array square matrices
    # -----------------------------------------------------------
    def remove_nan(self):
        
        """
        --> Inputs:      - Order object
        
        --> Outputs:     - Boolean: 1 --> order empty as NaNs everywhere; 0 otherwise 
        """

        ### Remove blaze
        I_bl = self.I_raw/self.blaze
        W_bl = self.W_raw

        
        ### Spot the NaNs:
        ### In "*t.fits" files, regions of high telluric absorptions are replaced by NaNs
        ### as no precise estimation of the flux could be carried out
        ### Here we build a vector 'ind' stroring the position of the NaNs in every spectrum
        ind   = []
        for nn in range(len(I_bl)):
            i = np.where((np.isfinite(I_bl[nn])==True)&(np.isfinite(W_bl[nn])==True))[0]
            ind.append(i)
        r  = np.array(list(set.intersection(*map(set,ind))),dtype=int)
        r  = np.sort(np.unique(r))

        ### remove the NaNs
        I_ini = []
        W_ini = []
        B_ini = []
        A_ini = []
        for nn in range(len(I_bl)):
            I_ini.append(I_bl[nn,r])
            W_ini.append(self.W_raw[nn,r])
            A_ini.append(self.I_atm[nn,r])
            B_ini.append(self.blaze[nn,r])

        ### Convert into 2D array object
        self.I_raw  = np.array(I_ini,dtype=float)    
        self.W_raw  = np.array(W_ini,dtype=float)[0]   
        self.I_atm  = np.array(A_ini,dtype=float)
        self.B_raw  = np.array(B_ini,dtype=float) 
        self.W_mean = self.W_raw.mean()   ### Compute mean of the actual observations

        ### Remove the order if it contains only NaNs
        if len(self.I_raw[0]) == 0:
            tx = "\nOrder " + str(self.number) + " is empty and thus removed from the analysis"
            print(tx)
            return 1
        else:
            return 0




    # -----------------------------------------------------------
    # Add a synthetic planet signature to an input sequence of spectra
    # 1. Compute a synthetic sequence of spectra by:
    #    (i)  shifting a planetary template (normalized transit depth)
    #         according to the planet RV
    #    (ii) Weight by a transit window (1 at mid-transit, 0 outside
    #         of transit
    # 2. Interpolate each synthetic planet spectrum and multiply an
    #    input sequence of spectra by the synthetic sequence of transit
    #    depths
    # -----------------------------------------------------------
    def add_planet(self,type_obs,Wm,Im,window,planet_speed,Vc,ampl=1.0,pixel=np.linspace(-1.13,1.13,11)):

        """
        --> Inputs:      - Order object
                         - type_obs (str) emission or tranmission
                         - Wm:      Wavelength vector of the planet atmosphere template
                         - Im:      Template of wavelength-dependent transit depth (i.e., model) 
                         - window:  Transit window
                         - planet_speed : speed of the planet along the orbit
                         - Vc:      Velocimetric correction to move from Geocentric frame to stellar rest frame
                                    Typically: Vc = Stellar systemic vel. + planet-signature RV - Barycentric Earth RV 
                         - ampl:    Amplification factor: amplify the injected planetary signal

        --> Outputs:     - self.I_raw_pl
        """
        self.Wm    = Wm
        self.Im    = ampl*(Im-np.max(Im))+np.max(Im)        
        self.I_syn = np.zeros(self.I_raw.shape)

        Imm  = self.Im#/np.min(flux)
        if type_obs =="transmission":
            tdepth_interp     = PchipInterpolator(self.Wm,Imm)  # Interpolate model
        else:
            flux_interp = PchipInterpolator(self.Wm,Imm)  
            
        if type_obs =="transmission":
            for nn in range(len(self.I_raw)): # For each observation date
                if window[nn] != 0.0:
                    I_ttt = np.zeros(len(self.W_raw))
                    
                    # Shift model in the Geocentric frame
                    for pp in pixel: I_ttt += tdepth_interp(self.W_raw/(1.0+((planet_speed[nn]+Vc[nn]+pp)/(c0/1000.))))
                    self.I_syn[nn] = I_ttt/len(pixel)*window[nn]
                
        else:
            for nn in range(len(self.I_raw)): # For each observation date
                I_ttt = np.zeros(len(self.W_raw))
                # Shift model in the Geocentric frame
                for pp in pixel: I_ttt += flux_interp(self.W_raw/(1.0+((planet_speed[nn]+Vc[nn]+pp)/(c0/1000.))))
                self.I_syn[nn]  = I_ttt/len(pixel)
                
                

 

