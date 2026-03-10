#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 07/2022

@author: Baptiste & Florian
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.ticker as ticker
import matplotlib.patheffects as PathEffects
import warnings
import matplotlib.cbook

from scipy.optimize import minimize



def plot_transit(T_obs, phase, T0, flux, airmass, Vc, snr_mat,figure_name,type_obs):
        if type_obs == "transmission":
            TT     = 24.*(T_obs - T0)
        else:
            TT     = phase
        ypad   = 15  # pad of the y label
        plt.figure(figsize=(15,12))
        # Transit flux
        ax  = plt.subplot(411)
        ax.plot(TT,flux,"-+r",label="planet")
        plt.legend(loc=3,fontsize=16)
        ax.set_ylabel("Transit curve\n", labelpad=ypad)
        # Airmass
        ax = plt.subplot(412)
        plt.plot(TT,airmass,"-k")
        ax.set_ylabel("Airmass\n", labelpad=ypad)
        # RV correction between Geocentric frame and stellar rest frame
        ax = plt.subplot(413)
        plt.plot(TT,Vc,"-k")
        ax.set_ylabel("RV correction\n[km/s]", labelpad=ypad)
        # Maximum S/N
        ax = plt.subplot(414)
        plt.plot(TT,np.max(snr_mat,axis=1),"+k")
        plt.axhline(np.mean(np.max(snr_mat,axis=1)),ls="--",color="gray")
        if type_obs == "transmission":
            plt.xlabel("Time wrt transit [h]")
        else:
            plt.xlabel("Orbital Phase")
        ax.set_ylabel("Peak S/N\n", labelpad=ypad)
        plt.subplots_adjust(hspace=0.02)
        plt.savefig(figure_name,bbox_inches="tight")
        plt.close()


# -----------------------------------------------------------
# Compare the dispersion at the center of each spectrum (in 
# each order) to the photon noise provided by the SPIRou DRS
# -----------------------------------------------------------
def plot_spectrum_dispersion(lord,nam_fig):

    """
    --> Inputs:     - lord: list of Order objects

    --> Outputs:    - Plot displayed

     """

    # Initialization
    rms_sp     = np.zeros(len(lord))
    rms_sp_s   = np.zeros(len(lord))
    rms_drs    = np.zeros(len(lord))
    rms_drs_s  = np.zeros(len(lord))
    rms_pca    = np.zeros(len(lord))
    rms_pca_s  = np.zeros(len(lord))    
    wmean      = np.zeros(len(lord))
    npca       = np.zeros(len(lord))
    LO         = np.zeros(len(lord),dtype=int)

    for kk in range(len(lord)):
        O              = lord[kk]
        disp_mes       = 1./O.SNR_mes
        disp_drs       = 1./O.SNR
        disp_pca       = 1./O.SNR_mes_pca
        rms_sp[kk]     = np.mean(disp_mes)
        rms_sp_s[kk]   = np.std(disp_mes)
        rms_drs[kk]    = np.mean(disp_drs)
        rms_drs_s[kk]  = np.std(disp_drs)
        rms_pca[kk]    = np.mean(disp_pca)
        rms_pca_s[kk]  = np.std(disp_pca)
        npca[kk]       = O.n_com
        wmean[kk]      = O.W_mean
        LO[kk]         = O.number

    # Compute wavelength-order number correspondance
    WW,LO_pred,LO_predt = fit_order_wave(LO,wmean)
    plt.figure(figsize=(12,5))
    ax = plt.subplot(111)
    ax3 = ax.twinx()
    ax3.plot(LO, npca, marker="o", c="r", alpha=0.6)
    ax.plot([],[], marker="o", c="r", alpha=0.6, label="PCA components")
    ax.errorbar(LO,rms_sp,rms_sp_s,fmt="*",color="k",label="Reduced data",capsize=10.0,ms=10.)
    ax.errorbar(LO,rms_pca,rms_pca_s,fmt="^",color="g",label="After PCA",capsize=10.0,ms=7.5)
    ax.errorbar(LO,rms_drs,rms_drs_s,fmt="o",color="m",label="DRS",capsize=8.0)
    
    ax.legend(ncol=2)
    ax2 = ax.twiny()
    ax2.scatter(wmean,rms_sp,alpha=0)
    ax2.set_xticks(LO_pred)
    ax2.xaxis.set_major_locator(MultipleLocator(300))
    ax2.set_xlabel("Wavelength [nm]")
    ax2.set_xticklabels(WW)
    ax2.xaxis.set_minor_locator(MultipleLocator(100))
    ax.set_xlim(np.min(LO)-1,np.max(LO)+1)
    ax2.set_xlim(np.min(LO)-1,np.max(LO)+1)
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_ylabel("Spectrum dispersion")
    ax.set_xlabel("Order number")
    ax.set_yscale("log")
    ax3.set_ylim(0,np.max(npca)+1)
    ax3.set_ylabel("Components")
    plt.subplots_adjust(wspace=0.5,hspace = 0.)
    plt.savefig(nam_fig,bbox_inches="tight")
    plt.close()
    




def plot_reduction(phase,W0,W1,I1,W2,I2,W3,I3,W4,I4,Am,lab=["1","2","3","4"],filenam="reduc.png",lmin=-1,lmax=-1):

    ### Show the evolution of the sequence of spectra for a given order
    cmap = "viridis" # Another 'fancy' color map?
    I1  = I1/np.mean(I1)
    I1 -= np.mean(I1)

    #plt.figure()
    fig, axis = plt.subplots(6,2,height_ratios=[1,0.08,1,1,1,1],width_ratios=(30,1),figsize=(15,17))

    ax=axis[0,0]
    ax.plot(W0, Am, c="k")
    ax.set_ylabel("Transmission ")
    ax.tick_params(axis="x",labelbottom=False)

    #Trick for proper spacing between plots
    axis[0,1].set_visible(False)
    axis[1,0].set_visible(False)
    axis[1,1].set_visible(False)

    ###############################        
    ax   = axis[2,0]
    width = ax.get_window_extent().transformed(ax.get_figure().dpi_scale_trans.inverted()).width*ax.get_figure().dpi
    line_width = width /len(W0) / 2
    X,Y  = np.meshgrid(W1,phase)
    Z    = I1
    zmin = I1.mean() - 3.*I1.std()
    zmax = I1.mean() + 3.*I1.std()
    c    = ax.pcolor(X,Y,Z,cmap=cmap,vmin=zmin,vmax=zmax) 
    plt.colorbar(c,ax=ax,cax=axis[2,1])
    wl_rem = [value for value in W0 if value not in W1]
    ax.vlines(wl_rem, min(phase), max(phase), color="red", lw=line_width)
    ax.set_ylabel("Phase")  
    ax.set_xticks([])        
    if lmin>-1: ax.set_xlim(lmin,lmax)
    else: ax.set_xlim(W0.min(),W0.max())

    tx = lab[0]                           
    ax.text(np.min(W4)+1.5,0.9*np.max(phase),tx,color="w",fontsize=20,fontweight="bold")

    ###############################  
    ax   = axis[3,0]
    X,Y  = np.meshgrid(W2,phase)
    Z    = I2
    zmin = I2.mean() - 3.*I2.std()
    zmax = I2.mean() + 3.*I2.std()
    c    = ax.pcolor(X,Y,Z,cmap=cmap,vmin=zmin,vmax=zmax)
    wl_rem = [value for value in W0 if value not in W2]
    ax.vlines(wl_rem, min(phase), max(phase), color="red", lw=line_width)   
    ax.set_ylabel("Phase")  
    ax.set_xticks([])        
    if lmin>-1: ax.set_xlim(lmin,lmax)
    else: ax.set_xlim(W0.min(),W0.max())
    
    plt.colorbar(c,ax=ax,cax=axis[3,1])

    tx = lab[1]                           
    ax.text(np.min(W4)+1.5,0.9*np.max(phase),tx,color="w",fontsize=20,fontweight="bold")


    ###############################  
    ax   = axis[4,0]
    X,Y  = np.meshgrid(W3,phase)
    Z    = I3
    zmin = I3.mean() - 3.*I3.std()
    zmax = I3.mean() + 3.*I3.std()
    c    = ax.pcolor(X,Y,Z,cmap=cmap,vmin=zmin,vmax=zmax)  
    wl_rem = [value for value in W0 if value not in W3]
    ax.vlines(wl_rem, min(phase), max(phase), color="red", lw=line_width)
    ax.set_ylabel("Phase")  
    ax.set_xticks([])        
    if lmin>-1: ax.set_xlim(lmin,lmax)
    else: ax.set_xlim(W0.min(),W0.max())
    
    plt.colorbar(c,ax=ax,cax=axis[4,1])
    tx = lab[2]                           
    ax.text(np.min(W4)+1.5,0.9*np.max(phase),tx,color="w",fontsize=20,fontweight="bold")


    ###############################  
    ax   = axis[5,0]
    X,Y  = np.meshgrid(W4,phase)
    Z    = I4
    zmin = I4.mean() - 3.*I4.std()
    zmax = I4.mean() + 3.*I4.std()
    c    = ax.pcolor(X,Y,Z,cmap=cmap,vmin=zmin,vmax=zmax)  
    wl_rem = [value for value in W0 if value not in W4]
    ax.vlines(wl_rem, min(phase), max(phase), color="red", lw=line_width)
    ax.set_ylabel("Phase")  
    if lmin>-1: ax.set_xlim(lmin,lmax)
    else: ax.set_xlim(W0.min(),W0.max())
    
    plt.colorbar(c,ax=ax,cax=axis[5,1])
    ax.set_xlabel("Wavelength [nm]")

    tx = lab[3]                           
    ax.text(np.min(W4)+1.5,0.9*np.max(phase),tx,color="w",fontsize=20,fontweight="bold")
    plt.subplots_adjust(hspace=0.02,wspace=0.02)
    
    plt.savefig(filenam,bbox_inches="tight")
    #plt.show()
    plt.close()


def plot_reduction_tot(list_ord, phase, phase_rem, filename):
    phase_mask = np.ones_like(phase, dtype=bool)
    for pp in phase_rem:
        phase_mask[pp] = 0
    
    N_raw, N_pca = 0, 0
    W_raw, I_raw, W_pca, I_pca = [], [], [], []
    for nn in range(len(list_ord)):
        N_raw += len(list_ord[nn].W_raw)
        N_pca += len(list_ord[nn].W_fin)
        W_raw.append(list_ord[nn].W_raw)
        I_raw.append(list_ord[nn].I_raw0/np.mean(list_ord[nn].I_raw0))
        W_pca.append(list_ord[nn].W_fin)
        I_pca.append(list_ord[nn].I_pca)
        for pp in phase_rem:
            I_pca[nn] = np.insert(I_pca[nn], pp, np.zeros_like(I_pca[nn][0,:]),axis=0)
    
    std1, std2 = 0, 0
    for nn in range(len(list_ord)):
        std1 += np.std(I_raw[nn][phase_mask]) * len(W_raw[nn]) / N_raw
        std2 += np.std(I_pca[nn][phase_mask]) * len(W_pca[nn]) / N_pca

    zmin1,zmax1 = -3.*std1, 3*std1
    zmin2,zmax2 = -3.*std2, 3*std2
    fig, ((ax1,ax3),(ax2,ax4)) = plt.subplots(2,2, figsize=(15,8),width_ratios=(30,1))
    width = ax1.get_window_extent().transformed(ax1.get_figure().dpi_scale_trans.inverted()).width*ax1.get_figure().dpi
    height = ax1.get_window_extent().transformed(ax1.get_figure().dpi_scale_trans.inverted()).height*ax1.get_figure().dpi
    for nn in range(len(list_ord)):
        X1,Y1 = np.meshgrid(W_raw[nn],phase)
        Z1    = I_raw[nn] - I_raw[nn][phase_mask].mean() 
        c1    = ax1.pcolor(X1,Y1,Z1,cmap="viridis",vmin=zmin1,vmax=zmax1)
    
        X2,Y2 = np.meshgrid(W_pca[nn],phase)
        Z2    = I_pca[nn] - I_pca[nn][phase_mask].mean() 
        c2    = ax2.pcolor(X2,Y2,Z2,cmap="viridis",vmin=zmin2,vmax=zmax2)
    
    wmin, wmax = min(np.hstack((W_raw[0],W_raw[-1]))), max(np.hstack((W_raw[0],W_raw[-1])))
    for pp in phase_rem:
        ax2.hlines(phase[pp], xmin=wmin, xmax=wmax, colors="red", lw=height/4/len(phase))

    plt.colorbar(c1,ax=ax1,cax=ax3,label="Normalised \n Flux")
    plt.colorbar(c2,ax=ax2,cax=ax4,label="Normalised \n Flux")
    
    ax1.set_xlim((wmin,wmax))
    ax2.set_xlim((wmin,wmax))

    txt1 = ax1.text(wmin+20,0.9*np.max(phase),"Before Post-Processing",color="w",fontsize=40,fontweight="bold")
    txt1.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='k')])
    txt2 = ax2.text(wmin+20,0.9*np.max(phase),"After Post-Processing",color="w",fontsize=40,fontweight="bold")
    txt2.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='k')])

    plt.subplots_adjust(hspace=0.02,wspace=0.02)
    ax1.tick_params(axis="x",labelbottom=False)

    ax1.set_ylabel("Phase")
    ax2.set_ylabel("Phase")
    ax2.set_xlabel("Wavelength (nm)")

    plt.savefig(filename,bbox_inches="tight")
    plt.close()

# -----------------------------------------------------------
# Compute Order to mean wavelength equivalence 
# Usage: Plot order number as X-axis and mean wavelengths as Y axis
# In practice: fits an hyperbola between order nb and mean wavelength
# See function plots.plot_orders for more information
# -----------------------------------------------------------
def fit_order_wave(LO,wm_fin):

    """
    --> Inputs:     - LO: list of order numbers
                    - wm_fin: list of the mean wavelengths corresponding to LO

    --> Outputs:    - WW: Wavelength ticks for the plot
                    - LO_pred: order numbers corresponding to WW
                    - LO_predt: densely-sampled list of orders for minor ticks locators
    """

    par0    = np.array([100000,200.0],dtype=float) 
    res     = minimize(crit_hyp,par0,args=(LO,wm_fin))
    p_best  = res.x 
    LO_tot  = np.arange(np.min(LO),np.max(LO)+1)
    pp      = hyp(p_best,LO_tot)
    WWT      = np.linspace(np.max(wm_fin),np.min(wm_fin),16)
    WW       = np.array([2400.0,2100,1800,1500,1200,1000],dtype=int)
    LO_predt = hyp_inv(p_best,WWT)
    LO_pred  = hyp_inv(p_best,WW) 
    return WW,LO_pred,LO_predt


# -----------------------------------------------------------
# Simple hyperbola
# -----------------------------------------------------------
def hyp(par,xx):
    return par[0]/xx + par[1]

# --------------------
# Simple inverse hyperbola
# -----------------------------------------------------------
def hyp_inv(par,yy):
    return par[0]/(yy-par[1])

# ----------------------------------# -----------------------------------------------------------
# Return least-square difference between a hyperbola for 'par' 
# parameters and data yy.
# xx is the X-axis vector 
# -----------------------------------------------------------
def crit_hyp(par,xx,yy):
    y_pred = hyp(par,xx)
    return np.sum((yy-y_pred)**(2))  
    