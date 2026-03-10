import numpy as np
from numpy import ma

import sys
import os
import time

from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit
import scipy.stats as stats


#import astropy
#from astropy.io import fits
#from astropy.io import fits
#from astropy.stats import sigma_clip
#from astropy.modeling import models, fitting, polynomial

import matplotlib.pyplot as plt

import convolve as conv


def strided_app(a, L, S ):  
    # This function is an updated moving median 
    # used in the normalization process
    # Window len = L, Stride len/stepsize = S
    nrows = ((a.size-L)//S)+1
    n = a.strides[0]
    return np.lib.stride_tricks.as_strided(a, shape=(nrows,L), strides=(S*n,n))

class reduced:

    def __init__(self,nord,W,Rp,R_s,emission=False):

        self.ord  = nord  ## List of orders to be used -- after preselection

        self.Wm0 = W
        self.Wm = W
        self.Rp = Rp
        self.R_s = R_s
        self.emission = emission

    def convolve(self,rot_speed,superrot):
        #although it works really well, we have a mean issue due to the numerical
        #limits. Be careful with the mean !

         R = self.Rp
         init_mean = np.mean(R**2)
         rot_wl,rot_Rp = conv.rotate(R,self.Wm,rot_speed,superrot)
         final_mean = np.mean(rot_Rp**2)
         rot_Rp = rot_Rp**2-final_mean+init_mean
         rot_Rp = np.sqrt(rot_Rp)
         self.Wm = rot_wl
         self.Rp = rot_Rp

    def rotate(self, rot_speed):
        #Convolve with rotation kernel from rot_int_cmj
        self.Rp = conv.rot_int_cmj(self.Wm, self.Rp, rot_speed)
         
    def normalize(self,nf=501):     
        #renor
        out = np.zeros((2,np.shape(self.Wm)[0]))
        out[0] = self.Wm
        if self.emission:
            out[1] = self.Rp
        else:
            out[1] = (self.Rp/self.R_s)**2


        win = np.percentile(strided_app(out[1],nf,1),0.5, axis=-1)
        try :
            out_final = np.zeros((2,len(out[1])-nf+1))
        except:
            print ("Removed order : ",self.ord)
            exit
        out_final[0]  = out[0][int((nf-1)/2):-int((nf-1)/2)]
        out_final[1] = win - out[1][int((nf-1)/2):-int((nf-1)/2)]

        self.Wm = out_final[0]
        if self.emission:
            self.Rp = out_final[1]*-1
        else:
            self.Rp = out_final[1]


        



def prepare(model_dic,R_s,orderstot,winds=False,rot_speed=0.0,superrot=0.0,emission=False):

    models = []
    flux_star = []
    wavelength = model_dic["wavelength_nm"]
    radius = model_dic["radius_transm"]
    star_flux = model_dic["star_flux"]
    for i in range(len(orderstot)):
        ### First we have to select only the portion of  the model that is of interest for us
        # limits = [self.lambdas[no][0]*0.995,self.lambdas[no][1]*1.005]
        no = orderstot[i]
        M  = reduced(no,wavelength[i],radius[i],R_s,emission)
        
        if winds:
            M.convolve(rot_speed,superrot)
        elif rot_speed>0:
            M.rotate(rot_speed)
        
        M.normalize()
        if emission:
            flux_star.append([wavelength[i],star_flux[i]])

        models.append(M)
#### LRS if there is any need
    # wavelength_LR_microns = model_dic["wavelength_LR_microns"]
    # tdepth_LR = model_dic["tdepth_LR"]
    # fLR = interpolate.CubicSpline(wavelength_LR_microns,tdepth_LR)


    return {
       "models": models,
       "flux_star":flux_star
       # "interp_LR": fLR, 
       }











