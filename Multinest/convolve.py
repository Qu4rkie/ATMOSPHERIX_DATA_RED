#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 13:31:37 2021

@author: florian

"""

#Everything is in SI, of course !
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal


def rotate(R,wl,vrot,superrot,angle_super=25.0*np.pi/180,sigma=1600.0) :

    if vrot<100: #we are not able to see the difference anyway , and that prevents from creating exceptions in 
    #the code
        vrot = 100
    
    c0 = 299792458.0

    #we take a kernel of size 50 000 km/s, it is exagerated but is safe
    vlim = 50000.0
    nv = 5000
    v = np.linspace(-vlim,vlim,nv)
    dv = v[1]-v[0]
    #we interpolate the model onto a regularly spaced speed array
    w0 = np.mean(wl)
    speed = c0*(w0/wl-1)
    speed_int = np.arange(0.995*np.min(speed),0.995*np.max(speed),step=dv)
    fmod = interpolate.interp1d(speed,R)
    mod_int = fmod(speed_int)

    
    #we prepare convolution by defining an appropriate speed array
    n = int(2*vrot/dv) 
    fraclim = dv*n/2/vrot
    vrot_array = np.arange(-vrot*fraclim,vrot,step=dv) #a trick to ensure symmetry in this array
    first_conv=  1/np.sqrt(1-(vrot_array/vrot)**2) #rotation kernel
    second_conv = np.exp(-v**2/2/sigma**2) #instrumental kernel
#    third_conv = np.zeros(nv)
#    for i in range(nv):
#        if np.abs(v[i])<3000:
#            third_conv[i] = 1./6000.
    if superrot>100.0: #if lower, we won't see the difference anyway
        #limits of convolution for superrotating parts
        cos1 = np.cos(angle_super)
        pos1 = np.where(vrot_array/vrot<-cos1)[0][-1]
        first_conv_1 = np.zeros(n+1) #to ensure same size, we fill the convoluant with zeros
        first_conv_2 = np.zeros(n+1) 
        first_conv_1[:pos1+1] = 1/np.sqrt(1-((vrot_array[:pos1+1])/vrot)**2)
        cos2 = -cos1
        pos2 = np.where(vrot_array/vrot>-cos2)[0][0]
        first_conv_2[pos2:] = 1/np.sqrt(1-((vrot_array[pos2:])/vrot)**2)
        nsuper_conv = int(superrot/dv)
        if nsuper_conv%2 ==1:
            nsuper_conv+=1
        nsuper_half = int(nsuper_conv/2) #in the worst case this leads to an error of 5 m/s ... that's allright


        conv1_mod = 1/vrot/np.pi*((vrot-cos1*(vrot))/(pos1+1)*signal.oaconvolve(mod_int[2*nsuper_conv:]**2, first_conv_1,mode="same") + \
                                  (2*cos1*vrot)/(len(vrot_array)-2*pos1-2)*signal.oaconvolve(mod_int[2*nsuper_half:-2*nsuper_half]**2, first_conv[pos1+1:pos2],mode="same") + \
                                  (vrot-cos1*(vrot))/(pos1+1)*signal.oaconvolve(mod_int[:-2*nsuper_conv]**2, first_conv_2,mode="same"))
        conv2_mod = 1/sigma/np.sqrt(2*np.pi)*2*vlim/nv*signal.oaconvolve(conv1_mod,second_conv,mode="same")  
#        conv3_mod = 2*vlim/nv*signal.oaconvolve(conv2_mod,third_conv,mode="same") 
        convtot_mod = np.sqrt(conv2_mod)

        wl_int = w0/(1+speed_int[nsuper_conv:-nsuper_conv]/c0)

        
    else:
        conv1_mod = 1/vrot/np.pi*((np.max(vrot_array)-np.min(vrot_array))/len(vrot_array)*signal.oaconvolve(mod_int**2, first_conv,mode="same"))#+(0.0447251+0.051)*vrot*np.mean(mod_int**2))
        conv2_mod = 1/sigma/np.sqrt(2*np.pi)*2*vlim/nv*signal.oaconvolve(conv1_mod,second_conv,mode="same")  
#        conv3_mod = 2*vlim/nv*signal.oaconvolve(conv2_mod,third_conv,mode="same")
        convtot_mod = np.sqrt(conv2_mod)

        wl_int = w0/(1+speed_int/c0)


    #the model is largely oversampled and we don't want that, it is going to kill the calculation time
    diff = len(convtot_mod)/len(wl)
    spacing = max(1,round(5.*diff/6))


    return (wl_int[n+nv:-n-nv:spacing],convtot_mod[n+nv:-n-nv:spacing])
#np.savetxt("/home/florian/Bureau/Atmosphere_SPIRou/Transit_2D/Convolution/test_model/templates/allday_superrot2/lambdasallday_superrot2.txt",c0/freq_tot[::100][::-1]*1e9)
#np.savetxt("/home/florian/Bureau/Atmosphere_SPIRou/Transit_2D/Convolution/test_model/templates/allday_superrot2/Rpallday_superrot2.txt",R_tot[::100][::-1])

def rot_int_cmj(w, s, vsini, eps=0, nr=10, ntheta=100, dif = 0.0):        #Rotational Broading Integration from https://arxiv.org/pdf/2305.09693.pdf
    '''
    A routine to quickly rotationally broaden a spectrum in linear time.
    INPUTS:
    s - input spectrum
    w - wavelength scale of the input spectrum
    
    vsini (km/s) - projected rotational velocity
    
    OUTPUT:
    ns - a rotationally broadened spectrum on the wavelength scale w
    OPTIONAL INPUTS:
    eps (default = 0) - the coefficient of the limb darkening law
    
    nr (default = 10) - the number of radial bins on the projected disk
    
    ntheta (default = 100) - the number of azimuthal bins in the largest radial annulus
                            note: the number of bins at each r is int(r*ntheta) where r < 1
    
    dif (default = 0) - the differential rotation coefficient, applied according to the law
    Omeg(th)/Omeg(eq) = (1 - dif/2 - (dif/2) cos(2 th)). Dif = .675 nicely reproduces the law 
    proposed by Smith, 1994, A&A, Vol. 287, p. 523-534, to unify WTTS and CTTS. Dif = .23 is 
    similar to observed solar differential rotation. Note: the th in the above expression is 
    the stellar co-latitude, not the same as the integration variable used below. This is a 
    disk integration routine.
    '''

    ns = np.copy(s)*0.0
    tarea = 0.0
    dr = 1./nr
    for j in range(0, nr):
        r = dr/2.0 + j*dr
        area = ((r + dr/2.0)**2 - (r - dr/2.0)**2)/int(ntheta*r) * (1.0 - eps + eps*np.cos(np.arcsin(r)))
        for k in range(0,int(ntheta*r)):
            th = np.pi/int(ntheta*r) + k * 2.0*np.pi/int(ntheta*r)
            if dif != 0:
                vl = vsini * r * np.sin(th) * (1.0 - dif/2.0 - dif/2.0*np.cos(2.0*np.arccos(r*np.cos(th))))
                ns += area * np.interp(w + w*vl/2.9979e5, w, s)
                tarea += area
            else:
                vl = r * vsini * np.sin(th)
                ns += area * np.interp(w + w*vl/2.9979e5, w, s)
                tarea += area
          
    return ns/tarea





