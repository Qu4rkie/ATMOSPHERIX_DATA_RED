import matplotlib.pyplot as plt
import numpy as np

from scipy import stats
import pandas as pd


def calc_likelihood_HR(corr,like_type):

    try :
        (len(corr["data"]) == len(corr["model"]))
        (len(corr["data"]) == len(corr["std"]))
    except:
        raise NameError("data and model not equal")
        exit()

    if like_type=="Brogi":
        like = np.zeros(len(corr["data"]))
        for i  in range(len(corr["data"])):
            N = len(corr["data"][i])
            sf  = np.var(corr["data"][i])
            sg  = np.var(corr["model"][i])

            dat = np.array(corr["data"][i])-np.mean(np.array(corr["data"][i]))
            mod = np.array(corr["model"][i])-np.mean(np.array(corr["model"][i]))

            Rs = 1./N*np.sum(dat*mod)
            like[i]= -N/2.*np.log(sf+sg-2.*Rs)
            
        return np.sum(like)

    
           
    elif like_type=="Gibson":
        like = np.zeros(len(corr["data"]))
        for i  in range(len(corr["data"])):
            #print("Data :" +str(len(corr["data"][i])) +
            #" / Model :" + str(len(corr["model"][i])) +
            #" / Std :" +str(len(corr["std"][i])), flush=True)
            N = len(corr["data"][i])
            dat = np.array(corr["data"][i])-np.mean(np.array(corr["data"][i]))
            mod = np.array(corr["model"][i])-np.mean(np.array(corr["model"][i]))
            tolog = np.sum((dat-mod)**2/np.array(corr["std"][i])**2)/N
            
            if tolog<=0.0:
                print(tolog)
            like[i]= -N/2.*(np.log(tolog))
            
        return np.sum(like)
            
    elif like_type == "Gibson_global":
        like = np.zeros(len(corr["data"]))
        Ntot = 0
        for i  in range(len(corr["data"])):
            N = len(corr["data"][i])
            Ntot += N
            dat = np.array(corr["data"][i])-np.mean(np.array(corr["data"][i]))
            mod = np.array(corr["model"][i])-np.mean(np.array(corr["model"][i]))
            like[i] = np.sum((dat-mod)**2/np.array(corr["std"][i])**2)
            if like[i]<=0.0:
                print(i)
                exit()
            
        liketot= -Ntot/2.*(np.log(np.sum(like)/Ntot))
        return liketot      
                



    elif like_type == "Gibson_transit":
        liketot = 0
        k = 0
        for j in range(len(corr["number"])):
            Ntot = 0
            like = np.zeros(corr["number"][j])
            for i  in range(corr["number"][j]):       
                N = len(corr["data"][k])
                Ntot += N
                dat = np.array(corr["data"][k])-np.mean(np.array(corr["data"][k]))
                mod = np.array(corr["model"][k])-np.mean(np.array(corr["model"][k]))
                like[i] = np.sum((dat-mod)**2/np.array(corr["std"][k])**2)
                if like[i]<=0.0:
                    print(i)
                    exit()
                k+=1
            liketot += -Ntot/2.*(np.log(np.sum(like)/Ntot))
                
            
        
        return liketot        
        
def TP_prior_smooth(sigma_smooth,N,TProfile,logp_bot,logp_top):
    """
    Function to penalize the second derivative of the free TP profile, taken from Pelletier+2021.
    Limitation: knots have to be uniformly distributed in log10 pressure

    Parameters
    ----------
    sigma_smooth : float
        Smoothing parameter in units of Kelvin per pressure dex squared.  Low for rigid TP (strong penalty).  High for flexible TP (weak penalty).
    N : int
        Number of TP knots fitted.
    TProfile : array
        Numpy array of temperature points fitted at each knot. Ordered from top of atmosphere to bottom
    logp_bot : float
        log_10 pressure in bars at the bottom most (highest pressure) knot fitted.
    logp_top : float
        log_10 pressure in bars at the upper most (lowest pressure) knot fitted.

    Returns
    -------
    float
        log prior penalty

    """
    
    deltalogp = np.linspace(logp_top,logp_bot,num=N,retstep=True)[1]
    return (-1.0/(2.0*sigma_smooth**2)) * (1/(logp_bot-logp_top)) * np.sum(((TProfile[2:] - 2*TProfile[1:-1] + TProfile[:-2])**2)/(deltalogp**3)) - 0.5 * np.log(2*np.pi*sigma_smooth**2)
