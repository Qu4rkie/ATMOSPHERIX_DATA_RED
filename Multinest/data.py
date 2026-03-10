#!/usr/bin/env python3

import numpy as np

dir_global = "/user/home/yarivv/ATMOSPHERIX_DATA_RED/"
pkl_dir = dir_global + "pickle/MCMC/"

#Choose planet and import parameters to read data
planet_name = "TOI-2109b"
planet_file = dir_global + "Planets/" + planet_name + "_params.py"

G     = 6.67e-11       # [m^3.kg-1.s-2]
with open(planet_file) as file:
        exec(file.read())

def make_data(args):
    """Make data dictionary for your planet"""

     # CAREFUL !!! It must be a list of lists (I know I have to change it)
    #. If there is only one observation, put a comma
    # after the first obs
    pkl =  [[pkl_dir+"TOI2109_06-05-25_NIRPS_5pc_reduced_MCMC.pkl"],
    [pkl_dir+"TOI2109_08-05-25_NIRPS_5pc_reduced_MCMC.pkl"],
    [pkl_dir+"TOI2109_10-05-25_NIRPS_5pc_reduced_MCMC.pkl"],
    [pkl_dir+"TOI2109_06-06-25_NIRPS_5pc_reduced_MCMC.pkl"]]
    num_transit = len(pkl)

    #Choose instrument
    instrument="NIRPS"  #allowed values "SPIROU", "NIRPS", "HARPS", "IGRINS" (soon tm)
    instrument_data = np.genfromtxt(dir_global + "Instruments/wlen_" + instrument + ".dat")

    #you can use anything in the test case, where LRS is actually commented everywhere
    LRS_file = "example.txt"

    radius_RJ = Rp / Rjup
    Rs_Rsun = Rs / Rsun
    gravity_SI = G*Mp/(Rp*1e3)**2 #m/s2
    print("Derived Gravity for " + planet_name + " : " + str(gravity_SI))
    
    if instrument == "SPIROU":
        ### List of SPIRou absolute orders -- Reddest: 31; Bluest: 79
        orderstot   =  np.arange(31,80)[::-1].tolist() 
    elif instrument == "NIRPS":
        ### List of NIRPS absolute orders -- Reddest: 0; Bluest: 74
        orderstot   =  np.arange(0,75)[::-1].tolist()
    elif instrument == "IGRINS":
        ### List of IGRINS orders -- Reddest: 53; Bluest: 0
        orderstot   =  np.arange(0,54).tolist()
    elif instrument == "HARPS":
        ### List of HARPS orders -- Reddest: 71; Bluest: 0
        orderstot   =  np.arange(0,71).tolist()

    lambdas = np.zeros((130,2))
    for ind,order in enumerate(orderstot):
        lambdas[order,0] = instrument_data[ind,1]
        lambdas[order,1] = instrument_data[ind,2]



    return dict(
        radius_RJ=radius_RJ,
	    gravity_SI = gravity_SI,
	    Rs_Rsun = Rs_Rsun,
        system_distance = sy_dist,
        vsini = Vsini,
	    lambdas = lambdas,
        orderstot=orderstot,
        num_transit=num_transit,
        pkl = pkl,
        LRS_file=LRS_file,
		    )
