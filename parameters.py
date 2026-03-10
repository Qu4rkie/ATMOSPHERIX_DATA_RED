
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal



type_obs = "emission"
#type_obs = "transmission"

READ_DATA = False   #do you want to read some t.fits files ?
INJ_PLANET = False  #do you want to inject a planet ?
REDUCE_DATA = True #do you want to reduce one or several pkl file that has been read beforehand ?
CORREL_DATA = True  #do you want to perform correlation from reduced pkl files ? 

dir_data = "/group/exoplanetes/Observations/"
dir_global = "/user/home/yarivv/ATMOSPHERIX_DATA_RED/"

### Directory to save figures if plot = True
dir_figures = dir_global+"Figures/"

num_obs = 4 #Number of observing nights that will be treated independently
#before being added up in the correlation

#Choose instrument
instrument="SPIROU"  #allowed values "SPIROU", "NIRPS", "HARPS", "IGRINS" (soon tm)

#Choose planet and import parameters to read data
planet_name = "TOI-2109b"
planet_file = dir_global + "Planets/" + planet_name + "_params.py"

with open(planet_file) as file:
    exec(file.read())

###########################################################################
###########################################################################
################### PARAMETERS TO READ DATA
###########################################################################
###########################################################################


### Directory where all the "t.fits" files are stores 
dir_data = [
    #dir_data+"Data/NIRPS/APERO/TOI-2109/2023-04-30/",
    #dir_data + "Data/Data_Challenge/WASP76b_H2O_CO_in_WASP107/data_with_injected_synthetic",
    dir_data+"Data/NIRPS/APERO/TOI-2109/2023-05-14/",
    #dir_data+"Data/NIRPS/APERO/TOI-2109/2025-04-29/",
    dir_data+"Data/NIRPS/APERO/TOI-2109/2025-05-06/",
    dir_data+"Data/NIRPS/APERO/TOI-2109/2025-05-08/",
    dir_data+"Data/NIRPS/APERO/TOI-2109/2025-05-10/",
    dir_data+"Data/NIRPS/APERO/TOI-2109/2025-06-06/",
    dir_data+"Data/NIRPS/APERO/TOI-2109/2025-06-29/",

    #dir_data+"Data/SPIROU/TOI-2109/2025-06-04",
    #dir_data+"Data/SPIROU/TOI-2109/2025-06-07",   #Day
    #dir_data+"Data/SPIROU/TOI-2109/2025-06-08",
    #dir_data+"Data/SPIROU/TOI-2109/2025-06-09",   #Day
    #dir_data+"Data/SPIROU/TOI-2109/2025-06-10",
    #dir_data+"Data/SPIROU/TOI-2109/2025-07-09",
    #dir_data+"Data/SPIROU/TOI-2109/2025-07-10",   #Day
    #dir_data+"Data/SPIROU/TOI-2109/2025-07-13",
    #dir_data+"Data/SPIROU/TOI-2109/2025-07-14",   #Day

    #dir_data+"Data/HARPS/TOI-2109/2023-04-30",
    #dir_data+"Data/HARPS/TOI-2109/2023-05-14",
    #dir_data+"Data/HARPS/TOI-2109/2025-04-28",
    #dir_data+"Data/HARPS/TOI-2109/2025-05-05",
    #dir_data+"Data/HARPS/TOI-2109/2025-05-07",
    #dir_data+"Data/HARPS/TOI-2109/2025-05-09",
    #dir_data+"Data/HARPS/TOI-2109/2025-06-05",
    #dir_data+"Data/HARPS/TOI-2109/2025-06-28",
        ]

### Name of the pickle file to store the info in 
dir_save_read = dir_global+"pickle/read/"
read_name_fin = [
    #"TOI2109_30-04-23_NIRPS_read.pkl",
    #"WASP76_in_WASP107_read.pkl",
    "TOI2109_14-05-23_NIRPS_read.pkl",
    #"TOI2109_29-04-25_NIRPS_read.pkl",
    "TOI2109_06-05-25_NIRPS_read.pkl",
    "TOI2109_08-05-25_NIRPS_read.pkl",
    "TOI2109_10-05-25_NIRPS_read.pkl",
    "TOI2109_06-06-25_NIRPS_read.pkl",
    "TOI2109_29-06-25_NIRPS_read.pkl",

    #"TOI2109_04-06-25_SPIROU_read.pkl",
    #"TOI2109_07-06-25_SPIROU_read.pkl",
    #"TOI2109_08-06-25_SPIROU_read.pkl",
    #"TOI2109_09-06-25_SPIROU_read.pkl",
    #"TOI2109_10-06-25_SPIROU_read.pkl",
    #"TOI2109_09-07-25_SPIROU_read.pkl",
    #"TOI2109_10-07-25_SPIROU_read.pkl",
    #"TOI2109_13-07-25_SPIROU_read.pkl",
    #"TOI2109_14-07-25_SPIROU_read.pkl",

    #"TOI2109_14-05-23_HARPS_read.pkl",
    #"TOI2109_29-04-25_HARPS_read.pkl",
    #"TOI2109_06-05-25_HARPS_read.pkl",
    #"TOI2109_08-05-25_HARPS_read.pkl",
    #"TOI2109_10-05-25_HARPS_read.pkl",
    #"TOI2109_06-06-25_HARPS_read.pkl",
    #"TOI2109_29-06-25_HARPS_read.pkl",
    ]


if instrument == "SPIROU":
    ### List of SPIRou absolute orders -- Reddest: 31; Bluest: 79
    orders   =  np.arange(31,80)[::-1].tolist() 
elif instrument == "NIRPS":
    ### List of NIRPS absolute orders -- Reddest: 0; Bluest: 74
    orders   =  np.arange(0,75)[::-1].tolist()
elif instrument == "IGRINS":
    ### List of IGRINS orders -- Reddest: 53; Bluest: 0
    orders   =  np.arange(0,54).tolist()
elif instrument == "HARPS":
    ### List of HARPS orders -- Reddest: 71; Bluest: 0
    orders   =  np.arange(0,71).tolist()
nord = len(orders)

### Plots
plot_read     = True     # If True, plot transit info
figure_name_transit = []
for pkl_name in read_name_fin:
    figure_name_transit.append(dir_figures+"Transit/transit_"+pkl_name[:-4]+".png")


###########################################################################
###########################################################################
################### PARAMETERS TO INJECT PLANET 
###########################################################################
###########################################################################


planet_wavelength_nm_file = "/home/florian/Bureau/Atmosphere_SPIRou/Pipeline_v2/uptodate/ATMOSPHERIX_DATA_RED/Models/Results/lambdastest-GL15A_onlyH2O.txt"

# Radius (in meters) for transit. CAREFUL : it is not the transit depth, it is indeed radius
planet_radius_m_file = "/home/florian/Bureau/Atmosphere_SPIRou/Pipeline_v2/uptodate/ATMOSPHERIX_DATA_RED/Models/Results/RpGL15A_HD189_onlyH2O-VMR3-T900.txt"

#planetary flux in J/m-2/s-1
planet_flux_SI_file = "/home/florian/Bureau/Atmosphere_SPIRou/Pipeline_v2/uptodate/ATMOSPHERIX_DATA_RED/Models/Results/fluxtest-GL15A_onlyH2O.txt"
K_inj = 120.
V_inj = 10.
amp_inj = 1.





###########################################################################
###########################################################################
################### PARAMETERS TO REDUCE DATA
###########################################################################
###########################################################################

dir_reduce_in = dir_global+"pickle/read/"
dir_reduce_out = dir_global+"pickle/reduced/"
reduce_name_in = [    
    #"TOI2109_30-04-23_NIRPS_read.pkl", 
    #"WASP76_in_WASP107_read.pkl",                 
    #"TOI2109_14-05-23_NIRPS_read.pkl",
    #"TOI2109_29-04-25_NIRPS_read.pkl",
    #"TOI2109_06-05-25_NIRPS_read.pkl",
    #"TOI2109_08-05-25_NIRPS_read.pkl",
    #"TOI2109_10-05-25_NIRPS_read.pkl",
    #"TOI2109_06-06-25_NIRPS_read.pkl",
    #"TOI2109_29-06-25_NIRPS_read.pkl",

    #"TOI2109_04-06-25_SPIROU_read.pkl",
    "TOI2109_07-06-25_SPIROU_read.pkl",
    #"TOI2109_08-06-25_SPIROU_read.pkl",
    "TOI2109_09-06-25_SPIROU_read.pkl",
    #"TOI2109_10-06-25_SPIROU_read.pkl",
    #"TOI2109_09-07-25_SPIROU_read.pkl",
    "TOI2109_10-07-25_SPIROU_read.pkl",
    #"TOI2109_13-07-25_SPIROU_read.pkl",
    "TOI2109_14-07-25_SPIROU_read.pkl",

    #"TOI2109_14-05-23_HARPS_read.pkl",
    #"TOI2109_29-04-25_HARPS_read.pkl",
    #"TOI2109_06-05-25_HARPS_read.pkl",
    #"TOI2109_08-05-25_HARPS_read.pkl",
    #"TOI2109_10-05-25_HARPS_read.pkl",
    #"TOI2109_06-06-25_HARPS_read.pkl",
    #"TOI2109_29-06-25_HARPS_read.pkl",
    ]
#output file names
reduce_name_out  = [    
    #"TOI2109_30-04-23_NIRPS_reduced.pkl", 
    #"WASP76_in_WASP107_reduced.pkl",                 
    #"TOI2109_14-05-23_NIRPS_5pc_reduced.pkl",
    #"TOI2109_29-04-25_NIRPS_5pc_reduced.pkl",
    #"TOI2109_06-05-25_NIRPS_5pc_reduced.pkl",
    #"TOI2109_08-05-25_NIRPS_5pc_reduced.pkl",
    #"TOI2109_10-05-25_NIRPS_5pc_reduced.pkl",
    #"TOI2109_06-06-25_NIRPS_5pc_reduced.pkl",
    #"TOI2109_29-06-25_NIRPS_5pc_reduced.pkl",

    #"TOI2109_04-06-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_07-06-25_SPIROU_5pc_nomask_APERO.pkl",
    #"TOI2109_08-06-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_09-06-25_SPIROU_5pc_nomask_APERO.pkl",
    #"TOI2109_10-06-25_SPIROU_5pc_reduced.pkl",
    #"TOI2109_09-07-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_10-07-25_SPIROU_5pc_nomask_APERO.pkl",
    #"TOI2109_13-07-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_14-07-25_SPIROU_5pc_nomask_APERO.pkl",

    #"TOI2109_14-05-23_HARPS_5pc_reduced.pkl",
    #"TOI2109_29-04-25_HARPS_reduced.pkl",
    #"TOI2109_06-05-25_HARPS_reduced.pkl",
    #"TOI2109_08-05-25_HARPS_reduced.pkl",
    #"TOI2109_10-05-25_HARPS_reduced.pkl",
    #"TOI2109_06-06-25_HARPS_reduced.pkl",
    #"TOI2109_29-06-25_HARPS_reduced.pkl",
    ]
#information file
reduce_info_file = dir_figures+"info.dat"

#list of orders we don't want to reduce for whatever reason

orders_rem =np.append(31,np.arange(34,80)).tolist()

### Correction of stellar contamination
### Only used if synthetic spectrum available
corr_star  = False
WC_name    = ""            ### Input wavelength for synthetic stellar spectra
IC_name    = ""            ### Input flux for synthetic stellar spectra

#Detect and remove spectra with outlyingly low median intensity? - experimental
select_frames = False
sigma_frames  = 2      #sigma-value to remove low-intensity spectra

### Additional Boucher correction. If dep_min >=1, not included. 
dep_min  = 0.5  # remove all data when telluric relative absorption > 1 - dep_min (default - 0.7)
thres_up = 0.1      # Remove the line until reaching 1-thres_up (default - 0.1)
Npt_lim  = 600      # If the order contains less than Npt_lim points, it is discarded from the analysis

#Do we delete a master spectrum ?  
delete_master = True

#Do we want this master spectrum to be from the APERO Solution ? 
#Careful, if  we do, we cannot have a prenormalisation hence first_norm_type must
#be set to "none"
master_from_file = True
#master_W_ref_file = "~/ATMOSPHERIX_DATA_RED/Templates/C7A164F31A_pp_e2dsff_A_wavesol_ref_A.fits"
master_W_ref_file = "~/ATMOSPHERIX_DATA_RED/Templates/REF_WAVE_2400416c_AB.fits"
master_I_ref_file = "~/ATMOSPHERIX_DATA_RED/Templates/Template_TOI2109_tellu_obj_AB.fits"

#If not from apero, we can then manually decide the where is the transit in the phase direction,
#and exclude it for the calculation of the mean stellar spectrum.
#If set_window = False, the transit window defines n_ini and n_end
set_window = False
n_ini_fix,n_end_fix = 0,-1   ### Get transits start and end indices


### Interpolation parameters
if instrument=="SPIROU":
    pixel = np.linspace(-1.14,1.14,15) ### Sampling a SPIRou pixel in velocity space -- Width ~ 2.28 km/s
    sig_g    = 2.28                    ### STD of one SPIRou px in km/s
elif instrument=="NIRPS":
    pixel = np.linspace(-0.48,0.48,15) 
    sig_g = 0.96                     ### STD of one NIRPS px in km/s -- Width ~0.96 km/s for 3.9 px per FWHM
elif instrument=="HARPS":
    pixel = np.linspace(-0.41,0.41,15)
    sig_g = 0.82                       ### STD of one HARPS px in km/s -- Width ~0.82 km/s for 3.3 px per FWHM at R~110K
elif instrument=="IGRINS":
    pixel = np.linspace(-1.0,1.0,15)
    sig_g = 2.0                       ### STD of one IGRINS px in km/s -- Width ~2.0 km/s


N_bor    = 15                           ### Nb of pts removed at each extremity (twice)





### Normalisation parameters
first_norm_type = "percentile"
if master_from_file:
    first_norm_type = "none"
second_norm_type = "old"

N_med    = 150                          ### Nb of points used in the median filter for the inteprolation
sig_out  = 5.0                          ### Threshold for outliers identification during normalisation process 
N_adj = 2 ### Number of adjacent pixel removed with outliers
deg_px   = 2                            ### Degree of the polynomial fit to the distribution of pixel STDs

### Parameters for detrending with airmass
det_airmass = False
deg_airmass = 2

### Parameters PCA. Auto-tune automatically decides the number of component 
#to remove by comparing with white noise map.
mode_pca    = "pca"                     ### "pca" or "autoencoder"
wpca = False   #Use weighted pca
auto_tune   = False                 ### Automatic tuning of number of components
factor_pca = 1 #Factor in the auto tune: every PC above factor*white_noise_mean_eigenvalue is suppressed
min_pca  = 1 #minimim number of removed components
mode_norm_pca = "per_obs" #how to remove mean and std in the data before PCA. Four possibilities:
                         # "none" : data untouched.
                         # "global" : suppression of mean and division by the std of the whole data set 
                         # 'per_pix': same as global but column by colum (per pixel)
                         # 'per_obs': same as global but line by line (per observation)

 ### Nb of removed components if auto tune is false
npca        = np.array(5*np.ones(nord),dtype=int)     


### Plot info
plot_red    = True
numb        = 33 #e.g. 46 for SPIROU, 15 for NIRPS - 65 for H Band
    

#If you want to remove some orders, put them here
orders_rem     = [[55], [55], [55,56], [55,56],[],[],[],[],[]]
#[25,24,23,22,21,20] NIRPS

### Size of the estimation of the std of the order for final metrics
N_px          = 200



###########################################################################
###########################################################################
################### PARAMETERS FOR CORRELATION
###########################################################################
############################################################################

parallel = True
#This is just for the integration over a pixel
if instrument=="SPIROU":
    pixel_correl = np.linspace(-1.14,1.14,15)
elif instrument=="NIRPS":
    pixel_correl = np.linspace(-0.48,0.48,15)
elif instrument=="HARPS":
    pixel_correl = np.linspace(-0.41,0.41,15)
elif instrument=="IGRINS":
    pixel_correl = np.linspace(-1.0,1.0,15)
weights= np.ones(15)

#Kp intervals
Kpmin = 150.0
Kpmax = 350.0
Nkp = 101
Kp_array = np.linspace(Kpmin,Kpmax,Nkp)

#Vsys intervals
Vmin = -100.
Vmax= 100
Nv = 201
Vsys_array = np.linspace(Vmin,Vmax,Nv)

#Number of pkl observations files and their names

dir_correl_in = dir_global+"pickle/reduced/"


correl_name_in = [    
    #"TOI2109_30-04-23_reduced.pkl", 
    #"WASP76_in_WASP107_reduced.pkl",                 
    #"TOI2109_14-05-23_NIRPS_5pc_reduced.pkl",
    #"TOI2109_29-04-25_NIRPS_5pc_reduced.pkl",
    #"TOI2109_06-05-25_NIRPS_5pc_apero_mask_reduced.pkl",
    #"TOI2109_08-05-25_NIRPS_5pc_apero_mask_reduced.pkl",
    #"TOI2109_10-05-25_NIRPS_5pc_apero_mask_reduced.pkl",
    #"TOI2109_06-06-25_NIRPS_5pc_apero_mask_reduced.pkl",
    #"TOI2109_29-06-25_NIRPS_5pc_reduced.pkl",

    #"TOI2109_04-06-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_07-06-25_SPIROU_5pc_nomask_APERO.pkl",
    #"TOI2109_08-06-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_09-06-25_SPIROU_5pc_nomask_APERO.pkl",
    #"TOI2109_10-06-25_SPIROU_5pc_reduced.pkl",
    #"TOI2109_09-07-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_10-07-25_SPIROU_5pc_nomask_APERO.pkl",
    #"TOI2109_13-07-25_SPIROU_5pc_reduced.pkl",
    "TOI2109_14-07-25_SPIROU_5pc_nomask_APERO.pkl",

    #"TOI2109_14-05-23_HARPS_5pc_reduced.pkl",
    #"TOI2109_29-04-25_HARPS_reduced.pkl",
    #"TOI2109_06-05-25_HARPS_reduced.pkl",
    #"TOI2109_08-05-25_HARPS_reduced.pkl",
    #"TOI2109_10-05-25_HARPS_reduced.pkl",
    #"TOI2109_06-06-25_HARPS_reduced.pkl",
    #"TOI2109_29-06-25_HARPS_reduced.pkl",
    ]



#Do we save the correlation file ? If yes, put as much files as there are observations
save_ccf = True
dir_correl_out = dir_global+"pickle/correlated/"


correl_name_out = [    
    #"TOI2109_30-04-23__Fe_correlated.pkl",    #.   (bad SNR) - Failed Tcorr GC          
    #"TOI2109_14-05-23_NIRPS_5pc_Fe_scarlet_broad_correlated.pkl", #pre-eclipse (good SNR)
    #"TOI2109_29-04-25_NIRPS_5pc_Fe_broad_correlated.pkl", #nightside (best SNR)
    #"TOI2109_06-05-25_NIRPS_5pc_OH_apero_mask_correlated.pkl", #eclipse (best SNR: 40-48)
    #"TOI2109_08-05-25_NIRPS_5pc_OH_apero_mask_correlated.pkl", #eclipse (best SNR post eclipse: ~45) pre-eclipse SNR: 25-50
    #"TOI2109_10-05-25_NIRPS_5pc_OH_apero_mask_correlated.pkl", #eclipse (good SNR)
    #"TOI2109_06-06-25_NIRPS_5pc_OH_apero_mask_correlated.pkl", #eclipse (good-ish SNR: 20-40)
    #"TOI2109_29-06-25_NIRPS_5pc_Fe_scarlet_broad_correlated.pkl",  #post-eclipse (good-ish SNR) - Failed Tcorr QC

    #"TOI2109_04-06-25_SPIROU_5pc_OH_correlated.pkl",  # ~quad pre-T, good SNR~85 except 1 frame, but missing orders post reduction
    "TOI2109_07-06-25_SPIROU_5pc_APERO_CO_correlated.pkl",  #day pre-E, best SNR: 90-95
    #"TOI2109_08-06-25_SPIROU_5pc_Fe_broad_correlated.pkl",  #night pre-T, best SNR
    "TOI2109_09-06-25_SPIROU_5pc_APERO_CO_correlated.pkl",  #day pre-E, good SNR
    #"TOI2109_10-06-25_SPIROU_5pc_Fe_broad_correlated.pkl",  #night pre-T, good SNR
    #"TOI2109_09-07-25_SPIROU_5pc_Fe_broad_correlated.pkl",  #quad/night pre-T, good SNR - weird reduction ?
    "TOI2109_10-07-25_SPIROU_5pc_APERO_CO_correlated.pkl",  #day pre-E, good SNR
    #"TOI2109_13-07-25_SPIROU_5pc_Fe_broad_correlated.pkl",  #quad/night pre-T, best SNR
    "TOI2109_14-07-25_SPIROU_5pc_APERO_CO_correlated.pkl",  #day pre-E, good SNR

    #"TOI2109_14-05-23_HARPS_5pc_reduced.pkl",
    #"TOI2109_29-04-25_HARPS_Fe_correlated.pkl",
    #"TOI2109_06-05-25_HARPS_maroonx_correlated.pkl",
    #"TOI2109_08-05-25_HARPS_maroonx_correlated.pkl",
    #"TOI2109_10-05-25_HARPS_maroonx_correlated.pkl",
    #"TOI2109_06-06-25_HARPS_maroonx_correlated.pkl",
    #"TOI2109_29-06-25_HARPS_maroonx_correlated.pkl",
    ]

dir_correl_mod = dir_global+"Templates/TOI-2109b_CO_only_SPIROU/"
#dir_correl_mod = dir_global+"Templates/WASP33_13CO_SPIROU_broad/"

#DO we select orders or take them all ? If True, provide your order selection
# for each observation. If an order does not exist in the pkl file, it will 
# obivously not be used but will not trigger an error.
select_ord = True
list_ord_correl = np.arange(31,34) #not used if select_ord = False



#If false, the calculation is performed over the whole dataset. If 
#True, we only select observation that have a transit window > min_window
select_phase = True
min_window = 0.1

#Interpolation factor for the speed array. If you d'ont know what that means, choose something between 1 and 10
int_speed = 5

#Number of pixels to discard at the borders. 
nbor_correl = 10

#Do we include the projector from Gibson 2022 ?
use_proj = True
#If we just removed the mean and std of the whole map, we can use a fast verion of the projector
#Else, it will be even longer
proj_fast = True
mode_norm_pca_correl = "none" #if proj_fast is not used, we can choose
                          #how to remove mean and std in the data before PCA. Four possibilities:
                          # "none" : data untouched.
                          # "global" : suppression of mean and division by the std of the whole data set 
                          # 'per_pix': same as global but column by colum (per pixel)
                          # 'per_obs': same as global but line by line (per observation)

#Do we select only certain orders for the plot ? 
#if yes, lili is the list oforders to select
select_plot = False
list_ord_plot_correl = np.array([48,47,46,34,33,32])

#In order to calculate the std of the map,we need to exclude 
#a zone of the Kp-Vsys map around the planet. These are the limits 
#of this rectangular zone.
Kp_min_std = 290
Kp_max_std = 230
Vsys_min_std = -25
Vsys_max_std = 25

#number of levels in the contour plot
nlevels = 15

#Do we plot the correlation map at each obs ?
plot_ccf_indiv = True
#Do we save the individual CCF maps
save_ccf_indiv = True
#Do we plot the global correlation map ? 
plot_ccf_tot = True
#Do we save the global CCF map?
save_ccf_tot = True

save_path_indiv = [
    #dir_figures + "CCF/TOI2109_30-04-23_Fe.png",
    #dir_figures + "CCF/WASP76_in_WASP107.png",
    #dir_figures + "CCF/TOI2109_14-05-23_NIRPS_5pc_Fe_scarlet_broad.png",
    #dir_figures + "CCF/TOI2109_29-04-25_NIRPS_5pc_Fe_broad.png",
    #dir_figures + "CCF/TOI2109_06-05-25_NIRPS_5pc_OH_apero_mask.png",
    #dir_figures + "CCF/TOI2109_08-05-25_NIRPS_5pc_OH_apero_mask.png",
    #dir_figures + "CCF/TOI2109_10-05-25_NIRPS_5pc_OH_apero_mask.png",
    #dir_figures + "CCF/TOI2109_06-06-25_NIRPS_5pc_OH_apero_mask.png",
    #dir_figures + "CCF/TOI2109_29-06-25_NIRPS_5pc_Fe_scarlet_broad.png",

    #dir_figures + "CCF/TOI2109_04-06-25_SPIROU_OH_5pc.png",
    dir_figures + "CCF/TOI2109_07-06-25_SPIROU_CO_5pc_APERO.png",
    #dir_figures + "CCF/TOI2109_08-06-25_SPIROU_Fe_broad_5pc.png",
    dir_figures + "CCF/TOI2109_09-06-25_SPIROU_CO_5pc_APERO.png",
    #dir_figures + "CCF/TOI2109_10-06-25_SPIROU_Fe_broad_5pc.png",
    #dir_figures + "CCF/TOI2109_09-07-25_SPIROU_Fe_broad_5pc.png",
    dir_figures + "CCF/TOI2109_10-07-25_SPIROU_CO_5pc_APERO.png",
    #dir_figures + "CCF/TOI2109_13-07-25_SPIROU_Fe_broad_5pc.png",
    dir_figures + "CCF/TOI2109_14-07-25_SPIROU_CO_5pc_APERO.png",

    #dir_figures + "CCF/WASP33_darveau_SPIROU_13CO.png",
    #dir_figures + "CCF/WASP33_darveau2_SPIROU_13CO.png",

    #dir_figures + "CCF/TOI2109_14-05-23_HARPS_Fe.png",
    #dir_figures + "CCF/TOI2109_29-04-25_HARPS_Fe.png",
    #dir_figures + "CCF/TOI2109_06-05-25_HARPS_maroonx.png",
    #dir_figures + "CCF/TOI2109_08-05-25_HARPS_maroonx.png",
    #dir_figures + "CCF/TOI2109_10-05-25_HARPS_maroonx.png",
    #dir_figures + "CCF/TOI2109_06-06-25_HARPS_maroonx.png",
    #dir_figures + "CCF/TOI2109_29-06-25_HARPS_maroonx.png",    
    ]
save_path_tot = dir_figures +"CCF/TOI2109_SPIROU_all_daysides_CO_5pc_APERO.png"
#save_path_tot = dir_figures +"CCF/WASP33_all_13CO_broad_dayside.png"

#Do we add white lines at the planet position ? 
white_lines = True
Kp_planet = 263.
Vsys_planet = 0.


###########################################################################
###########################################################################
################### PARAMETERS FOR PLOTS
###########################################################################
############################################################################
SMALL_SIZE = 28
MEDIUM_SIZE = 32
BIGGER_SIZE = 34
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title




