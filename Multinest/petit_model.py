import numpy as np
import sys 
import prepare_model as prep_mod
from petitRADTRANS.radtrans import Radtrans
from petitRADTRANS import physical_constants as cst
from petitRADTRANS.physics import temperature_profile_function_guillot_global
from petitRADTRANS.spectral_model import SpectralModel
from scipy.interpolate import CubicSpline
from astropy.io import fits
from convolve import rot_int_cmj
import contextlib

def my_t_profile(pressures, temperature_points, pressure_points, **kwargs):
    TP =  CubicSpline(np.log(pressure_points),temperature_points, bc_type="natural",extrapolate=True)
    temperatures = TP(np.log(pressures))
    mask_low, mask_high = pressures < pressure_points[0], pressures > pressure_points[-1]
    temperatures[mask_low] = temperature_points[0]
    temperatures[mask_high] = temperature_points[-1]
    return temperatures

class Model(object):
    def __init__(self,config_dict):

        self.p_minbar = config_dict["p_minbar"]
        self.p_maxbar= config_dict["p_maxbar"]
        self.n_pressure=config_dict["n_pressure"]
        self.P0=config_dict["P0_bar"]
        self.radius=config_dict["radius_RJ"]*cst.r_jup_mean
        self.Rs=config_dict["Rs_Rsun"]*6.9634e10  #cm
        self.gravity_SI=config_dict["gravity_SI"]
        self.HHe_ratio=config_dict["HHe_ratio"]
        self.emission = config_dict["emission"]
        if self.emission:
            self.system_distance = config_dict["system_distance"]*1e5  #[cm]
            self.vsini = config_dict["vsini"] #[km/s]
        self.line_species = config_dict["line_species"]

        self.profile = config_dict["temperature_profile"]
        if self.profile == "guillot":
            self.kappa_IR = config_dict["kappa_IR"]
            self.gamma = config_dict["gamma"]
            self.T_int = config_dict["T_int"]
        elif self.profile == "parametric":
            self.pressure_points = config_dict['pressure_points']
        elif self.profile != "isothermal":
            print("Temperature Profile type "+self.profile+" unknown")

        self.lambdas = config_dict["lambdas"]
        self.orderstot = config_dict["orderstot"]

        self.num_transit= config_dict["num_transit"]
        self.winds = config_dict["winds"]
        self.atmospheres= []
        self.pressures=np.logspace(self.p_minbar,self.p_maxbar,self.n_pressure)

        for i in self.orderstot: 
            with contextlib.redirect_stdout(None):
                atmosphere = Radtrans(self.pressures,
                            line_species = list(self.line_species.keys()),
                            rayleigh_species = ['H2', 'He'],
                            gas_continuum_contributors=['H2--H2', 'H2--He'],
                            wavelength_boundaries=[self.lambdas[i][0]/1000.0*0.995,self.lambdas[i][1]/1000.0*1.005],
                            line_opacity_mode='lbl',
                            line_by_line_opacity_sampling=4,)
            self.atmospheres.append(atmosphere)

        if self.emission:
            print("Loading Star Spectrum...", flush=True)
            wav_pheonix = fits.open(config_dict["wav_pheonix"])
            star_pheonix = fits.open(config_dict["star_pheonix"])

            self.wl_star = wav_pheonix[0].data / 1e8 #[cm]
            I_pheonix = star_pheonix[0].data
            I_broadened = rot_int_cmj(self.wl_star,I_pheonix,80)
            self.I_star = I_broadened * (self.Rs/self.system_distance)**2

        self.abundances = {}


#### LRS if there is any need
#        self.atmosphere_LRS = Radtrans(pressures=self.pressures,
#                                       line_species=[ 'H2O', 'CO-NatAbund', 'CO2'],
#                                       rayleigh_species=['H2', 'He'],
#                                       gas_continuum_contributors=['H2--H2', 'H2--He'],
#                                       wavelength_boundaries=[2.5, 5.5])
#




    def compute_petit(self, para_dic): # creates an atmospheric model
        if "Vrot" in para_dic:
            self.rot_speed = para_dic["Vrot"]

        if self.profile == "isothermal":
            temperatures = para_dic["T_eq"]*np.ones_like(self.pressures)
        elif self.profile == "guillot_global":
            # Guillot temperature profile 
            temperatures = temperature_profile_function_guillot_global(
            pressures=self.pressures,
            infrared_mean_opacity=self.kappa_IR,
            gamma=self.gamma,
            gravities=self.gravity_SI*100,
            intrinsic_temperature=self.T_int,
            equilibrium_temperature=para_dic["T_eq"]
        )
        elif self.profile == "parametric":
            temperature_points = np.array([para_dic["T"+str(X+1)] for X in range(len(self.pressure_points))])
            temperatures = my_t_profile(self.pressures, temperature_points, self.pressure_points)
        
        Z,MW_sum=0,0
        mass_fractions = {}
        for key in list(self.line_species.keys()):
            try:
                Z += 10**para_dic["MMR_"+key]
            except:
                print(key+" MMR not in parameters!")
            MW_sum += 10**para_dic["MMR_"+key]/self.line_species[key]
            mass_frac = {key:10**para_dic["MMR_"+key]*np.ones_like(temperatures)}
            mass_fractions.update(mass_frac)

        MMR_H2 = (1.0-Z)*(1-self.HHe_ratio)
        MMR_He = self.HHe_ratio*(1.0-Z)
        MW_sum += MMR_H2/2.0+MMR_He/4.0
        MMW = 1.0/MW_sum*np.ones_like(temperatures)

        mass_fractions.update(H2= MMR_H2*np.ones_like(temperatures), He= MMR_He*np.ones_like(temperatures))
        self.mass_fractions = mass_fractions

#### LRS if there is any need
        # self.mass_fractions_LRS = {'H2':MMR_H2* np.ones_like(temperatures),
        #                        'He': MMR_He * np.ones_like(temperatures),
        #                        'H2O': 10.0**para_dic["MMR_H2O"] * np.ones_like(temperatures),
        #                        'CO-NatAbund': 10.0**para_dic["MMR_CO"] * np.ones_like(temperatures),
        #                        'CO2': 10.0**para_dic["MMR_CO2"] * np.ones_like(temperatures)}

        wavelength_nm = []
        radius_transm = []
        star_flux = []
        #Calculate the radius order by order
        if self.emission:
            for atmo in self.atmospheres:
                wavelengths, flux, _ = atmo.calculate_flux(temperatures=temperatures,
                                                            mass_fractions=self.mass_fractions,
                                                            mean_molar_masses = MMW,
                                                            reference_gravity = self.gravity_SI*100,
                                                            frequencies_to_wavelengths=True)
                #save it in nm and Fp/F*
                flux_planet = flux * (self.radius/self.system_distance)**2
                flux_star = np.interp(wavelengths,self.wl_star, self.I_star)
                wavelength_nm.append(wavelengths*1.0e7)
                radius_transm.append(flux_planet)
                star_flux.append(flux_star)
        else:
            for atmo in self.atmospheres:
                wavelengths, transit_radii, _ = atmo.calculate_transit_radii(temperatures=temperatures,
                                                                               mass_fractions=self.mass_fractions,
                                                                               mean_molar_masses=MMW,
                                                                               reference_gravity=self.gravity_SI*100,
                                                                               planet_radius=self.radius,
                                                                               reference_pressure=self.P0)

                #save it in nm and cm
                wavelength_nm.append(wavelengths*1.0e7)  #[cm]
                radius_transm.append(transit_radii)      #[cm]

        if self.winds:
            self.superrot = 0.0
        


#### LRS if there is any need
        # wavelengths_LRS, transit_radii_LRS, _ = self.atmosphere_LRS.calculate_transit_radii(
        # temperatures=temperature,
        # mass_fractions=self.mass_fractions_LRS,
        # mean_molar_masses=MMW,
        # reference_gravity=self.gravity*100*para_dic["dg"],
        # planet_radius=self.radius,
        # reference_pressure=self.P0)



        # wavelength_LR_microns = wavelengths_LRS*1.e4
        # radius_LR = transit_radii_LRS/100.0
        # tdepth_LR = (radius_LR/self.Rs)**2*100.
        
        

        return {
            "wavelength_nm": wavelength_nm,
            "radius_transm": radius_transm,
            "star_flux": star_flux
#### LRS if there is any need
            # "wavelength_LR_microns": wavelength_LR_microns,
            # "tdepth_LR": tdepth_LR,
        }

    def reduce_model(self, model_dic): #renormalizes the atmospheric mode
        #print(self.winds)
        if self.winds:
            return prep_mod.prepare(model_dic,self.Rs,self.orderstot,winds=self.winds,rot_speed=
                              self.rot_speed,superrot=self.superrot,emission=self.emission)
        else:
            return prep_mod.prepare(model_dic,self.Rs,self.orderstot,rot_speed=
                              self.rot_speed,emission=self.emission)


    def return_reduced_model(self,para_dic):
        model_dic = self.compute_petit(para_dic)
        return self.reduce_model(model_dic)


