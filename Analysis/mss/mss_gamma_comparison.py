############################################################
#TODO
##############################################################
import numpy as np
import scipy.io as sio
import geopy.distance as geo
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
#import numpy.testing as testing
import matplotlib.pyplot as plt
#matplotlib preferences:
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('savefig', dpi=300)
import warnings
#warnings.filterwarnings('ignore')
    


datafile_path = "/home/ole/share-windows/processed_mss/emb217/TR1-4.npz"

transect_name = datafile_path[-9:-4]
    
print("\n",transect_name)

#########################################
#Load data###############################
#########################################
data = np.load(datafile_path)

print(data.files)

number_of_profiles = data["number_of_profiles"] #
lat = data["lat"] #Latitude of the profiles
lon = data["lon"] #Longitude of the profiles
distance = data["distance"] #distance from the starting profile (monotonically increasing)

"""
interp_pressure = data["interp_pressure"]
oxygen_grid = data["oxygen_grid"]
oxygen_grid = data["oxygen_grid"] = data["oxygen_grid = data["oxygen_grid"]"]
#salinity_grid = data["salinity_grid"]
#consv_temperature_grid = data["consv_temperature_grid"]
#density_grid = data["density_grid"]
"""

eps_pressure = data["eps_pressure"]
eps_grid = data["eps_grid"]
corrected_eps_grid = data["corrected_eps_grid"]
corrected_eps_wiki_grid = data["corrected_eps_wiki_grid"]
eps_consv_temperature_grid = data["eps_consv_temperature_grid"]
eps_oxygen_grid = data["eps_oxygen_grid"] 
eps_oxygen_sat_grid = data["eps_oxygen_sat_grid"]   

eps_N_squared_grid = data["eps_N_squared_grid"]
eps_density_grid = data["eps_density_grid"]
#eps_viscosity_grid = data["eps_viscosity_grid"]
eps_Reynolds_bouyancy_grid = data["eps_Reynolds_bouyancy_grid"]
corrected_eps_Reynolds_bouyancy_grid = data["corrected_eps_Reynolds_bouyancy_grid"]
eps_wiki_Reynolds_bouyancy_grid = data["eps_wiki_Reynolds_bouyancy_grid"]
corrected_eps_wiki_Reynolds_bouyancy_grid = data["corrected_eps_wiki_Reynolds_bouyancy_grid"]


print("Number of profiles:",number_of_profiles)

print(min(eps_pressure),max(eps_pressure),len(eps_pressure))

#conversion from pressure coordinates to depth
eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
 
    
#########################################
#Calculate fluxes########################
#########################################

Gamma_Osborn_eps_grid = thesis.Osborn(eps_Reynolds_bouyancy_grid)
turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)
#remove negative diffusivity 
turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid<0] = np.nan
oxygen_flux_osborn_grid = - turbulent_diffusivity_Osborn_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
#convert from m*micromol/(kg*s) to mmol/(m^2*d)
oxygen_flux_osborn_grid = oxygen_flux_osborn_grid*86400*(1000/eps_density_grid)        

Gamma_BB_eps_grid = thesis.BB(eps_Reynolds_bouyancy_grid)
turbulent_diffusivity_BB_grid = Gamma_BB_eps_grid * eps_grid / (eps_N_squared_grid)
#remove negative diffusivity    
turbulent_diffusivity_BB_grid[turbulent_diffusivity_BB_grid<0] = np.nan
oxygen_flux_BB_grid = - turbulent_diffusivity_BB_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
#convert from m*micromol/(kg*s) to mmol/(m^2*d)
oxygen_flux_BB_grid = oxygen_flux_BB_grid*86400*(1000/eps_density_grid)

Gamma_Skif_eps_grid = thesis.Skif(eps_Reynolds_bouyancy_grid)
turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
#remove negative diffusivity
turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid<0] = np.nan
oxygen_flux_Skif_grid = - turbulent_diffusivity_Skif_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
#convert from m*micromol/(kg*s) to mmol/(m^2*d)
oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400*(1000/eps_density_grid)

#compare with the imported functions
oxygen_flux_osborn_grid2 = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
np.testing.assert_equal(oxygen_flux_osborn_grid,oxygen_flux_osborn_grid2)
oxygen_flux_BB_grid2 = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
np.testing.assert_equal(oxygen_flux_BB_grid,oxygen_flux_BB_grid2)
oxygen_flux_Skif_grid2 = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
np.testing.assert_equal(oxygen_flux_Skif_grid,oxygen_flux_Skif_grid2)
        
#########################################
#Plotting###############################
#########################################

f1,axarray1 = plt.subplots(nrows = 1, ncols = 6, sharey = True)

profile_index = -5

axarray1[0].plot(eps_oxygen_sat_grid[profile_index,:],eps_pressure,)
axarray1[1].plot(eps_N_squared_grid[profile_index,:],eps_pressure)
axarray1[2].plot(np.log10(eps_grid[profile_index,:]),eps_pressure)
axarray1[3].plot(eps_Reynolds_bouyancy_grid[profile_index,:],eps_pressure)
axarray1[4].plot(Gamma_Osborn_eps_grid[profile_index,:],eps_pressure, label = "Osborn")
axarray1[4].plot(Gamma_Skif_eps_grid[profile_index,:],eps_pressure, label = "Shih et al")
axarray1[4].plot(Gamma_BB_eps_grid[profile_index,:],eps_pressure, label = "Bouffard et al")

axarray1[5].plot(oxygen_flux_osborn_grid[profile_index,:],eps_pressure, label = "Osborn")
axarray1[5].plot(oxygen_flux_Skif_grid[profile_index,:],eps_pressure, label = "Shih et al")
axarray1[5].plot(oxygen_flux_BB_grid[profile_index,:],eps_pressure, label = "Bouffard et al")

axarray1[0].set_ylim(0,125)
axarray1[0].invert_yaxis()
axarray1[0].set_xlabel("oxygen sat [%]")
axarray1[0].set_ylabel("pressure [dbar]")
axarray1[1].set_xlabel(r"$N^2$ $[1/s^2]$")
axarray1[2].set_xlabel(r"log10(eps) $[m^2 s^{-3}]$")
axarray1[3].set_xlabel(r"Re$_b$")
axarray1[3].set_xlim(-10,500)
axarray1[4].set_xlabel(r"$\Gamma$")
axarray1[5].set_xlabel(r"oxygen flux [mmol/(m$^2$*d]")
axarray1[5].set_xlim(-15,15)

axarray1[4].legend(loc = "upper center")
axarray1[5].legend(loc = "upper center")

f1.set_size_inches(18,10.5)
f1.tight_layout()
plt.show()

