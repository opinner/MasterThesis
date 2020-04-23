############################################################
#this program loads all profiles from a cruise, removes outliers
#and makes a scatter plot of dissipation in the lowermost 10m 
#and the slope of the basin. 

#TODO plot mean dissipation per transect
##############################################################
import numpy as np
import scipy.io as sio
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

import geopy.distance as geo
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
#warnings.filterwarnings('ignore')
    


datafile_path = "/home/ole/share-windows/emb217_mss_data/TR1-4.mat"

transect_name = DATAFILENAME[:-4]
    
print("\n",transect_name)

#########################################
#Load data###############################
#########################################
data = np.load(datafile_path)

number_of_profiles = data["number_of_profiles"] #
lat = data["lat"] #Latitude of the profiles
lon = data["lon"] #Longitude of the profiles
distance = data["distance"] #distance from the starting profile (monotonically increasing)

interp_pressure = data["interp_pressure"]
oxygen_grid = data["oxygen_grid"]
salinity_grid = data["salinity_grid"]
consv_temperature_grid = data["consv_temperature_grid"]
density_grid = data["density_grid"]

eps_pressure = data["eps_pressure"]
eps_grid = data["eps_grid"]
corrected_eps_grid = data["corrected_eps_grid"]
corrected_eps_wiki_grid = data["corrected_eps_wiki_grid"]
eps_consv_temperature_grid = data["eps_consv_temperature_grid"]
eps_oxygen_grid = data["eps_oxygen_grid"] 

eps_N_squared_grid = data["eps_N_squared_grid"]
eps_density_grid = data["eps_density_grid"]
#eps_viscosity_grid = data["eps_viscosity_grid"]
eps_Reynolds_bouyancy_grid = data["eps_Reynolds_bouyancy_grid"]
corrected_eps_Reynolds_bouyancy_grid = data["corrected_eps_Reynolds_bouyancy_grid"]
eps_wiki_Reynolds_bouyancy_grid = data["eps_wiki_Reynolds_bouyancy_grid"]
corrected_eps_wiki_Reynolds_bouyancy_grid = data["corrected_eps_wiki_Reynolds_bouyancy_grid"]

"""
number_of_profiles              number of profiles/casts in the transect
lat                             latitude in degrees (as a float) of the casts
lon                             longitude in degrees (as a float) of the casts
distance                        distance in km from the starting point of the transect

interp_pressure                 equidistant 1D pressure array between the highest and the lowest measured pressure value
oxygen_grid                     oxygen concentration in in micromol per litre as a grid (number_of_profiles x len(interp_pressure))
salinity_grid                   salinity in g/kg as a grid (number_of_profiles x len(interp_pressure)) 
consv_temperature_grid          conservative temperature in degrees Celsius as a grid (number_of_profiles x len(interp_pressure))
density_grid                    density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))

eps_pressure                    pressure values to the dissipation rate values (the pressure distance between points is bigger than in interp_pressure) 
eps_grid                        measured dissipation rate values (number_of_profiles x len(eps_pressure))
eps_consv_temperature_grid      conservative temperature as a grid (number_of_profiles x len(eps_pressure))
eps_oxygen_grid                 oxygen concentration in micromol per litre as a grid (number_of_profiles x len(eps_pressure))
eps_N_squared_grid              N^2, the Brunt-Vaisala frequency in 1/s^2 as a grid (number_of_profiles x len(eps_pressure))
eps_density_grid                density in kg/m^3 as a grid (number_of_profiles x len(eps_pressure))

eps_viscosity_grid
eps_Reynolds_bouyancy_grid
corrected_eps_Reynolds_bouyancy_grid 
eps_wiki_Reynolds_bouyancy_grid
corrected_eps_wiki_Reynolds_bouyancy_grid 


"""

print("Number of profiles:",number_of_profiles)

print(min(eps_pressure),max(eps_pressure),len(eps_pressure))

#########################################
#Get bathymetry##########################
#########################################
#calculate the idices of the bottom and some meters above that
results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid,height_above_ground = height_above_ground)
"""
bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
"""
bathymetrie,list_of_bathymetrie_indices = results[0]
#BBL,list_of_BBL_indices = results[1] #not needed here
BBL_range,list_of_BBL_range_indices = results[2]

eps_N_grid = np.sqrt(eps_N_squared_grid)
#ozmidov scale
ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))

#conversion from pressure coordinates to depth
eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
bathymetrie_in_m = gsw.z_from_p(bathymetrie,np.mean(lat))

eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))

distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetrie_in_m,(-1,1))
boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)
    
    
#########################################
#Calculate fluxes########################
#########################################

Gamma_Osborn_eps_grid = thesis.Osborn(eps_Reynolds_bouyancy_grid)
turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)
#remove negative diffusivity 
turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid<0] = np.nan
oxygen_flux_osborn_grid = turbulent_diffusivity_Osborn_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
#convert from m*micromol/(kg*s) to mmol/(m^2*d)
oxygen_flux_osborn_grid = oxygen_flux_osborn_grid*86400*(1000/eps_density_grid)        

Gamma_BB_eps_grid = thesis.BB(eps_Reynolds_bouyancy_grid)
turbulent_diffusivity_BB_grid = Gamma_BB_eps_grid * eps_grid / (eps_N_squared_grid)
#remove negative diffusivity    
turbulent_diffusivity_BB_grid[turbulent_diffusivity_BB_grid<0] = np.nan
oxygen_flux_BB_grid = turbulent_diffusivity_BB_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
#convert from m*micromol/(kg*s) to mmol/(m^2*d)
oxygen_flux_BB_grid = oxygen_flux_BB_grid*86400*(1000/eps_density_grid)

Gamma_Skif_eps_grid = thesis.Skif(eps_Reynolds_bouyancy_grid)
turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
#remove negative diffusivity
turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid<0] = np.nan
oxygen_flux_Skif_grid = turbulent_diffusivity_Skif_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
#convert from m*micromol/(kg*s) to mmol/(m^2*d)
oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400*(1000/eps_density_grid)


#########################################
#PLotting###############################
#########################################

f1,axarray1 = plt.subplots(nrows = 1, ncols = 4, sharey = True)


axarray1[0].plot(oxygen)
axarray1[1].plot(N_squared)
axarray1[1].plot(Reb)
axarray1[1].plot(Gamma)

