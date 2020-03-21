#TODO Sort density profile
#TODO First rough thorpe scale calculation

#---------------------------------------------------------#
#---------------------------------------------------------#
#Plots one mss transect 
#plus as an example one profile from that transect
#---------------------------------------------------------#

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw 
import mss_functions as thesis

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)



#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################


datafile_path = "/home/ole/share-windows/emb217_mss_data/TR1-4.mat"

splitted_filename = datafile_path.split("/")
cruisename = splitted_filename[4][0:6]
DATAFILENAME = splitted_filename[-1]
print("cruisename",cruisename)  

#define the pictures
f1, axarr1 = plt.subplots(nrows = 1, ncols = 3, sharey = True)



    
print("Filename:",sio.whosmat(datafile_path))


results = thesis.load_clean_and_interpolate_data(datafile_path)

try:
    number_of_profiles,lat,lon,distance = results[0]
    interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid = results[1]
    eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid = results[2]
except TypeError:
    print(cruisename,DATAFILENAME[:-4],"is skipped!")
            
results = thesis.find_bottom_and_bottom_currents(number_of_profiles,interp_pressure,density_grid,oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02)
bathymetrie,list_of_bathymetrie_indices = results[0]          


"""
#calculate the viscosity (2 different formula)
eps_wiki_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_density_grid,"Wikipedia")
eps_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_density_grid)

#calculate the Reynolds bouyancy number defined on eps_pressure
eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared_grid)
eps_wiki_Reynolds_bouyancy_grid = eps_grid/(eps_wiki_viscosity_grid*eps_N_squared_grid)

#find the BBL
results = thesis.find_bottom_and_bottom_currents(number_of_profiles,interp_pressure,density_grid,oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02)
bathymetrie,list_of_bathymetrie_indices = results[0]
BBL,list_of_BBL_indices = results[1]
BBL_range,list_of_BBL_range_indices = results[2]
"""


eps_N_grid = np.sqrt(eps_N_squared_grid)
#ozmidov scale
ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))

#conversion from pressure coordinates to depth
eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
bathymetrie_in_m = gsw.z_from_p(bathymetrie,np.mean(lat))

eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))

distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetrie_in_m,(-1,1))
boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)
        
        
#append the last distance plus the last difference (for plotting all the n profiles we need a distance array of size n+1 
plotmesh_distance = np.append(distance,2*distance[-1]-distance[-2])
plotmesh_longitude = np.append(lon,2*lon[-1]-lon[-2])
        
    
#Plotting
transect_index = -10
print("Profile at Longitude",lon[transect_index])

img1_1a = axarr1[0].plot(density_grid[transect_index,:],interp_pressure[:])
img1_1b = axarr1[0].plot(sorted(density_grid[transect_index,:]),interp_pressure[:])

img1_2 = axarr1[1].plot(density_grid[transect_index,:]-sorted(density_grid[transect_index,:]),interp_pressure[:])

axarr1[2].plot(ozmidov_scale_grid[transect_index,:],eps_pressure, label = "Ozmidov scale")
axarr1[2].plot(distance_from_ground_grid[transect_index,:],eps_pressure, label = "distance from ground")
axarr1[2].set_xlim(-10,10)

f1.set_size_inches(18,10.5)

axarr1[0].invert_yaxis()
axarr1[2].legend()

f1.suptitle("Thorpe scale")

f1.tight_layout()  
plt.show()


