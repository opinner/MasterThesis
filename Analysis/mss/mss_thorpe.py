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



def experimental_BBL(number_of_profiles,interp_pressure,density_grid,oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02):
    import numpy as np

    """
    procedural function    
    
    searches for the highermost nanvalues as the ground
    calculates the position of the bottom boundary layer    
    
    input:
       number_of_profiles               number of profiles/casts in the transect
       interp_pressure
       
       density_grid                     density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
       oxygen_grid                      oxygen saturation in percent as a grid (number_of_profiles x len(interp_pressure))
       
       height_above_ground              Default value 10m
       minimal_density_difference       Default value 0.02 kg/m^3
    
    
    return values:
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
    
    
    """    
    #search for bottom currents
    ###########################################################################################################################################################
    bathymetrie = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_bathymetrie_indices = np.zeros(number_of_profiles)

    halocline = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_halocline_indices = np.zeros(number_of_profiles)

    BBL = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_BBL_indices = np.zeros(number_of_profiles)

    BBL_range = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_BBL_range_indices = np.zeros(number_of_profiles)

    for i in range(number_of_profiles):

        #------------------search for bathymetrie values starts from below:-------------------------------
        #search is done in the fine grid

        #returns the pressure of the last nan value in a continuous row starting from high pressure (TODO:is it better to use the last index with data?)
        nan_index =  -np.argmax(np.flip(~np.isnan(density_grid[i,:]))) #at the moment the index is negative
        nan_index = density_grid[i,:].size + nan_index #now defined as positive index
        
       
        if nan_index == density_grid[i,:].size:
            if not np.isnan(density_grid[i,-1]): #if there are no NAN values towards the bottom
                nan_index = len(interp_pressure)-1
        
        
        assert(nan_index>=0)
                
        list_of_bathymetrie_indices[i] = nan_index 
        bathymetrie[i] = interp_pressure[nan_index]
        
      
        #TODO
        #------------------search for halocline values starts from above:-------------------------------
        


        #------------------search for BBL values starts from below:-------------------------------  
        
        #index of maximal distance bottom plus 15m 
        
        BBL_boundary_index = np.argmax(interp_pressure >= (bathymetrie[i]-height_above_ground))
        assert(interp_pressure[BBL_boundary_index]<bathymetrie[i]) #tests if the point 15m above the ground is really above
        
        #TODO get the index (and from that the pressure) where the density difference is bigger than 0.01
        #TODO still yields negative Indices
        #BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:density_grid[i,:].size + nan_index])>minimal_density_difference))
        
        #
        
        #get the index (and from that the pressure) where the density difference is at maximum (in the lowermost 15 m)
        BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:nan_index]))) -1 
        
        assert(BBL_index>=0)
        print(nan_index,BBL_index)
        ((density_grid[i,BBL_index]-density_grid[i,BBL_index+1]) < minimal_density_difference)
        
        
        #check if the maximum is at the edge of the intervall or if the maximum is too small
        if (BBL_index == BBL_boundary_index) or (BBL_index == (BBL_boundary_index+1)) or ((density_grid[i,BBL_index]-density_grid[i,BBL_index+1]) < minimal_density_difference):
            BBL_index = nan_index #equivalent to a BBL thickness of 0
        

        
        #print(BBL_index,nan_index)
        #print("BBL",interp_pressure[BBL_index])
        #print(bathymetrie[i])
       
        list_of_BBL_indices[i] = BBL_index 
        BBL[i] = interp_pressure[BBL_index]
        
        list_of_BBL_range_indices[i] = BBL_boundary_index 
        BBL_range[i] = interp_pressure[BBL_boundary_index]
        
        
    return [[bathymetrie,list_of_bathymetrie_indices],[BBL,list_of_BBL_indices],[BBL_range,list_of_BBL_range_indices]]


#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################
            
#Constants
rho_0 = 1000 #kg/m³
g = 9.81 #m/s² #could be replace by gsw.grav(lat,p)

datafile_path = "/home/ole/share-windows/emb217_mss_data/TR1-4.mat"

splitted_filename = datafile_path.split("/")
cruisename = splitted_filename[4][0:6]
DATAFILENAME = splitted_filename[-1]
print("cruisename",cruisename)  

#define the pictures
f1, axarr1 = plt.subplots(nrows = 1, ncols = 2, sharey = True)



    
print("Filename:",sio.whosmat(datafile_path))


results = thesis.load_clean_and_interpolate_data(datafile_path)

try:
    number_of_profiles,lat,lon,distance = results[0]
    interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid = results[1]
    eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid = results[2]
except TypeError:
    print(cruisename,DATAFILENAME[:-4],"is skipped!")
            
          


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

#append the last distance plus the last difference (for plotting all the n profiles we need a distance array of size n+1 
plotmesh_distance = np.append(distance,2*distance[-1]-distance[-2])
plotmesh_longitude = np.append(lon,2*lon[-1]-lon[-2])
        
    
#Plotting
transect_index = -10
print("Profile at Longitude",lon[transect_index])

img1_1a = axarr1[0].plot(density_grid[transect_index,:],interp_pressure[:])
img1_1b = axarr1[0].plot(sorted(density_grid[transect_index,:]),interp_pressure[:])

img1_2 = axarr1[1].plot(density_grid[transect_index,:]-sorted(density_grid[transect_index,:]),interp_pressure[:])

f1.set_size_inches(18,10.5)

axarr1[0].invert_yaxis()

f1.suptitle("Thorpe scale")

f1.tight_layout()  
plt.show()


