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


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)


def get_viscosity(T):
#% T is temperature in degC; vis in m2/s
#% vis=(1.792747-(.05126103*T)+(0.0005918645*T*T))*1e-6;
#% Ilker
    return (1.792747-(0.05126103*T)+(0.0005918645*T*T))*1e-6

def wiki_viscosity(T):
    A = 29.39*1e-3
    B = 507.88 
    C = 149.3
    
    return A * np.exp(B/(T-C))



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
                
        list_of_bathymetrie_indices[i] = nan_index 
        bathymetrie[i] = interp_pressure[nan_index]
        
      
        #TODO
        #------------------search for halocline values starts from above:-------------------------------
        


        #------------------search for BBL values starts from below:-------------------------------  
        
        #index of maximal distance bottom plus 15m 
        
        BBL_boundary_index = np.argmax(interp_pressure >= (bathymetrie[i]-height_above_ground))
        assert(interp_pressure[BBL_boundary_index]<bathymetrie[i]) #tests if the point 15m above the ground is really above
        
        #TODO get the index (and from that the pressure) where the density difference is bigger than 0.01
        BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:density_grid[i,:].size + nan_index])>minimal_density_difference))
        
        #get the index (and from that the pressure) where the density difference is at maximum (in the lowermost 15 m)
        assert(nan_index>=0)
        #BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:nan_index]))) -1 
        
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

FILENAME = "/home/ole/share-windows/emb217_mss_data/TR1-4.mat"

splitted_filename = FILENAME.split("/")
cruisename = splitted_filename[4][0:6]
print("cruisename",cruisename)  

#define the pictures
f1, axarr1 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
f2, axarr2 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
f3, axarr3 = plt.subplots(nrows = 3, ncols = 1)
f4, axarr4 = plt.subplots(nrows = 1, ncols = 8, sharey = True)#, sharex = True, 
f5, axarr5 = plt.subplots(2)


    
print("Filename:",sio.whosmat(FILENAME))


results = thesis.load_clean_and_interpolate_data(datafile_path)

try:
    number_of_profiles,lat,lon,distance = results[0]
    interp_pressure,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid = results[1]
    eps_pressure,eps_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid = results[2]
except TypeError:
    print(cruisename,DATAFILENAME[:-4],"is skipped!")
    continue
            
          

#calculate the viscosity (2 different formula)
eps_wiki_viscosity_grid = wiki_viscosity(eps_consv_temperature_grid)/eps_density_grid
eps_viscosity_grid = get_viscosity(eps_consv_temperature_grid)

#calculate the Reynolds bouyancy number defined on eps_pressure
eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared_grid)
eps_wiki_Reynolds_bouyancy_grid = eps_grid/(eps_wiki_viscosity_grid*eps_N_squared_grid)

results = thesis.find_bottom_and_bottom_currents(number_of_profiles,interp_pressure,density_grid,oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02)
bathymetrie,list_of_bathymetrie_indices = results[0]
BBL,list_of_BBL_indices = results[1]
BBL_range,list_of_BBL_range_indices = results[2]
        
#Test of different BBL conditions
exp_results = experimental_BBL(number_of_profiles,interp_pressure,density_grid,oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02)
exp_bathymetrie,exp_list_of_bathymetrie_indices = results[0]
exp_BBL,exp_list_of_BBL_indices = results[1]
exp_BBL_range,exp_list_of_BBL_range_indices = results[2]

#append the last distance plus the last difference (for plotting all the n profiles we need a distance array of size n+1 
plotmesh_distance = np.append(distance,2*distance[-1]-distance[-2])
       
axarr5[0].plot(distance,bathymetrie)
axarr5[0].invert_yaxis()
axarr5[1].plot(lon)    
    
#draw the calculated layers in the plot    
axarr1[0].plot(distance,bathymetrie)
axarr1[0].plot(distance,BBL)
    
#Plotting
transect_index = -5
print("Profile at Longitude",lon[transect_index])
#print(lon)
#print(np.all(np.diff(lon)>0))
img4_0 = axarr4[0].plot(oxygen_grid[transect_index,:],interp_pressure)
img4_0b = axarr4[0].plot(0,BBL[transect_index],"Dr")
img4_0c = axarr4[0].plot(0,bathymetrie[transect_index],"Dg")
img4_0d = axarr4[0].plot(0,BBL_range[transect_index],"ok")
print("Bottom",bathymetrie[transect_index],"BBL",BBL[transect_index],"max BBL",interp_pressure[transect_index],)


img4_1 = axarr4[1].plot(salinity_grid[transect_index,:],interp_pressure)

img4_2 = axarr4[2].plot(consv_temperature_grid[transect_index,:],interp_pressure)
img4_2b = axarr4[2].plot(eps_consv_temperature_grid[transect_index,:],eps_pressure)

#img4_3 = axarr4[3].plot(BV_freq_squared_grid_gsw[transect_index,:],mid_point_pressure, label = "fine grid")
img4_3b = axarr4[3].plot(np.log10(eps_N_squared[transect_index,:]),eps_pressure, label = "eps grid")

img4_4 = axarr4[4].plot(eps_viscosity_grid[transect_index,:]*10**6,eps_pressure,label = "Ilker")
img4_4b = axarr4[4].plot(eps_wiki_viscosity_grid[transect_index,:]*10**6,eps_pressure,"--",label = "Wikipedia")
img4_5 = axarr4[5].plot(np.log10(eps_grid[transect_index,:]),eps_pressure)
img4_6 = axarr4[6].plot(eps_Reynolds_bouyancy_grid[transect_index,:],eps_pressure,label = "Ilker")
img4_6b = axarr4[6].plot(eps_wiki_Reynolds_bouyancy_grid[transect_index,:],eps_pressure,"--",label = "Wikipedia")

print(np.shape(interp_pressure))
print(np.shape(salinity_grid))
img4_7 = axarr4[7].plot(np.diff(density_grid[transect_index,:]),interp_pressure[1:])
#img4_7b = axarr4[7].plot(BV_freq_squared_grid_gsw[transect_index,BBL_boundary_index:density_grid[transect_index,:].size + nan_index],mid_point_pressure[BBL_boundary_index:density_grid[transect_index,:].size + nan_index])

axarr4[0].set_xlabel("oxygen")
axarr4[0].set_ylabel("pressure [dbar]")
axarr4[1].set_xlabel("SA")
axarr4[2].set_xlabel("CT")
axarr4[3].set_xlabel("N^2")
axarr4[4].set_xlabel("viscosity / 10^(-6)")
axarr4[5].set_xlabel("log10(eps)")
axarr4[6].set_xlabel("Reynolds_bouyancy")
axarr4[3].legend()
axarr4[4].legend()
axarr4[6].legend()

axarr4[6].set_xlim([-15,100])

axarr4[0].set_title("Measurement 02")
axarr4[1].set_title("Measurement SA")
axarr4[2].set_title("Measurement CT")
axarr4[3].set_title("Calculation N^2")
axarr4[4].set_title("Calculation nu")
axarr4[5].set_title("Measurement eps")
axarr4[6].set_title("Calculation Re_b")



img3_0 = axarr3[0].pcolormesh(lon,eps_pressure,eps_Reynolds_bouyancy_grid.T, vmin = -10, vmax = 1000)
img3_1 = axarr3[1].hist(np.log10(eps_grid.flatten()),bins = 300)
img3_1 = axarr3[1].hist(np.log10(eps_viscosity_grid.flatten()),bins = 300)
img3_1 = axarr3[1].hist(np.log10(eps_N_squared.flatten()),bins = 300)
img3_2 = axarr3[2].hist(np.log10(eps_Reynolds_bouyancy_grid.flatten()),bins = 300) 
 
#Plot the data   
img1_0 = axarr1[0].pcolormesh(plotmesh_distance,interp_pressure,oxygen_grid.T)
img1_1 = axarr1[1].pcolormesh(plotmesh_distance,interp_pressure,salinity_grid.T)
img1_2 = axarr1[2].pcolormesh(plotmesh_distance,interp_pressure,consv_temperature_grid.T)
img1_3 = axarr1[3].pcolormesh(plotmesh_distance,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -10)

img2_0 = axarr2[0].pcolormesh(plotmesh_distance,interp_pressure,density_grid.T)
img2_1 = axarr2[1].pcolormesh(plotmesh_distance,eps_pressure,eps_N_squared.T,vmin = 0, vmax = 0.015)
img2_2 = axarr2[2].pcolormesh(plotmesh_distance,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -10)
img2_3 = axarr2[3].pcolormesh(plotmesh_distance,eps_pressure,eps_Reynolds_bouyancy_grid.T, vmin = -5, vmax = 100)

axarr1[3].set_xlabel("distance [km]")
axarr2[3].set_xlabel("distance [km]")

    
    
f1.set_size_inches(18,10.5)
f2.set_size_inches(18,10.5)
f3.set_size_inches(18,10.5)
f4.set_size_inches(18,10.5)

colorbar(img1_0).set_label("Oxygen [??]") 
colorbar(img1_1).set_label("salinity [SA]") 
colorbar(img1_2).set_label("consv_temperature [C]")
colorbar(img1_3).set_label("log10(dissipation) [??]")

colorbar(img2_0).set_label(r"density [kg/$m^3$]")
colorbar(img2_1).set_label(r"$N^2$ $[1/s^2]$")
colorbar(img2_2).set_label(r"log10($\epsilon$) [??]")  
colorbar(img2_3).set_label(r"$Re_b$")

colorbar(img3_0)



axarr1[0].invert_yaxis()
axarr2[0].invert_yaxis()
axarr3[0].invert_yaxis() 
axarr4[0].invert_yaxis() 

f1.suptitle("TEST")
f2.suptitle("TEST")

f1.tight_layout() 
f2.tight_layout() 
f4.tight_layout()  
plt.show()


