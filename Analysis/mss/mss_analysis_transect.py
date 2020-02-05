#---------------------------------------------------------#
#Plots one mss transect 
#plus as an exampe one profile from that transect
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

data = sio.loadmat(FILENAME)


STA_substructure = data["STA"]
DATA_substructure = data["DATA"]
MIX_substructure = data["MIX"]
CTD_substructure = data["CTD"]

#print(STA_substructure.dtype)
#print(DATA_substructure.dtype)
#print(CTD_substructure.dtype)
#print(MIX_substructure.dtype)

lat = STA_substructure["LAT"][0]
lon = STA_substructure["LON"][0]

#print(lat)

pressure = CTD_substructure["P"][0]
oxygen = CTD_substructure["O2"][0]
absolute_salinity = CTD_substructure["SA"][0] #is this unit sufficient
consv_temperature = CTD_substructure["CT"][0] #TODO better use conservative temperature?
alpha = CTD_substructure["ALPHA"][0]
beta = CTD_substructure["BETA"][0]

#print("alpha",np.shape(alpha),np.shape(alpha[0]))

eps = MIX_substructure["eps"][0]
eps_pressure = MIX_substructure["P"][0]


    

"""
print("lat",np.shape(lat))

print(np.shape(eps_pressure))
print(np.shape(eps))

print(np.shape(np.asarray(pressure)))
"""

number_of_profiles = np.shape(pressure)[-1]

latitude = []
longitude = []

distance = np.zeros(number_of_profiles)
origin = (float(lat[0][0][0]),float(lon[0][0][0])) #lots of brackets to get a number, not an array (workaround)
for i in range(number_of_profiles):
    current_point = (float(lat[i][0][0]),float(lon[i][0][0]))
    latitude.append(float(lat[i][0][0]))
    longitude.append(float(lon[i][0][0]))
    distance[i] = geo.geodesic(origin,current_point).km #Distance in km, change to nautical miles?

lat = np.asarray(latitude)
lon = np.asarray(longitude)


print("lat",max(lat),min(lat))
print("lon",max(lon),min(lon))

#Data Cleaning
#-----------------------------------------------------------------------------------------------------------------
#remove data from a file, with overlapping positional points
if (FILENAME == "/home/ole/share-windows/emb217_mss_data/TR1-8.mat"):
    lat = np.delete(lat,np.s_[33:47])
    lon = np.delete(lon,np.s_[33:47])
    distance = np.delete(distance,np.s_[33:47])
    
    pressure = np.delete(pressure,np.s_[33:47],axis=0)
    oxygen = np.delete(oxygen,np.s_[33:47],axis=0)
    absolute_salinity =  np.delete(absolute_salinity,np.s_[33:47],axis=0)
    consv_temperature = np.delete(consv_temperature,np.s_[33:47],axis=0)
    alpha = np.delete(alpha,np.s_[33:47],axis=0)
    beta = np.delete(beta,np.s_[33:47],axis=0)
    
    
    eps = np.delete(eps,np.s_[33:47],axis=0)
    eps_pressure = np.delete(eps_pressure,np.s_[33:47],axis=0)
    
    number_of_profiles = np.shape(pressure)[-1]

#remove excess data that already belongs to TS119
if cruisename == "emb169" and  DATAFILENAME[:-4] == "TS118":
    lat = lat[:21]
    lon = lon[:21]
    distance = distance[:21]
    
    pressure = pressure[:21]
    oxygen = oxygen[:21]
    absolute_salinity =  absolute_salinity[:21]
    consv_temperature = consv_temperature[:21]
    alpha = alpha[:21]
    beta = beta[:21]
    
    
    eps = eps[:21]
    eps_pressure = eps_pressure[:21]
    number_of_profiles = np.shape(pressure)[-1]

#removes the last data point, that dont seem to belong to the transect
if cruisename == "emb169" and  DATAFILENAME[:-4] == "TRR109":
    lat = lat[:-1]
    lon = lon[:-1]
    distance = distance[:-1]
    
    pressure = pressure[:-1]
    oxygen = oxygen[:-1]
    absolute_salinity =  absolute_salinity[:-1]
    consv_temperature = consv_temperature[:-1]
    alpha = alpha[:-1]
    beta = beta[:-1]
    
    
    eps = eps[:-1]
    eps_pressure = eps_pressure[:-1]
    number_of_profiles = np.shape(pressure)[-1]
            
                        
if cruisename == "emb169" and  DATAFILENAME[:-4] == "TRR109":
    lat = lat[:-1]
    lon = lon[:-1]
    distance = distance[:-1]
    
    pressure = pressure[:-1]
    oxygen = oxygen[:-1]
    absolute_salinity =  absolute_salinity[:-1]
    consv_temperature = consv_temperature[:-1]
    alpha = alpha[:-1]
    beta = beta[:-1]
    
    
    eps = eps[:-1]
    eps_pressure = eps_pressure[:-1]
    number_of_profiles = np.shape(pressure)[-1]
#-----------------------------------------------------------------------------------------------------------------
            
#test if distance is monotonically increasing
assert(np.all(np.diff(distance)>0))

#initial values 
min_pressure = 10
max_pressure = 60
max_size = 1000
min_size = 3000

#select the start and end point for the
for i in range(number_of_profiles):
    if np.nanmin(pressure[i]) < min_pressure:
        min_pressure = np.nanmin(pressure[i])
    if np.nanmax(pressure[i]) > max_pressure:
        max_pressure = np.nanmax(pressure[i])
    if pressure[i].size > max_size:
        max_size = pressure[i].size       
    if pressure[i].size < min_size:
        min_size = pressure[i].size  
"""
    if np.nanmin(eps_pressure[i]) < min_eps_pressure:
        min_eps_pressure = np.nanmin(eps_pressure[i])
    if np.nanmax(eps_pressure[i]) > max_eps_pressure:
        max_eps_pressure = np.nanmax(eps_pressure[i])
    if eps_pressure[i].size > max_eps_size:
        max_eps_size = eps_pressure[i].size       
    if eps_pressure[i].size < min_eps_size:
        min_eps_size = eps_pressure[i].size  
"""
print("High resolution",min_pressure,max_pressure,min_size,max_size)
#print("Coarse resolution",min_eps_pressure,max_eps_pressure,min_eps_size,max_eps_size)

#check if that worked correctly
assert(max_size>= min_size)
#assert(max_eps_size>= min_eps_size)


#creates a pressure axis between the maximum and minimum pressure    
interp_pressure = np.linspace(min_pressure,max_pressure,min_size)

#creates pressure grid, where every column is equal to interp_pressure
pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],min_size))

#create grids that changes in distance on x and in depth on y-axis
oxygen_grid = np.zeros((np.shape(pressure)[-1],min_size))
salinity_grid = np.copy(oxygen_grid)
consv_temperature_grid = np.copy(oxygen_grid)
alpha_grid = np.copy(consv_temperature_grid)
beta_grid = np.copy(salinity_grid)

#check if the pressure points for every eps profile are the same
for i in range(number_of_profiles):  
    assert(np.all(eps_pressure[i].flatten() == eps_pressure[0].flatten()))
#if yes the pressure can be a 1D array instead of a 2D array    
eps_pressure = eps_pressure[0].flatten()
    
#averaged of approx 5 depth bins (???)
eps_grid = np.zeros((number_of_profiles,eps_pressure.size))

#needed for the interpolation of S and T to the same grid as the eps
eps_salinity_grid = np.ones((np.shape(eps_grid)))
eps_consv_temperature_grid = np.ones(np.shape(eps_grid))

#vector times matrix multiplication to get a 2D array, where every column is equal to eps_pressure
eps_pressure_grid = np.reshape(eps_pressure,(1,-1))*np.ones(np.shape(eps_grid))

#create a pressure axis where every point is shifted by half the distance to the next one
shifted_pressure = eps_pressure + np.mean(np.diff(eps_pressure))/2

#prepend a point at the beginning to be a point longer than the pressure axis we want from the Nsquared function
shifted_pressure = np.append(eps_pressure[0]-np.mean(np.diff(eps_pressure))/2, shifted_pressure)

#from that create a grid, where we have just n-times the pressure axis with n the number of profiles 
shifted_pressure_grid = np.reshape(shifted_pressure,(1,-1))*np.ones((number_of_profiles,shifted_pressure.size))
shifted_salinity_grid = np.ones((np.shape(shifted_pressure_grid)))
shifted_consv_temperature_grid = np.ones((np.shape(shifted_pressure_grid)))




for i in range(number_of_profiles): 
    #interpolation to a common fine grid
    oxygen_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),oxygen[i].flatten(), left = np.nan, right = np.nan)
    salinity_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
    consv_temperature_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
    alpha_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),alpha[i].flatten(), left = np.nan, right = np.nan)
    beta_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),beta[i].flatten(), left = np.nan, right = np.nan)
        
    #interpolation to a shifted grid for the Nsquared function
    shifted_salinity_grid[i] = np.interp(shifted_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
    shifted_consv_temperature_grid[i] = np.interp(shifted_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
    
    #just changes the format of eps slightly
    assert(eps[i].flatten().size == eps_pressure.size)
    eps_grid[i] = eps[i].flatten()

    #interpolate S and T to the same grid as eps
    eps_salinity_grid[i] = np.interp(eps_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
    eps_consv_temperature_grid[i] = np.interp(eps_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
    
#fine density grid
density_grid = gsw.rho(salinity_grid,consv_temperature_grid,pressure_grid)

#density grid on the same points as the eps grid
eps_density_grid = gsw.rho(eps_salinity_grid,eps_consv_temperature_grid,eps_pressure_grid)

#TODO compare with density_grid?
density_grid_check = (1 - alpha_grid * (consv_temperature_grid) + beta_grid * (salinity_grid))*rho_0
difference = density_grid-density_grid_check

#TODO calculate N^2 with the gsw toolbox 
#BV_freq_squared_grid_gsw, midpoint_pressure_grid = gsw.Nsquared(salinity_grid,consv_temperature_grid,pressure_grid, lat = np.mean(lat), axis = 1)

#calculate N^2 with the gsw toolbox, by using teh shifted grid we should get a N^2 grid that is defined at teh same points as eps_grid
eps_N_squared, crosscheck_pressure_grid = gsw.Nsquared(shifted_salinity_grid,shifted_consv_temperature_grid,shifted_pressure_grid, lat = np.mean(lat), axis = 1)
crosscheck_pressure = np.mean(crosscheck_pressure_grid, axis = 0)

#test if we really calculated N^2 for the same points as the dissipation measurement
assert(np.all(crosscheck_pressure == eps_pressure))

print("------------------TEST--------------------------")
print("eps_pressure\n",np.shape(eps_pressure))

print("shifted_pressure\n",np.shape(shifted_pressure))
print("crosscheck_pressure\n",np.shape(crosscheck_pressure))



#mid_point_pressure = np.mean(midpoint_pressure_grid, axis = 0)
#print(np.mean(np.std(midpoint_pressure_grid, axis = 0)))
#assert(np.all(np.std(midpoint_pressure_grid, axis = 0) <10**(-10))) #midpoints should be every time the same

#interpolate the temperature and to the eps pressure grid
#coarse_consv_temperature_grid = np.zeros(np.shape(eps_grid)) #TODO DO I need the in-situ temperature here?
#coarse_N_squared_grid = np.zeros(np.shape(eps_grid))
#coarse_density_grid = np.copy(coarse_N_squared_grid)

#for i in range(number_of_profiles):
     
    #midpoint grid to coarse grid 
    #coarse_N_squared_grid[i] = np.interp(interp_coarse_pressure,mid_point_pressure,BV_freq_squared_grid_gsw[i].flatten(), left = np.nan, right = np.nan)
     
    #fine grid to coarse grid
    #coarse_consv_temperature_grid[i] = np.interp(interp_coarse_pressure,interp_pressure,consv_temperature_grid[i].flatten(), left = np.nan, right = np.nan)
    #coarse_density_grid[i] = np.interp(interp_coarse_pressure,interp_pressure,density_grid[i].flatten(), left = np.nan, right = np.nan)


eps_wiki_viscosity_grid = wiki_viscosity(eps_consv_temperature_grid)/eps_density_grid
eps_viscosity_grid = get_viscosity(eps_consv_temperature_grid)

print("eps_viscosity_grid\n",np.shape(eps_viscosity_grid))

eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared)
eps_wiki_Reynolds_bouyancy_grid = eps_grid/(eps_wiki_viscosity_grid*eps_N_squared)

#TODO delete?
#np.testing.assert_equal(eps_grid,shifted_eps_grid)
#assert(np.all(interp_coarse_pressure == eps_pressure))



#print(np.max(interp_pressure),np.min(interp_pressure))

#search for bottom currents
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
    nan_index =  -np.argmax(np.flip(~np.isnan(salinity_grid[i,:]))) #at the moment the index is negative
    nan_index = density_grid[i,:].size + nan_index #now defined as positive index
    
   
    if nan_index == density_grid[i,:].size:
        if not np.isnan(salinity_grid[i,-1]): #if there are no NAN values towards the bottom
            nan_index = len(interp_pressure)-1
            
    list_of_bathymetrie_indices[i] = nan_index 
    bathymetrie[i] = interp_pressure[nan_index]
    
    #print(nan_index)
    #print("bathymetrie check\n")
    #print(interp_pressure[nan_index-1],salinity_grid[i,nan_index-1])
    #print(interp_pressure[nan_index],salinity_grid[i,nan_index])
    #print(interp_pressure[nan_index+1],salinity_grid[i,nan_index+1])
    #while interp_coarse_pressure[j]
    #if 
  
    #TODO
    #------------------search for halocline values starts from above:-------------------------------
    


    #------------------search for BBL values starts from below:-------------------------------  
    
    #index of bottom plus 15m 
    #print("Surface",interp_pressure[0],"Bottom",bathymetrie[i])
    height_above_ground = 10
    BBL_boundary_index = np.argmax(interp_pressure >= (bathymetrie[i]-height_above_ground))
    assert(interp_pressure[BBL_boundary_index]<bathymetrie[i]) #tests if the point 15m above the ground is really above
    #COMMENT
    #print("Max BBL layer",interp_pressure[BBL_boundary_index],bathymetrie[i], density_grid[i,:].size + nan_index)
    #print(np.shape(density_grid[i,BBL_boundary_index:density_grid[i,:].size + nan_index]))
    #get the index (and from that the pressure) where the density difference is bigger than 0.01
    #BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:density_grid[i,:].size + nan_index])>0.01))
    
    #get the index (and from that the pressure) where the density difference is at maximum (in the lowermost 15 m)
    assert(nan_index>=0)
    BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:nan_index]))) -1 
    
    #check if the maximum is at the edge of the intervall or if the maximum is too small
    minimal_density_difference = 0.02
    if (BBL_index == BBL_boundary_index) or (BBL_index == (BBL_boundary_index+1)) or ((density_grid[i,BBL_index]-density_grid[i,BBL_index-1]) < minimal_density_difference):
        BBL_index = nan_index #equivalent to a BBL thickness of 0
    
    #print(BBL_index,nan_index)
    #print("BBL",interp_pressure[BBL_index])
    #print(bathymetrie[i])
   
    list_of_BBL_indices[i] = BBL_index 
    BBL[i] = interp_pressure[BBL_index]
    
    list_of_BBL_range_indices[i] = BBL_boundary_index 
    BBL_range[i] = interp_pressure[BBL_boundary_index]

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


