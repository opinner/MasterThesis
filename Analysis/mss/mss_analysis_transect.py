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

#diff plus adding a row of column of nan values (by diff in axis = 1, the size of axis 1 is conserved) 
def shape_preserving_diff(array,axis):
    concatenated_array = np.diff(array,axis)
    
    nan_array = np.nan*np.ones((np.shape(array)[0],1))

    new_array = np.concatenate((nan_array,concatenated_array),axis = axis)
    
    return new_array

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

FILENAME = "/home/ole/share-windows/emb217_mss_data/TR1-1.mat"

#define the pictures
f1, axarr1 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
f2, axarr2 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
f3, axarr3 = plt.subplots(nrows = 3, ncols = 1)
f4, axarr4 = plt.subplots(nrows = 1, ncols = 7, sharey = True)#, sharex = True, 

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

number_of_transects = np.shape(pressure)[-1]

latitude = []
longitude = []

distance = np.zeros(number_of_transects)
origin = (float(lat[0][0][0]),float(lon[0][0][0])) #lots of brackets to get a number, not an array (workaround)
for i in range(number_of_transects):
    current_point = (float(lat[i][0][0]),float(lon[i][0][0]))
    latitude.append(float(lat[i][0][0]))
    longitude.append(float(lon[i][0][0]))
    distance[i] = geo.geodesic(origin,current_point).km #Distance in km, change to nautical miles?

lat = np.asarray(latitude)
lon = np.asarray(longitude)


print("lat",max(lat),min(lat))
print("lon",max(lon),min(lon))
#initial values 
min_pressure = 10
max_pressure = 60
max_size = 1000
min_size = 3000

#eps profile has a coarser resolution
min_eps_pressure = 10
max_eps_pressure = 60
max_eps_size = 100
min_eps_size = 400

#select the start and end point for the
for i in range(number_of_transects):
    if np.nanmin(pressure[i]) < min_pressure:
        min_pressure = np.nanmin(pressure[i])
    if np.nanmax(pressure[i]) > max_pressure:
        max_pressure = np.nanmax(pressure[i])
    if pressure[i].size > max_size:
        max_size = pressure[i].size       
    if pressure[i].size < min_size:
        min_size = pressure[i].size  

    if np.nanmin(eps_pressure[i]) < min_eps_pressure:
        min_eps_pressure = np.nanmin(eps_pressure[i])
    if np.nanmax(eps_pressure[i]) > max_eps_pressure:
        max_eps_pressure = np.nanmax(eps_pressure[i])
    if eps_pressure[i].size > max_eps_size:
        max_eps_size = eps_pressure[i].size       
    if eps_pressure[i].size < min_eps_size:
        min_eps_size = eps_pressure[i].size  

print("High resolution",min_pressure,max_pressure,min_size,max_size)
print("Coarse resolution",min_eps_pressure,max_eps_pressure,min_eps_size,max_eps_size)

#check if that worked correctly
assert(max_size>= min_size)
assert(max_eps_size>= min_eps_size)


    
interp_pressure = np.linspace(min_pressure,max_pressure,min_size)

#creates pressure grid, where every column is equal to interp_pressure
pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],min_size))

interp_coarse_pressure = np.linspace(min_eps_pressure,max_eps_pressure,min_eps_size)


#create grids with distance on x and depth on y-axis
oxygen_grid = np.zeros((np.shape(pressure)[-1],min_size))
salinity_grid = np.copy(oxygen_grid)
consv_temperature_grid = np.copy(oxygen_grid)
alpha_grid = np.copy(consv_temperature_grid)
beta_grid = np.copy(salinity_grid)

#averaged of approx 5 depth bins (???)
eps_grid = np.zeros((np.shape(pressure)[-1],min_eps_size))


for i in range(number_of_transects): 
    #print(np.shape(oxygen_grid),np.shape(interp_pressure),np.shape(pressure[i]),np.shape(oxygen[i]),np.shape(eps[i].flatten()))
    oxygen_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),oxygen[i].flatten(), left = np.nan, right = np.nan)
    salinity_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
    consv_temperature_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
    alpha_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),alpha[i].flatten(), left = np.nan, right = np.nan)
    beta_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),beta[i].flatten(), left = np.nan, right = np.nan)
    
    eps_grid[i] = np.interp(interp_coarse_pressure,eps_pressure[i].flatten(),eps[i].flatten(), left = np.nan, right = np.nan)
    

density_grid = gsw.rho(salinity_grid,consv_temperature_grid,pressure_grid)
#density_grid_check = (1 - alpha_grid * (consv_temperature_grid) + beta_grid * (salinity_grid))*rho_0
#difference = density_grid-density_grid_check


#creates pressure grid, where every column is equal to interp_pressure
pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],min_size))

#calculate N^2 with the gsw toolbox 
BN_freq_squared_grid_gsw, midpoint_pressure_grid = gsw.Nsquared(salinity_grid,consv_temperature_grid,pressure_grid, lat = np.mean(lat), axis = 1)

print("pre",np.shape(salinity_grid),"post",np.shape(BN_freq_squared_grid_gsw))
print(np.shape(midpoint_pressure_grid))

mid_point_pressure = np.mean(midpoint_pressure_grid, axis = 0)
print(np.mean(np.std(midpoint_pressure_grid, axis = 0)))
assert(np.all(np.std(midpoint_pressure_grid, axis = 0) <10**(-10)))

#interpolate the temperature and to the eps pressure grid
coarse_consv_temperature_grid = np.zeros(np.shape(eps_grid)) #TODO DO I need the in-situ temperature here?
coarse_N_squared_grid = np.zeros(np.shape(eps_grid))
coarse_density_grid = np.copy(coarse_N_squared_grid)

for i in range(number_of_transects):
     
    #midpoint grid to coarse grid 
    coarse_N_squared_grid[i] = np.interp(interp_coarse_pressure,mid_point_pressure,BN_freq_squared_grid_gsw[i].flatten(), left = np.nan, right = np.nan)
     
    #fine grid to coarse grid
    coarse_consv_temperature_grid[i] = np.interp(interp_coarse_pressure,interp_pressure,consv_temperature_grid[i].flatten(), left = np.nan, right = np.nan)
    coarse_density_grid[i] = np.interp(interp_coarse_pressure,interp_pressure,density_grid[i].flatten(), left = np.nan, right = np.nan)

coarse_viscosity_T_grid = get_viscosity(coarse_consv_temperature_grid)

coarse_viscosity_TS_grid = wiki_viscosity(coarse_consv_temperature_grid)/coarse_density_grid


#all used grids should be defined on the coarse pressure grid
Reynolds_bouyancy_grid = eps_grid/(coarse_viscosity_T_grid*coarse_N_squared_grid) #from Maffioli 2014
Reynolds_bouyancy_TS_grid = eps_grid/(coarse_viscosity_TS_grid*coarse_N_squared_grid)

print(np.max(interp_pressure),np.min(interp_pressure))

#Plotting
transect_index = 30
img4_0 = axarr4[0].plot(oxygen_grid[transect_index,:],interp_pressure)
img4_1 = axarr4[1].plot(salinity_grid[transect_index,:],interp_pressure)

img4_2 = axarr4[2].plot(consv_temperature_grid[transect_index,:],interp_pressure)
img4_2b = axarr4[2].plot(coarse_consv_temperature_grid[transect_index,:],interp_coarse_pressure)

img4_3 = axarr4[3].plot(BN_freq_squared_grid_gsw[transect_index,:],mid_point_pressure, label = "fine grid")
img4_3b = axarr4[3].plot(coarse_N_squared_grid[transect_index,:],interp_coarse_pressure, label = "coarse grid")

img4_4 = axarr4[4].plot(coarse_viscosity_T_grid[transect_index,:]*10**6,interp_coarse_pressure,label = "Ilker")
img4_4b = axarr4[4].plot(coarse_viscosity_TS_grid[transect_index,:]*10**6,interp_coarse_pressure,"--",label = "Wikipedia")
img4_5 = axarr4[5].plot(np.log10(eps_grid[transect_index,:]),interp_coarse_pressure)
img4_6 = axarr4[6].plot(Reynolds_bouyancy_grid[transect_index,:],interp_coarse_pressure,label = "Ilker")
img4_6b = axarr4[6].plot(Reynolds_bouyancy_TS_grid[transect_index,:],interp_coarse_pressure,"--",label = "Wikipedia")

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



img3_0 = axarr3[0].pcolormesh(distance,interp_coarse_pressure,Reynolds_bouyancy_grid.T, vmin = -10, vmax = 1000)
img3_1 = axarr3[1].hist(np.log10(eps_grid.flatten()),bins = 300)
img3_1 = axarr3[1].hist(np.log10(coarse_viscosity_TS_grid.flatten()),bins = 300)
img3_1 = axarr3[1].hist(np.log10(coarse_N_squared_grid.flatten()),bins = 300)
img3_2 = axarr3[2].hist(np.log10(Reynolds_bouyancy_grid.flatten()),bins = 300) 
 
#Plot the data   
img1_0 = axarr1[0].pcolormesh(distance,interp_pressure,oxygen_grid.T)
img1_1 = axarr1[1].pcolormesh(distance,interp_pressure,salinity_grid.T)
img1_2 = axarr1[2].pcolormesh(distance,interp_pressure,consv_temperature_grid.T)
img1_3 = axarr1[3].pcolormesh(distance,interp_coarse_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -10)

img2_0 = axarr2[0].pcolormesh(distance,interp_pressure,density_grid.T)
img2_1 = axarr2[1].pcolormesh(distance,interp_pressure,BN_freq_squared_grid_gsw.T,vmin = 0, vmax = 0.015)
img2_2 = axarr2[2].pcolormesh(distance,interp_coarse_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -10)
img2_3 = axarr2[3].pcolormesh(distance,interp_coarse_pressure,Reynolds_bouyancy_grid.T, vmin = -5, vmax = 100)

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


