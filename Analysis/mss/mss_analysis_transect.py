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
    
#Constants
rho_0 = 1000 #kg/m³
g = 9.81 #m/s² #could be replace by gsw.grav(lat,p)

FILENAME = "/home/ole/Thesis/emb217_mss_data/TR1-1.mat"

#define the pictures
f1, axarr1 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
f2, axarr2 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
f3, axarr3 = plt.subplots(nrows = 3, ncols = 1)
f4, axarr4 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)

print(sio.whosmat(FILENAME))

data = sio.loadmat(FILENAME)


STA_substructure = data["STA"]
DATA_substructure = data["DATA"]
MIX_substructure = data["MIX"]
CTD_substructure = data["CTD"]

print(STA_substructure.dtype)
print(DATA_substructure.dtype)
print(CTD_substructure.dtype)
print(MIX_substructure.dtype)

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
MIX_pressure = MIX_substructure["P"][0]


print("lat",np.shape(lat))

print(np.shape(MIX_pressure))
print(np.shape(eps))

print(np.shape(np.asarray(pressure)))

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

min_pressure = 10
max_pressure = 60
max_size = 1000
for i in range(number_of_transects):
    if np.nanmin(pressure[i]) < min_pressure:
        min_pressure = np.nanmin(pressure[i])
    if np.nanmax(pressure[i]) > max_pressure:
        max_pressure = np.nanmax(pressure[i])
    if pressure[i].size > max_size:
        max_size = pressure[i].size       

for i in range(40):
    pass
    #axarr2.plot(np.log10(eps[i]),MIX_pressure[i])
    
print(min_pressure,max_pressure,max_size)
interp_pressure = np.linspace(min_pressure,max_pressure,max_size)

#creates pressure grid, where every column is equal to interp_pressure
pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],max_size))

eps_pressure = np.linspace(min_pressure,max_pressure,319) #TODO no hardcoding


#create grids with distance on x and depth on y-axis
oxygen_grid = np.zeros((np.shape(pressure)[-1],max_size))
salinity_grid = np.copy(oxygen_grid)
consv_temperature_grid = np.copy(oxygen_grid)
alpha_grid = np.copy(consv_temperature_grid)
beta_grid = np.copy(salinity_grid)

#averaged of approx 5 depth bins (???)
eps_grid = np.zeros((np.shape(pressure)[-1],319)) #no hardcoding


for i in range(number_of_transects): 
    #print(np.shape(oxygen_grid),np.shape(interp_pressure),np.shape(pressure[i]),np.shape(oxygen[i]),np.shape(eps[i].flatten()))
    oxygen_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),oxygen[i].flatten(), left = np.nan, right = np.nan)
    salinity_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
    consv_temperature_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
    alpha_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),alpha[i].flatten(), left = np.nan, right = np.nan)
    beta_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),beta[i].flatten(), left = np.nan, right = np.nan)
    
    eps_grid[i] = np.interp(eps_pressure,MIX_pressure[i].flatten(),eps[i].flatten(), left = np.nan, right = np.nan)
    

density_grid = gsw.rho(salinity_grid,consv_temperature_grid,pressure_grid)
density_grid_check = (1 - alpha_grid * (consv_temperature_grid) + beta_grid * (salinity_grid))*rho_0
difference = density_grid-density_grid_check


#creates pressure grid, where every column is equal to interp_pressure
pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],max_size))

#TODO WHat is the best way to calculate the discrete derivation? King used centered scheme
#discrete derivation of the density in respect to pressure
drho_dp = shape_preserving_diff(density_grid,axis = 1)/shape_preserving_diff(pressure_grid,axis = 1)

print(np.shape(drho_dp))

BN_freq_squared_grid = g/rho_0 * drho_dp

#why does this give bad results????
#BN_freq_squared_grid_check = (g**2)*(shape_preserving_diff(density_grid_check,axis = 1) - 1/(gsw.sound_speed(salinity_grid,consv_temperature_grid,pressure_grid)**2))

BN_freq_squared_grid_check = g*((1/(pressure_grid)*(shape_preserving_diff(density_grid_check,axis = 1))) + g/(gsw.sound_speed(salinity_grid,consv_temperature_grid,pressure_grid)**2))

BN_freq_squared_grid_gsw, dump = gsw.Nsquared(salinity_grid,consv_temperature_grid,pressure_grid, lat = np.mean(lat), axis = 1)
nan_array = np.nan*np.ones((np.shape(BN_freq_squared_grid_gsw)[0],1))
BN_freq_squared_grid_gsw = np.concatenate((nan_array,BN_freq_squared_grid_gsw),axis = 1)
    
    
difference_grid = BN_freq_squared_grid-BN_freq_squared_grid_gsw

#TODO
BN_freq_squared_cleaned_grid = np.copy(BN_freq_squared_grid)
BN_freq_squared_cleaned_grid_check = np.copy(BN_freq_squared_grid_check)

#replace all negative values with 0 (is that good?) to compute the frequency
BN_freq_squared_cleaned_grid[BN_freq_squared_grid < 0] = np.nan #replace all negative values with 0 (is that good?)
BN_freq_squared_cleaned_grid_check[BN_freq_squared_grid_check < 0]
#draw the square root
BN_freq_cleaned_grid = np.sqrt(BN_freq_squared_cleaned_grid)
BN_freq_cleaned_grid_check = np.sqrt(BN_freq_squared_cleaned_grid_check)
#Reynolds_bouyancy_grid = eps_grid/(viscosity_grid*BN_freq_squared_grid) #from Maffioli 2014


img4_0 = axarr4[0].pcolormesh(distance,interp_pressure,BN_freq_squared_grid.T, vmin = 0, vmax = 0.01)
img4_1 = axarr4[1].pcolormesh(distance,interp_pressure,BN_freq_squared_grid_gsw.T, vmin = 0, vmax = 0.01)


img3_0 = axarr3[0].pcolormesh(distance,interp_pressure,difference_grid.T)#, vmin = 0, vmax = 0.06)

img3_1 = axarr3[1].hist(BN_freq_squared_grid.flatten(),bins = 300)
img3_2 = axarr3[2].hist(difference_grid.flatten(),bins = 300) 
img3_2b = axarr3[1].hist(BN_freq_squared_grid_gsw.flatten(),bins = 300) 
 
 
#Plot the data   
img1_0 = axarr1[0].pcolormesh(distance,interp_pressure,oxygen_grid.T)
img1_1 = axarr1[1].pcolormesh(distance,interp_pressure,density_grid.T)
#img1_2 = axarr1[2].pcolormesh(distance,interp_pressure,BN_freq_squared_grid.T)#,vmin = 0, vmax = 0.06)
img1_2 = axarr1[2].pcolormesh(distance,interp_pressure,BN_freq_squared_grid_gsw.T,vmin = 0, vmax = 0.015)
img1_3 = axarr1[3].pcolormesh(distance,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -10)


img2_0 = axarr2[0].pcolormesh(distance,interp_pressure,salinity_grid.T)
img2_1 = axarr2[1].pcolormesh(distance,interp_pressure,consv_temperature_grid.T)
img2_2 = axarr2[2].pcolormesh(distance,interp_pressure,alpha_grid.T)
img2_3 = axarr2[3].pcolormesh(distance,interp_pressure,beta_grid.T)


f1.set_size_inches(18,10.5)
f2.set_size_inches(18,10.5)
f3.set_size_inches(18,10.5)
f4.set_size_inches(18,10.5)

colorbar(img1_0).set_label("Oxygen [??]") 
colorbar(img1_1).set_label(r"density [kg/$m^3$]")
colorbar(img1_2).set_label(r"$N^2$ $[1/s^2]$")
colorbar(img1_3).set_label("log10(dissipation) [??]")  

colorbar(img2_0).set_label("salinity [SA]") 
colorbar(img2_1).set_label("consv_temperature [C]") 
colorbar(img2_2).set_label("alpha [??]")
colorbar(img2_3).set_label("beta [??]")

colorbar(img3_0)

colorbar(img4_0)
colorbar(img4_1)


axarr1[0].invert_yaxis()
axarr2[0].invert_yaxis()
axarr3[0].invert_yaxis() 
axarr4[0].invert_yaxis()   
plt.show()


