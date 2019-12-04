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

        

FILENAME = "/home/ole/share-windows/emb217_mss_data/TR1-1.mat"
#FILENAME = "/home/ole/share-windows/emb177_mss_data/TS1_1.mat"
#FILENAME = "/home/ole/share-windows/emb169_mss_data/MSS055/matlab/TS11.mat"

#define the pictures
f1, axarr1 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)
f2, axarr2 = plt.subplots(nrows = 1, ncols = 2, sharex = True, sharey = True)

print("Filename:",sio.whosmat(FILENAME))

data = sio.loadmat(FILENAME)

STA_substructure = data["STA"]
DATA_substructure = data["DATA"]
MIX_substructure = data["MIX"]
CTD_substructure = data["CTD"]
BOT_substructure = data["BOT"]

print("\nSTA\n",STA_substructure.dtype)
print("\nDATA\n",DATA_substructure.dtype)
print("\nCTD\n",CTD_substructure.dtype)
print("\nMIX\n",MIX_substructure.dtype)
print("\nBOT\n",BOT_substructure.dtype)

lat = STA_substructure["LAT"][0]
lon = STA_substructure["LON"][0]

oxygen_ctd_pressure = CTD_substructure["P"][0]
oxygen_ctd = CTD_substructure["O2"][0]

oxygen_data_pressure = DATA_substructure["P"][0]
oxygen_data = DATA_substructure["rawO2"][0]
#oxygen_data = DATA_substructure["O2"][0]

number_of_transects = np.shape(oxygen_ctd_pressure)[-1]

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

#select the start and end point for the
for i in range(number_of_transects):
    if np.nanmin(oxygen_ctd_pressure[i]) < min_pressure:
        min_pressure = np.nanmin(oxygen_ctd_pressure[i])
    if np.nanmax(oxygen_ctd_pressure[i]) > max_pressure:
        max_pressure = np.nanmax(oxygen_ctd_pressure[i])
    if oxygen_ctd_pressure[i].size > max_size:
        max_size = oxygen_ctd_pressure[i].size       
    if oxygen_ctd_pressure[i].size < min_size:
        min_size = oxygen_ctd_pressure[i].size  


print("Resolution",min_pressure,max_pressure,min_size,max_size)

#check if that worked correctly
assert(max_size>= min_size)


    
interp_pressure = np.linspace(min_pressure,max_pressure,min_size)

#creates pressure grid, where every column is equal to interp_pressure
pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(oxygen_ctd_pressure)[-1],min_size))



#create grids with distance on x and depth on y-axis
oxygen_ctd_grid = np.zeros((np.shape(oxygen_ctd_pressure)[-1],min_size))
oxygen_data_grid = np.zeros((np.shape(oxygen_ctd_pressure)[-1],min_size))


#interpolate the values to a common grid
for i in range(number_of_transects): 
    oxygen_ctd_grid[i] = np.interp(interp_pressure,oxygen_ctd_pressure[i].flatten(),oxygen_ctd[i].flatten(), left = np.nan, right = np.nan)
        
    oxygen_data_grid[i] = np.interp(interp_pressure,oxygen_data_pressure[i].flatten(),oxygen_data[i].flatten(), left = np.nan, right = np.nan)

    


print(np.max(interp_pressure),np.min(interp_pressure))

#Plotting
img1_0 = axarr1[0].pcolormesh(distance,interp_pressure,oxygen_ctd_grid.T)
img1_1 = axarr1[1].pcolormesh(distance,interp_pressure,oxygen_data_grid.T)

axarr1[0].set_title("oxygen ctd")
axarr1[1].set_title("oxygen data")

colorbar(img1_0).set_label("CTD Oxygen [??]") 
colorbar(img1_1).set_label("DATA Oxygen [??]") 

#figure2
transect_index = 30
img2_0 = axarr2[0].plot(oxygen_ctd_grid[transect_index,:],interp_pressure)
img2_1 = axarr2[1].plot(oxygen_data_grid[transect_index,:],interp_pressure)

axarr2[0].set_xlabel("oxygen ctd")
axarr2[1].set_xlabel("oxygen data")
axarr2[0].set_ylabel("pressure [dbar]")
 
   
f1.set_size_inches(18,10.5)
f2.set_size_inches(18,10.5)


axarr1[0].invert_yaxis()

f1.tight_layout() 

plt.show()


