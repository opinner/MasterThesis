#---------------------------------------------------------------------------#
#Creates a multidem#
#---------------------------------------------------------------------------#

import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import datetime as dt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit 
from mpl_toolkits.basemap import Basemap

def custom_colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)



cmap = plt.get_cmap('seismic')
cmap.set_bad(color = 'lightgrey')#, alpha = 0.5)
  
##########################################
# Load ADCD data
#########################################    

#Load TC-flach
#-------------------------    
#print(sio.whosmat(FILENAME))
datafile_path = "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP1200/data/EMB217_TC-flach_adcp1200_val.mat"
data = sio.loadmat(datafile_path)
data = data["adcpavg"]
substructure = data.dtype

rtc = data["rtc"][0][0].flatten()
curr = data["curr"][0][0]

trim_emb217_flach = [0,8086]
rtc = rtc[0:8086]
curr = curr[:,0:8086]
    
all_west_east_flach = np.real(curr).T
all_north_south_flach = np.imag(curr).T

#interpolate to fill in missing data (possible due to a constant timestep)
index_adcp_flach = 29
temp_x = all_west_east_flach[index_adcp_flach,:]
temp_y = all_north_south_flach[index_adcp_flach,:]
xi = np.arange(len(temp_x))
mask = np.isfinite(temp_x)
west_east_flach = np.interp(xi, xi[mask], temp_x[mask])
north_south_flach = np.interp(xi, xi[mask], temp_y[mask])

#convert matlab time to utc
utc_flach = np.asarray(pl.num2date(rtc-366))

depth_flach = (data["depth"][0][0]).flatten()
number_of_depth_bins_flach = depth_flach.size
assert(np.mean(depth_flach)>0)

#Load TC-tief
#-------------------------    
#print(sio.whosmat(FILENAME))
datafile_path = "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/adcp/data/EMB217_TC-tief_adcp300_val.mat"
data = sio.loadmat(datafile_path)
data = data["adcpavg"]
substructure = data.dtype

rtc = data["rtc"][0][0].flatten()
curr = data["curr"][0][0]

trim_emb217_tief = [100,8165]
rtc = rtc[100:8165]
    
all_west_east_tief = np.real(curr).T
all_north_south_tief = np.imag(curr).T

all_north_south_tief = all_north_south_tief[:,100:8165]
all_west_east_tief = all_west_east_tief[:,100:8165] 

#interpolate to fill in missing data (possible due to a constant timestep)
index_adcp_tief = 7
print(np.shape(all_west_east_tief))
temp_x = all_west_east_tief[index_adcp_tief,:]
temp_y = all_north_south_tief[index_adcp_tief,:]
xi = np.arange(len(temp_x))
mask = np.isfinite(temp_x)
west_east_tief = np.interp(xi, xi[mask], temp_x[mask])
north_south_tief = np.interp(xi, xi[mask], temp_y[mask])

print(np.shape(west_east_flach))
print(np.shape(west_east_tief))

#convert matlab time to utc
utc_tief = np.asarray(pl.num2date(rtc-366))

depth_tief = (data["depth"][0][0]).flatten()
number_of_depth_bins_tief = depth_tief.size
assert(np.mean(depth_tief)>0)

"""
for i,d in enumerate(depth_flach):
    print(i,d)
for i,d in enumerate(depth_tief):
    print(i,d) 
"""
              
#Load oxygen timeseries
#------------------------- 
data = np.load("/home/ole/windows/Preprocessing_TC_stations/PME/data/PME_emb217_TC_Flach_oxygen.npz",allow_pickle = True)
#print(data.files)
#concentration = data["concentration"]
saturation_flach = data["saturation"]
label_list_flach = data["label_list"]
PME_depth_flach = 78-label_list_flach[:,1].astype("float")
PME_utc_flach = data["utc"]
#print(depth_flach)
#print(np.shape(saturation_flach))
    
data = np.load("/home/ole/windows/Preprocessing_TC_stations/PME/data/PME_emb217_TC_Tief_oxygen.npz",allow_pickle=True)
#print(data.files)
#concentration = data["concentration"]
saturation_tief = data["saturation"]
label_list_tief = data["label_list"]
PME_depth_tief = 99-label_list_tief[:,1].astype("float")
PME_utc_tief = data["utc"]
#print(PME_depth_tief)  
#print(np.shape(saturation_tief))




#Load bathymetrie data
#-------------------------   
from PIL import Image
im = Image.open('/home/ole/windows/Area_of_interest.tif')
#im.show()
bathymetrie_array = np.array(im)
ny, nx = bathymetrie_array.shape


NCOLS = 320
NROWS = 103
XLLCORNER = 20.41776252494949
YLLCORNER = 57.30439748660123
CELLSIZE = 0.0010416666699999966
NODATA_VALUE = -32767.0
XRCORNER = XLLCORNER+nx*CELLSIZE
YRCORNER = YLLCORNER+ny*CELLSIZE

lons = np.arange(XLLCORNER,XRCORNER,CELLSIZE) 
lats  = np.arange(YLLCORNER,YRCORNER,CELLSIZE)  

lons = np.linspace(XLLCORNER,XLLCORNER+nx*CELLSIZE,nx) 
lats  = np.linspace(YLLCORNER,YLLCORNER+ny*CELLSIZE,ny)

emb217_flach = [57.3200,20.6150]
print("emb217_flach",emb217_flach)

emb217_tief = [57.3200,20.600]
print("emb217_tief",emb217_tief)

#for time_index in range(utc.size):


#Set up the picture:
output_picture, (curr_axis, oxygen_axis) = plt.subplots(nrows = 2, ncols = 1, gridspec_kw={'height_ratios': [4, 1]})





#figure 1
#set the xticks (dateformat and tick location)
hfmt = mdates.DateFormatter('%d %b')

"""
axarr1[row].xaxis.set_major_locator(mdates.DayLocator())
axarr1[row].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
axarr1[row].xaxis.set_major_formatter(hfmt)
"""

map_ax = Basemap(ax=curr_axis, llcrnrlon= XLLCORNER, llcrnrlat= YLLCORNER, urcrnrlon = XRCORNER, urcrnrlat = YRCORNER, suppress_ticks=False)
#map_ax.drawcoastlines()

#xin = np.linspace(map_ax.xmin,map_ax.xmax,nx) # nx is the number of x points on the grid
#yin = np.linspace(map_ax.ymin,map_ax.ymax,ny) # ny in the number of y points on the grid

lons, lats = np.meshgrid(lons,lats)

xx, yy = map_ax(lons,lats)

print(utc_flach[0],utc_flach[-1])
print(utc_tief[0],utc_tief[-1])

print(np.max(bathymetrie_array),np.min(bathymetrie_array))
levels = np.arange(50,110,5)
cmap_RdBu = plt.get_cmap('RdBu_r')
cntrf = map_ax.contourf(xx, yy,bathymetrie_array, levels,cmap="Blues", latlon = True)
cntr = map_ax.contour(xx, yy,bathymetrie_array, levels, colors ="black", latlon = True)

cbar = map_ax.colorbar(cntrf,location='right',pad='5%')
cbar.ax.invert_yaxis() 
cbar.set_label("depth [m]")

map_ax.plot(emb217_flach[1],emb217_flach[0],".",color = "red", latlon=True)
map_ax.plot(emb217_tief[1],emb217_tief[0],".", color = "green", latlon=True)


curr_axis.set_xlim(20.58,20.65)
curr_axis.set_ylim(YLLCORNER,57.34)

curr_axis.set_xlabel("longitude")
curr_axis.set_ylabel("latitude")



oxygen_axis.plot(PME_utc_flach,saturation_flach[0,:],color = "tab:red")
oxygen_axis.plot(PME_utc_tief,saturation_tief[0,:],color = "tab:green")

oxygen_axis.set_xlabel("2019")
oxygen_axis.set_ylabel("oxygen")
oxygen_axis.xaxis.set_major_locator(mdates.DayLocator())
oxygen_axis.xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
oxygen_axis.xaxis.set_major_formatter(hfmt)

output_picture.subplots_adjust(top=0.939,bottom=0.081,left=0.08,right=0.939,hspace=0.196,wspace=0.12)

#axarr3[0].set_title(title_fig1,fontsize=16)
#axarr3[1].set_title(title_fig2,fontsize=16) 

output_picture.suptitle("emb217",fontsize=20,fontweight="bold")

output_picture.set_size_inches(1.618*7.2,7.2)#(18,10.5)
    

# initialization function: plot the background of each frame
def init():
    scale = 0.1
    curr_axis.arrow(x=emb217_flach[1],y=emb217_flach[0], dx=scale*west_east_flach[0], dy = scale*north_south_flach[0], color = "red", alpha = 1)
    curr_axis.arrow(x=emb217_tief[1],y=emb217_tief[0], dx=scale*west_east_tief[0], dy = scale*north_south_tief[0], color = "green", alpha = 1)


# animation function.  This is called sequentially
def animate(i):
    scale = 0.1
    curr_axis.arrow(x=emb217_flach[1],y=emb217_flach[0], dx=scale*west_east_flach[0], dy = scale*north_south_flach[0], color = "red", alpha = 1)
    curr_axis.arrow(x=emb217_tief[1],y=emb217_tief[0], dx=scale*west_east_tief[0], dy = scale*north_south_tief[0], color = "green", alpha = 1)
    oxygen_axis.axvline(mdates.date2num(start_flach_2017))


plt.show()


