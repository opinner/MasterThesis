#---------------------------------------------------------------------------#
#Creates a multidem#

#TODO Change projection for scale display
#---------------------------------------------------------------------------#
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit 
from mpl_toolkits.basemap import Basemap
from matplotlib import animation
from matplotlib.gridspec import GridSpec
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

#parameters
offset = 3500 #4350#
max_frames = 4000 #3000
displayed_depth = 75

  
##########################################
# Load ADCD data
#########################################    

#Load TC-flach
#-------------------------    
#print(sio.whosmat(FILENAME))
if displayed_depth > 68:
    datafile_path = "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP1200/data/EMB217_TC-flach_adcp1200_val.mat"
else:    
    datafile_path = "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data/EMB217_TC-flach_adcp600_val.mat"
data = sio.loadmat(datafile_path)
data = data["adcpavg"]
substructure = data.dtype

depth_flach = (data["depth"][0][0]).flatten()
number_of_depth_bins_flach = depth_flach.size
assert(np.mean(depth_flach)>0)

rtc = data["rtc"][0][0].flatten()
curr = data["curr"][0][0]

trim_emb217_flach = [0,8086]
rtc = rtc[0:8086]
curr = curr[0:8086,:]
 
print(np.shape(curr))   
all_west_east_flach = np.real(curr).T
all_north_south_flach = np.imag(curr).T

index_adcp_flach = np.argmin(np.abs(depth_flach-displayed_depth))
print("index_adcp_flach",index_adcp_flach,depth_flach[index_adcp_flach])
#interpolate to fill in missing data (possible due to a constant timestep)
#index_adcp_flach = 29
temp_x = all_west_east_flach[index_adcp_flach,:]
temp_y = all_north_south_flach[index_adcp_flach,:]
xi = np.arange(len(temp_x))
mask = np.isfinite(temp_x)
west_east_flach = np.interp(xi, xi[mask], temp_x[mask])
north_south_flach = np.interp(xi, xi[mask], temp_y[mask])

#convert matlab time to utc
utc_flach = np.asarray(mdates.num2date(rtc-366))
print("utc flach:",np.shape(utc_flach),np.shape(west_east_flach))



#Load TC-tief
#-------------------------    
#print(sio.whosmat(FILENAME))
datafile_path = "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/adcp/data/EMB217_TC-tief_adcp300_val.mat"
data = sio.loadmat(datafile_path)
data = data["adcpavg"]
substructure = data.dtype

depth_tief = (data["depth"][0][0]).flatten()
number_of_depth_bins_tief = depth_tief.size
assert(np.mean(depth_tief)>0)

rtc = data["rtc"][0][0].flatten()
curr = data["curr"][0][0]

trim_emb217_tief = [100,8165]
rtc = rtc[100:8165]
    
all_west_east_tief = np.real(curr).T
all_north_south_tief = np.imag(curr).T

all_north_south_tief = all_north_south_tief[:,100:8165]
all_west_east_tief = all_west_east_tief[:,100:8165] 


index_adcp_tief = np.argmin(np.abs(depth_tief-displayed_depth))
print("index_adcp_tief",index_adcp_tief, depth_tief[index_adcp_tief])

#interpolate to fill in missing data (possible due to a constant timestep)
#index_adcp_tief = 7
#print(np.shape(all_west_east_tief))
temp_x = all_west_east_tief[index_adcp_tief,:]
temp_y = all_north_south_tief[index_adcp_tief,:]
xi = np.arange(len(temp_x))
mask = np.isfinite(temp_x)
west_east_tief = np.interp(xi, xi[mask], temp_x[mask])
north_south_tief = np.interp(xi, xi[mask], temp_y[mask])

#print(np.shape(west_east_flach))
#print(np.shape(west_east_tief))

#convert matlab time to utc
utc_tief = np.asarray(mdates.num2date(rtc-366))




equal_time_index = np.argmin(np.abs(utc_flach-utc_tief[0]))
assert(equal_time_index!=0)
utc_flach = utc_flach[equal_time_index:]
west_east_flach = west_east_flach[equal_time_index:]
north_south_flach = north_south_flach[equal_time_index:]

#print(np.shape(utc_flach),np.shape(utc_tief))

if (utc_tief.size-utc_flach.size) > 0:
    west_east_flach = np.pad(west_east_flach,(0,utc_tief.size - utc_flach.size),'constant', constant_values=np.nan)
    north_south_flach = np.pad(north_south_flach,(0,utc_tief.size - utc_flach.size),'constant', constant_values=np.nan)



print(np.shape(west_east_flach),np.shape(west_east_tief))
print(utc_flach[0],utc_flach[-1])
print(utc_tief[0],utc_tief[-1])


              
#Load oxygen timeseries
#------------------------- 
data = np.load("/home/ole/windows/Preprocessing_TC_stations/PME/data/PME_emb217_TC_Flach_oxygen.npz",allow_pickle = True)
#print(data.files)
#concentration = data["concentration"]
saturation_flach = data["saturation"]
label_list_flach = data["label_list"]  #contains sensor ID and hight above the sea floor
depth_at_TC_flach = 78
PME_depth_flach = depth_at_TC_flach-label_list_flach[:,1].astype("float")
print("PME_depth_flach")
print(PME_depth_flach)

PME_utc_flach = data["utc"]
PME_index_flach = np.argmin(np.abs(PME_depth_flach-displayed_depth))
print(PME_index_flach,PME_depth_flach[PME_index_flach])
#print(np.shape(saturation_flach))
    
data = np.load("/home/ole/windows/Preprocessing_TC_stations/PME/data/PME_emb217_TC_Tief_oxygen.npz",allow_pickle=True)
#print(data.files)
#concentration = data["concentration"]
saturation_tief = data["saturation"]
label_list_tief = data["label_list"] #contains sensor ID and hight above the sea floor
depth_at_TC_tief = 99
PME_depth_tief = depth_at_TC_tief-label_list_tief[:,1].astype("float")
print("PME_depth_tief")
print(PME_depth_tief)

PME_utc_tief = data["utc"]
PME_index_tief = np.argmin(np.abs(PME_depth_tief-displayed_depth))
print(PME_index_tief,PME_depth_tief[PME_index_tief])




#Load bathymetrie data
#-------------------------   
from PIL import Image
im = Image.open('/home/ole/windows/all_data/bathymetry_data/Area_of_interest.tif')
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
#print("emb217_flach",emb217_flach)

emb217_tief = [57.3200,20.600]
#print("emb217_tief",emb217_tief)

##############################################################################################################
##############################################################################################################
##############################################################################################################
#PLOTTING
##############################################################################################################
##############################################################################################################
##############################################################################################################

#Set up the picture:
output_picture, (curr_axis, oxygen_axis) = plt.subplots(nrows = 2, ncols = 1, gridspec_kw={'height_ratios': [4, 1]})

map_ax = Basemap(ax=curr_axis, llcrnrlon= XLLCORNER, llcrnrlat= YLLCORNER, urcrnrlon = XRCORNER, urcrnrlat = YRCORNER, suppress_ticks=False)
#map_ax.drawcoastlines()


lons, lats = np.meshgrid(lons,lats)

xx, yy = map_ax(lons,lats)



#print(np.max(bathymetrie_array),np.min(bathymetrie_array))
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



oxygen_axis.plot(PME_utc_flach,saturation_flach[PME_index_tief,:],color = "tab:red")
oxygen_axis.plot(PME_utc_tief,saturation_tief[PME_index_tief,:],color = "tab:green")
#oxygen_axis.plot(utc_tief,west_east_tief[:],color = "tab:red")
#oxygen_axis.plot(utc_tief,north_south_tief[:],color = "tab:green")

oxygen_axis.set_xlabel("2019")
oxygen_axis.set_ylabel(r"O$_2$ saturation [%]")
oxygen_axis.set_title("Oxygen measurements in "+str(displayed_depth)+" m depth", fontweight="bold")
oxygen_axis.xaxis.set_major_locator(mdates.DayLocator())
oxygen_axis.xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
hfmt = mdates.DateFormatter('%d %b')
oxygen_axis.xaxis.set_major_formatter(hfmt)
oxygen_axis.set_ylim(0,78)
#oxygen_axis.set_ylim(-0.3,+0.3)

output_picture.subplots_adjust(top=0.942,bottom=0.066,left=0.053,right=0.976,hspace=0.239,wspace=0.12)

#axarr3[0].set_title(title_fig1,fontsize=16)
#axarr3[1].set_title(title_fig2,fontsize=16) 

output_picture.suptitle("ADCP measurements from emb217 at "+str(displayed_depth)+" m depth",fontsize=20,fontweight="bold")

output_picture.set_size_inches(16,10.5)


scale = 1

quiver_flach = curr_axis.quiver([emb217_flach[1]],[emb217_flach[0]],[west_east_flach[0]],[north_south_flach[0]], color = "red", alpha = 1, scale = scale)
quiver_tief = curr_axis.quiver([emb217_tief[1]],[emb217_tief[0]],[west_east_tief[0]],[north_south_tief[0]], color = "green", alpha = 1, scale = scale)

current_time, = oxygen_axis.plot([mdates.date2num(utc_tief[offset]),mdates.date2num(utc_tief[offset])],[0,80],"k-",lw=2)

plt.show()

# initialization function: plot the background of each frame
def init():
    current_time.set_data([],[])
    return oxygen_axis,

# animation function.  This is called sequentially
def animate(i,quiver_flach,quiver_tief):

    quiver_flach.set_UVC(west_east_flach[offset+i],north_south_flach[offset+i])
    quiver_tief.set_UVC(west_east_tief[i],north_south_tief[offset+i])
    current_time.set_data([mdates.date2num(utc_tief[offset+i]),mdates.date2num(utc_tief[offset+i])],[0,80])

    
    return oxygen_axis,quiver_flach,quiver_tief,


assert((max_frames+offset) < west_east_flach.size)
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(output_picture, animate, fargs = (quiver_flach,quiver_tief), init_func=init,
                               frames=max_frames, interval=10, blit=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save("oxygen_animation_"+str(displayed_depth)+"m.mp4", fps=30, extra_args=['-vcodec', 'libx264'])

#plt.show()

