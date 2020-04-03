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
    pass    
    #datafile_path = "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data/EMB217_TC-flach_adcp600_val.mat"
data = sio.loadmat(datafile_path)
data = data["adcpavg"]
substructure = data.dtype

depth_flach = (data["depth"][0][0]).flatten()
number_of_depth_bins_flach = depth_flach.size
assert(np.mean(depth_flach)>0)

rtc = data["rtc"][0][0].flatten()
curr = data["curr"][0][0]
vert_v = data["vu"][0][0].T

trim_emb217_flach = [0,8086]
rtc = rtc[0:8086]
curr = curr[0:8086,:]
print("vert_v",np.shape(vert_v))

 
print(np.shape(curr))   
all_west_east_flach = np.real(curr).T
all_north_south_flach = np.imag(curr).T


index_adcp_flach = np.argmin(np.abs(depth_flach-displayed_depth))
index_adcp_flach = 35

for i,depth in enumerate(depth_flach):
    print(i,depth,"data:",not np.all(np.isnan(all_west_east_flach[i,:])))

print("index_adcp_flach",index_adcp_flach,depth_flach[index_adcp_flach])
vert_v_flach = vert_v[index_adcp_flach,0:8086]

print("all_west_east_flach",np.shape(all_west_east_flach))

west_east_flach = all_west_east_flach[index_adcp_flach,:]
north_south_flach = all_north_south_flach[index_adcp_flach,:]
print(west_east_flach)

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
vert_v = data["vu"][0][0].T

trim_emb217_tief = [100,8165]
rtc = rtc[100:8165]
print("vert_v",np.shape(vert_v))

    
all_west_east_tief = np.real(curr).T
all_north_south_tief = np.imag(curr).T

all_north_south_tief = all_north_south_tief[:,100:8165]
all_west_east_tief = all_west_east_tief[:,100:8165] 


index_adcp_tief = np.argmin(np.abs(depth_tief-displayed_depth))
index_adcp_tief = 5

print("index_adcp_tief",index_adcp_tief, depth_tief[index_adcp_tief])
vert_v_tief = vert_v[index_adcp_tief,100:8165]

#for i,depth in enumerate(depth_tief):
#    print(i,depth)
    
#interpolate to fill in missing data (possible due to a constant timestep)
#index_adcp_tief = 7
print("all_west_east_tief",np.shape(all_west_east_tief))
west_east_tief = all_west_east_tief[index_adcp_tief,:]
north_south_tief = all_north_south_tief[index_adcp_tief,:]


#print(np.shape(west_east_flach))
#print(np.shape(west_east_tief))

#convert matlab time to utc
utc_tief = np.asarray(mdates.num2date(rtc-366))





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
#PME_index_flach = np.argmin(np.abs(PME_depth_flach-displayed_depth))
PME_index_flach = 4
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
#PME_index_tief = np.argmin(np.abs(PME_depth_tief-displayed_depth))
PME_index_tief = 4
print(PME_index_tief,PME_depth_tief[PME_index_tief])


##############################################################################################################
##############################################################################################################
##############################################################################################################
#PLOTTING
##############################################################################################################
##############################################################################################################
##############################################################################################################

#Set up the picture:
output_picture, axarr = plt.subplots(nrows = 4, ncols = 1, sharex = True)

axarr[0].set_title("Flach: oxygen", fontweight="bold")
axarr[2].set_title("Tief: oxygen", fontweight="bold")

for i,label in enumerate(label_list_flach[:,1]):    
    axarr[0].plot(PME_utc_flach,saturation_flach[i,:], label ="flach "+label)
for i,label in enumerate(label_list_tief[:,1]):
    axarr[2].plot(PME_utc_tief,saturation_tief[i,:], label = "tief "+label)
        
axarr[3].plot(utc_tief,west_east_tief[:],color = "tab:green", label = "WE")
axarr[3].plot(utc_tief,north_south_tief[:],color = "tab:blue", label = "NS")
axarr[3].plot(utc_tief,vert_v_tief[:],color = "k", label = "up/down")
axarr[1].plot(utc_flach,west_east_flach[:],color = "tab:red", label = "WE")
axarr[1].plot(utc_flach,north_south_flach[:],color = "tab:orange", label = "NS")
axarr[1].plot(utc_flach,vert_v_flach[:],color = "k", label = "up/down")
axarr[1].plot(utc_flach,all_west_east_flach[index_adcp_flach-5,:], color = "k")
axarr[1].plot(utc_flach,all_north_south_flach[index_adcp_flach-5,:], color = "tab:blue")

axarr[0].xaxis.set_major_locator(mdates.DayLocator())
axarr[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
hfmt = mdates.DateFormatter('%d %b')
axarr[0].xaxis.set_major_formatter(hfmt)
#axarr[0].tick_params(which='both', width=2)


for axis in axarr:
    axis.tick_params(which='major', length=7, width = 2)
    axis.tick_params(which='minor', length=4)

axarr[1].set_title("shallow mooring "+str(depth_flach[index_adcp_flach])+" depth", fontweight="bold")
axarr[1].set_ylabel("velocity [m/s]")


axarr[3].set_xlabel("2019")
axarr[3].set_title("deep mooring "+str(depth_tief[index_adcp_tief])+" depth", fontweight="bold")
axarr[3].set_ylabel("velocity [m/s]")

leg_flach = axarr[3].legend(loc = "upper right")
leg_tief = axarr[2].legend(loc = "upper right")
leg_ox1 = axarr[1].legend(loc = "upper right")
leg_ox2 = axarr[0].legend(loc = "upper right")

for line_flach,line_tief,line_ox1, line_ox2 in zip(leg_flach.get_lines(),leg_tief.get_lines(),leg_ox1.get_lines(), leg_ox2.get_lines()):
    line_flach.set_linewidth(2.0)
    line_tief.set_linewidth(2.0)
    line_ox1.set_linewidth(2.0)
    line_ox2.set_linewidth(2.0)

#subplots_adjust(top=0.944,bottom=0.066,left=0.06,right=0.979,hspace=0.166,wspace=0.12)
output_picture.set_size_inches(18,10.5) 
output_picture.tight_layout()


import pickle
pickle.dump(output_picture, open('Oxygen_drop.fig.pickle', 'wb'))

plt.show()
