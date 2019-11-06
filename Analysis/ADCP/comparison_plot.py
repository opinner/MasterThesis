#---------------------------------------------------------------------------#
#https://www.sciencedirect.com/science/article/pii/0198014980900461         #
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


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)




LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/adcp/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/adcp/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP1200/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/adcp/data"]


#create ouput pictures, which will be filled later in the code
#figure 1 for TC_Flach
f1, axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = "row", sharey=True)

#figure 2 for TC_Tief
f2, axarr2 = plt.subplots(nrows = 3, ncols = 2, sharex = "row", sharey=True)

max_time_delta = dt.timedelta(minutes=1) #start value
    
for FOLDERNAME in LIST_OF_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    cruisename = splitted_foldername[5]
    flach_or_tief = splitted_foldername[8]
    
    if len(flach_or_tief) > 8:
        flach_or_tief = flach_or_tief[6:]
    
    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".mat":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        print("--------------------------------------\nfile:")
        print(datafile_path[25:])
        print(cruisename,flach_or_tief[3:])
        
        #print(sio.whosmat(FILENAME))
        data = sio.loadmat(datafile_path)

        
        #print(data.keys())
        data = data["adcpavg"]
        substructure = data.dtype
        #print(substructure)
        
        
        rtc = data["rtc"][0][0].flatten()
        curr = data["curr"][0][0]
        vertical_v = data["vu"][0][0].T
        
        
        #exceptional trim for the recovery process
        if FOLDERNAME == "/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/adcp/data":
            curr = curr[0:5488]
            rtc = rtc[0:5488]
            vertical_v = vertical_v[0:5488]
            
            
        if FOLDERNAME == "/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/adcp/data":
            curr = curr[0:12775]
            rtc = rtc[0:12775]
            vertical_v = vertical_v[0:12775]
            
        west_east = np.real(curr).T
        north_south = np.imag(curr).T

        #convert matlab time to utc
        utc = np.asarray(pl.num2date(rtc-366))

        depth = (data["depth"][0][0]).flatten()
        number_of_depth_bins = depth.size
        assert(np.mean(depth)>0)

        #calculate the timesteps of the measurments
        timestep = np.nanmean(np.diff(rtc))
        print("timestep in minutes",timestep*24*60) #result circa 900 s = 15min
        print("timestep in seconds",timestep*24*60*60)
        #standard deviation from that 
        #print(np.nanstd(np.diff(rtc))*86400)
          


        #set plot limit to either zero or the minimum depth
        bottom_limit = max(0,min(depth))
    

        #Einstellungen Plot:
        vmin = -0.3
        vmax = +0.3     

        if flach_or_tief[3:].lower() == "flach":
            
            if cruisename == "emb169":
                y_position = 0
                start_flach_2017 = utc[0]
            elif cruisename == "emb177":
                y_position = 1
                start_flach_2018 = utc[0]
            elif cruisename == "emb217":
                y_position = 2
                start_flach_2019 = utc[0]
                
                max_time_delta = max(max_time_delta,utc[-1]-utc[-0])
            else:
                print("ERROR: ID der Ausfahrt nicht bestimmbar")
                
                
            #fill figure 1 with data
            img1_1 = axarr1[y_position,0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = plt.cm.RdYlBu_r)
            img1_2 = axarr1[y_position,1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = plt.cm.RdYlBu_r)
            
            #colorbar(img1_1).set_label('velocity [m/s]')
            #colorbar(img1_2).set_label('velocity [m/s]')
            
        
        elif flach_or_tief[3:].lower() == "tief":
            
            if cruisename == "emb169":
                y_position = 0
                start_tief_2017 = utc[0]
            elif cruisename == "emb177":
                y_position = 1
                start_tief_2018 = utc[0]
            elif cruisename == "emb217":
                y_position = 2
                start_tief_2019 = utc[0]    
            else:
                print("ERROR: ID der Ausfahrt nicht bestimmbar")
                
            #fill figure 2 with data
            #fill figure 1 with data
            img2_1 = axarr2[y_position,0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = plt.cm.RdYlBu_r)
            img2_2 = axarr2[y_position,1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = plt.cm.RdYlBu_r)
        
            #colorbar(img2_1).set_label('velocity [m/s]')
            #colorbar(img2_2).set_label('velocity [m/s]')
    
        else:
            print("ERROR: Ort (flach_or_tief) nicht feststellbar")
    
    
#Preferences of the plots

#figure 1
#set the xticks (dateformat and tick location)
hfmt = mdates.DateFormatter('%d %b')
for row in range(3):
    for column in range(2):
        axarr1[row,column].xaxis.set_major_locator(mdates.DayLocator())
        axarr1[row,column].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
        axarr1[row,column].xaxis.set_major_formatter(hfmt)

        axarr1[row,column].set_ylim(bottom = 0, top = 80)
        print("xlim:",row,column,axarr1[row,column].get_xlim())
        axarr1[row,column].invert_yaxis()


        axarr2[row,column].xaxis.set_major_locator(mdates.DayLocator())
        axarr2[row,column].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
        axarr2[row,column].xaxis.set_major_formatter(hfmt)

        
        axarr2[row,column].set_ylim(bottom = 0, top = 90)
        axarr2[row,column].invert_yaxis()





for row in range(3):
    for column in range(2):
        print("xlim:",row,column,axarr1[row,column].get_xlim())

axarr1[0,0].set_xlabel("2017")
axarr1[0,1].set_xlabel("2017")
axarr1[1,0].set_xlabel("2018")
axarr1[1,1].set_xlabel("2018")
axarr1[2,0].set_xlabel("2019")        
axarr1[2,1].set_xlabel("2019")        

axarr2[0,0].set_xlabel("2017")
axarr2[0,1].set_xlabel("2017")
axarr2[1,0].set_xlabel("2018")
axarr2[1,1].set_xlabel("2018")
axarr2[2,0].set_xlabel("2019")        
axarr2[2,1].set_xlabel("2019") 

title_fig1 = "west component u"
title_fig2 = "north component v"

f1.subplots_adjust(right=0.9)
cbar_ax = f1.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
f1.colorbar(img1_2, cax=cbar_ax).set_label("velocity [m/s]") 

f2.subplots_adjust(right=0.9)
cbar_ax = f2.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
f2.colorbar(img2_2, cax=cbar_ax).set_label("velocity [m/s]") 
        
axarr1[0,0].set_title(title_fig1)
axarr1[0,1].set_title(title_fig2)

axarr2[0,0].set_title(title_fig1)
axarr2[0,1].set_title(title_fig2)
     
f1.suptitle("ADCP Comparison Flach")
f2.suptitle("ADCP Comparison Tief")


f1.set_size_inches(18,10.5)
f2.set_size_inches(18,10.5)
    
#Save the plot as png
plot1_name = "./pictures/"+"ADCP_comparison_flach" 
f1.savefig(plot1_name)

plot2_name = "./pictures/"+"ADCP_comparison_tief" 
f2.savefig(plot2_name)
    
#preferences of the x-limit
axarr1[0,0].set_xlim(left = mdates.date2num(start_flach_2017), right = mdates.date2num(start_flach_2017 + max_time_delta))
axarr1[0,1].set_xlim(left = mdates.date2num(start_flach_2017), right = mdates.date2num(start_flach_2017 + max_time_delta))
axarr1[1,0].set_xlim(left = mdates.date2num(start_flach_2018), right = mdates.date2num(start_flach_2018 + max_time_delta))
axarr1[1,1].set_xlim(left = mdates.date2num(start_flach_2018), right = mdates.date2num(start_flach_2018 + max_time_delta))
axarr1[2,0].set_xlim(left = mdates.date2num(start_flach_2019), right = mdates.date2num(start_flach_2019 + max_time_delta))
axarr1[2,1].set_xlim(left = mdates.date2num(start_flach_2019), right = mdates.date2num(start_flach_2019 + max_time_delta))

axarr2[0,0].set_xlim(left = mdates.date2num(start_tief_2017), right = mdates.date2num(start_tief_2017 + max_time_delta))
axarr2[0,1].set_xlim(left = mdates.date2num(start_tief_2017), right = mdates.date2num(start_tief_2017 + max_time_delta))
axarr2[1,0].set_xlim(left = mdates.date2num(start_tief_2018), right = mdates.date2num(start_tief_2018 + max_time_delta))
axarr2[1,1].set_xlim(left = mdates.date2num(start_tief_2018), right = mdates.date2num(start_tief_2018 + max_time_delta))
axarr2[2,0].set_xlim(left = mdates.date2num(start_tief_2019), right = mdates.date2num(start_tief_2019 + max_time_delta))
axarr2[2,1].set_xlim(left = mdates.date2num(start_tief_2019), right = mdates.date2num(start_tief_2019 + max_time_delta))

#Save the changend plot again as png
plot1_name = "./pictures/"+"ADCP_comparison_flach_same_scale" 
f1.savefig(plot1_name)

plot2_name = "./pictures/"+"ADCP_comparison_tief_same_scale" 
f2.savefig(plot2_name)
    
plt.show()
    
 
