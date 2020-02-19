##################################################################
#creates 3 horizonatal slices of the adcp measured velocity data per measuring station
#as an input for a matlab routine for computing the power spectral density

#Warning: This program is heavily hardcoded for emb169,emb177 and emb217
###################################################################
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import datetime as dt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit 
import numpy.ma as ma
import pandas as pd

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#TODO is centering necessary?
#TODO mark the interpolated data points
#TODO cut to the right size
#TODO which time formatting does fancy_psd need?
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)

#defines which fraction of missing data is acceptable
acceptance_limit = 0.2

#Einstellungen Plot:
vmin = -0.3
vmax = +0.3  
cmap = plt.get_cmap('seismic')
cmap.set_bad(color = 'lightgrey')#, alpha = 0.5)
hfmt = mdates.DateFormatter('%d %b')
        
#for the UBUNTU Laptop
#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/adcp/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/adcp/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP1200/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/adcp/data"]

#for the VM
LIST_OF_FOLDERS = ["/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_flach/adcp/data","/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_tief/adcp/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-flach/adcp/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-tief/adcp/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP1200/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/adcp/data"]

#Cuts close to the bottom, over and under the halocline
#Handpicked from the complete ADCP plots
desired_depths_emb169_flach = [58,65,72]
desired_depths_emb177_flach = [30, 50, 61]
desired_depths_emb217_flach = [18, 50, 75]

desired_depths_emb169_tief = [20,58,85]
desired_depths_emb177_tief = [40, 70, 85]
desired_depths_emb217_tief = [18, 50, 85]

#cuts to avoid edge effects (manually picked) 
cuts_emb169_flach = [0,-1]
cuts_emb177_flach = [87,12778]
cuts_emb217_flach = [0,8086]

cuts_emb169_tief = [0,-1]
cuts_emb177_tief = [60,5491]
cuts_emb217_tief = [99,8165]

#moorings:
emb169_flach_original_cuts_u = [[],[],[],]
emb169_flach_original_cuts_v = [[],[],[],]
emb169_flach_centered_cuts_u = [[],[],[],]
emb169_flach_centered_cuts_v = [[],[],[],]
emb169_flach_depths_of_cuts = [[],[],[],]
                
emb169_tief_original_cuts_u = [[],[],[],]
emb169_tief_original_cuts_v = [[],[],[],]
emb169_tief_centered_cuts_u = [[],[],[],]
emb169_tief_centered_cuts_v = [[],[],[],]
emb169_tief_depths_of_cuts = [[],[],[],]

emb177_flach_original_cuts_u = [[],[],[],]
emb177_flach_original_cuts_v = [[],[],[],]
emb177_flach_centered_cuts_u = [[],[],[],]
emb177_flach_centered_cuts_v = [[],[],[],]
emb177_flach_depths_of_cuts = [[],[],[],]

emb177_tief_original_cuts_u = [[],[],[],]
emb177_tief_original_cuts_v = [[],[],[],]
emb177_tief_centered_cuts_u = [[],[],[],]
emb177_tief_centered_cuts_v = [[],[],[],]
emb177_tief_depths_of_cuts = [[],[],[],]

emb217_flach_original_cuts_u = [[],[],[],]
emb217_flach_original_cuts_v = [[],[],[],]
emb217_flach_centered_cuts_u = [[],[],[],]
emb217_flach_centered_cuts_v = [[],[],[],]
emb217_flach_depths_of_cuts = [[],[],[],]

emb217_tief_original_cuts_u = [[],[],[],]
emb217_tief_original_cuts_v = [[],[],[],]
emb217_tief_centered_cuts_u = [[],[],[],]
emb217_tief_centered_cuts_v = [[],[],[],]
emb217_tief_depths_of_cuts = [[],[],[],]


#define the pictures
picture_emb169_flach_cuts, emb169_flach_axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
picture_emb169_flach, emb169_flach_axarr2 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)

picture_emb169_tief_cuts, emb169_tief_axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
picture_emb169_tief, emb169_tief_axarr2 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)

picture_emb177_flach_cuts, emb177_flach_axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
picture_emb177_flach, emb177_flach_axarr2 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)

picture_emb177_tief_cuts, emb177_tief_axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
picture_emb177_tief, emb177_tief_axarr2 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)

picture_emb217_flach_cuts, emb217_flach_axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
picture_emb217_flach, emb217_flach_axarr2 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)

picture_emb217_tief_cuts, emb217_tief_axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
picture_emb217_tief, emb217_tief_axarr2 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)        

            
for FOLDERNAME in LIST_OF_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    cruisename = splitted_foldername[5]
    flach_or_tief = splitted_foldername[8][3:]
    
    #Exception to account for a different naming convention
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
        
        #print("--------------------------------------\nfile:")
        #print(datafile_path[25:])
        #print(cruisename,flach_or_tief[3:])
        
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
        matlab_time = rtc
        utc = np.asarray(pl.num2date(rtc-366))

        depth = (data["depth"][0][0]).flatten()
        number_of_depth_bins = depth.size
        assert(np.mean(depth)>0)

        #calculate the timesteps of the measurments
        timestep = np.nanmean(np.diff(rtc))
        #print("timestep in minutes",timestep*24*60) #result circa 900 s = 15min
        #print("timestep in seconds",timestep*24*60*60)

              
        #dfm = difference from mean 
        try:
            dfm_west_east = west_east[:,:]-np.reshape(np.nanmean(west_east[:,:],axis=1),(depth.size,1))
            dfm_north_south = north_south[:,:]-np.reshape(np.nanmean(north_south[:,:],axis=1),(depth.size,1))
        except RuntimeWarning:
            print("depth with no data")
        #TODO Why has this no effect? How to catch the warning
        #https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing  


        #check for missing data and fill them partly through linear interpolation 
        for i in np.arange(depth.size):
            Number_of_Nans_west_east = np.count_nonzero(np.isnan(dfm_west_east[i,:]))  
            Number_of_Nans_north_south = np.count_nonzero(np.isnan(dfm_north_south[i,:])) 

            #Only accepted if missing data is less than the set fraction "acceptance limit" of the total data
            Is_amount_of_NaNs_in_west_east_tolerable = Number_of_Nans_west_east < acceptance_limit*dfm_west_east[i,:].size
            Is_amount_of_NaNs_in_north_south_tolerable = Number_of_Nans_north_south < acceptance_limit*dfm_north_south[i,:].size
                
            if (Is_amount_of_NaNs_in_west_east_tolerable and Number_of_Nans_west_east != 0) :
            
                #if i == desired_depth_index:
                    #print("Interpolation at following depths:")
                    #print("row index\t depth\t amount of NaN\t ")
                    #print("\t",i,"\t",np.around(depth[i],3),"\t",Number_of_Nans_west_east) #Depths are the same for west_east and north_south
            
                #Linear Interpolation of missing data
                x = dfm_west_east[i,:]
                xi = np.arange(len(x)) #this is only possible because of a constant timestep

                mask = np.isfinite(x)
                dfm_west_east[i,:] = np.interp(xi, xi[mask], x[mask])


            if (Is_amount_of_NaNs_in_north_south_tolerable and Number_of_Nans_north_south != 0) :
            
                x = dfm_north_south[i,:]
                xi = np.arange(len(x)) #this is only possible because of a constant timestep

                mask = np.isfinite(x)
                dfm_north_south[i,:] = np.interp(xi, xi[mask], x[mask])


        
        if cruisename == "emb169":
            if flach_or_tief == "tief":
                horizontal_cut = desired_depths_emb169_tief
                
                img1_1 = emb169_tief_axarr2[0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = cmap)
                img1_2 = emb169_tief_axarr2[1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = cmap)
                colorbar(img1_1)
                colorbar(img1_2)
                
                emb169_tief_axarr2[0].xaxis.set_major_locator(mdates.DayLocator())
                emb169_tief_axarr2[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb169_tief_axarr2[0].xaxis.set_major_formatter(hfmt)

                emb169_tief_axarr2[0].set_ylim(bottom = 0, top = 90)
                emb169_tief_axarr2[0].invert_yaxis()
            
                emb169_tief_axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
                emb169_tief_axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb169_tief_axarr2[1].xaxis.set_major_formatter(hfmt)

                emb169_tief_axarr2[1].set_ylim(bottom = 0, top = 90)
                emb169_tief_axarr2[1].invert_yaxis()
                
            elif  flach_or_tief == "flach":
                horizontal_cut = desired_depths_emb169_flach
                
                img1_1 = emb169_flach_axarr2[0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = cmap)
                img1_2 = emb169_flach_axarr2[1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = cmap)
                colorbar(img1_1)
                colorbar(img1_2)
                
                emb169_flach_axarr2[0].xaxis.set_major_locator(mdates.DayLocator())
                emb169_flach_axarr2[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb169_flach_axarr2[0].xaxis.set_major_formatter(hfmt)

                emb169_flach_axarr2[0].set_ylim(bottom = 0, top = 90)
                emb169_flach_axarr2[0].invert_yaxis()
            
                emb169_flach_axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
                emb169_flach_axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb169_flach_axarr2[1].xaxis.set_major_formatter(hfmt)

                emb169_flach_axarr2[1].set_ylim(bottom = 0, top = 90)
                emb169_flach_axarr2[1].invert_yaxis()
                
            else:
                print(flach_or_tief)
                assert(False)
                
        if cruisename == "emb177":
            if flach_or_tief == "tief":
                horizontal_cut = desired_depths_emb177_tief
                
                emb177_tief_axarr2[0].axvline(x = utc[cuts_emb177_tief[0]],color='k')
                emb177_tief_axarr2[0].axvline(x = utc[cuts_emb177_tief[1]],color='k')
                emb177_tief_axarr2[1].axvline(x = utc[cuts_emb177_tief[0]],color='k')
                emb177_tief_axarr2[1].axvline(x = utc[cuts_emb177_tief[1]],color='k')
                emb177_tief_axarr2[0].set_xlim(left = mdates.date2num(utc[cuts_emb177_tief[0]]), right = mdates.date2num(utc[cuts_emb177_tief[1]]))
                emb177_tief_axarr2[1].set_xlim(left = mdates.date2num(utc[cuts_emb177_tief[0]]), right = mdates.date2num(utc[cuts_emb177_tief[1]]))
                
                img1_1 = emb177_tief_axarr2[0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = cmap)
                img1_2 = emb177_tief_axarr2[1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = cmap)
                colorbar(img1_1)
                colorbar(img1_2)
                
                emb177_tief_axarr2[0].xaxis.set_major_locator(mdates.DayLocator())
                emb177_tief_axarr2[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb177_tief_axarr2[0].xaxis.set_major_formatter(hfmt)

                emb177_tief_axarr2[0].set_ylim(bottom = 0, top = 90)
                emb177_tief_axarr2[0].invert_yaxis()
            
                emb177_tief_axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
                emb177_tief_axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb177_tief_axarr2[1].xaxis.set_major_formatter(hfmt)

                emb177_tief_axarr2[1].set_ylim(bottom = 0, top = 90)
                emb177_tief_axarr2[1].invert_yaxis()
                
            elif  flach_or_tief == "flach":
                horizontal_cut = desired_depths_emb177_flach
                
                img1_1 = emb177_flach_axarr2[0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = cmap)
                img1_2 = emb177_flach_axarr2[1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = cmap)
                colorbar(img1_1)
                colorbar(img1_2)
                
                emb177_flach_axarr2[0].xaxis.set_major_locator(mdates.DayLocator())
                emb177_flach_axarr2[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb177_flach_axarr2[0].xaxis.set_major_formatter(hfmt)

                emb177_flach_axarr2[0].set_ylim(bottom = 0, top = 90)
                emb177_flach_axarr2[0].invert_yaxis()
            
                emb177_flach_axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
                emb177_flach_axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb177_flach_axarr2[1].xaxis.set_major_formatter(hfmt)

                emb177_flach_axarr2[1].set_ylim(bottom = 0, top = 90)
                emb177_flach_axarr2[1].invert_yaxis()
                
            else:
                assert(False)
                
                
        if cruisename == "emb217":
        
            if flach_or_tief == "tief" or flach_or_tief == "Tief":
                horizontal_cut = desired_depths_emb217_tief
                
                img1_1 = emb217_tief_axarr2[0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = cmap)
                img1_2 = emb217_tief_axarr2[1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = cmap)
                
                colorbar(img1_1)
                colorbar(img1_2)
                
                emb217_tief_axarr2[0].xaxis.set_major_locator(mdates.DayLocator())
                emb217_tief_axarr2[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb217_tief_axarr2[0].xaxis.set_major_formatter(hfmt)

                emb217_tief_axarr2[0].set_ylim(bottom = 0, top = 90)
                emb217_tief_axarr2[0].invert_yaxis()
            
                emb217_tief_axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
                emb217_tief_axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb217_tief_axarr2[1].xaxis.set_major_formatter(hfmt)

                emb217_tief_axarr2[1].set_ylim(bottom = 0, top = 90)
                emb217_tief_axarr2[1].invert_yaxis()
                
        
            elif  flach_or_tief == "flach" or flach_or_tief == "Flach":
                horizontal_cut = desired_depths_emb217_flach
                
                img1_1 = emb217_flach_axarr2[0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = cmap)
                img1_2 = emb217_flach_axarr2[1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = cmap)
                colorbar(img1_1)
                colorbar(img1_2)
                
                emb217_flach_axarr2[0].xaxis.set_major_locator(mdates.DayLocator())
                emb217_flach_axarr2[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb217_flach_axarr2[0].xaxis.set_major_formatter(hfmt)

                emb217_flach_axarr2[0].set_ylim(bottom = 0, top = 90)
                emb217_flach_axarr2[0].invert_yaxis()
            
                emb217_flach_axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
                emb217_flach_axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
                emb217_flach_axarr2[1].xaxis.set_major_formatter(hfmt)

                emb217_flach_axarr2[1].set_ylim(bottom = 0, top = 90)
                emb217_flach_axarr2[1].invert_yaxis()
                
            else:
                print(flach_or_tief)
                assert(False)

        
        
        
        
        
        print("\n",datafile_path)
        print("from ",depth[0]," to ",depth[-1])
        
        
        for desired_depth in horizontal_cut:
        
            #finds the index of the depth closely to the prior determined desired depth value
            slice_index = np.argmin(np.absolute(depth-desired_depth))
            
            assert(np.count_nonzero(np.isnan(dfm_north_south[slice_index,:])) == 0)
            assert(np.count_nonzero(np.isnan(dfm_west_east[slice_index,:])) == 0)           
            
            #prevents that the nearest edge gets chosen
            if (abs(depth[slice_index]-desired_depth)) < 1:

                
                #which of the 3 slices the current is (used for plotting)
                order_index = np.argmin(np.absolute(np.asarray(horizontal_cut)-desired_depth))
                nan_plot_west_east = ma.masked_where(~np.isnan(west_east[slice_index,:]),dfm_west_east[slice_index,:])
                nan_plot_north_south = ma.masked_where(~np.isnan(west_east[slice_index,:]),dfm_north_south[slice_index,:])
                nan_plot_utc = ma.masked_where(np.isnan(west_east[slice_index,:]),utc)
                                
                cut_mask_array = np.zeros(west_east[slice_index,:].size, dtype=bool)
                
                print(cruisename,flach_or_tief,depth[slice_index],desired_depth,order_index)

                if cruisename == "emb169":
                    if flach_or_tief == "tief":
                    
                        cut_mask_array[cuts_emb169_tief[0]:cuts_emb169_tief[1]] = True
                
                        cut_plot_west_east  = ma.masked_where(cut_mask_array,dfm_west_east[slice_index,:])
                        cut_plot_north_south = ma.masked_where(cut_mask_array,dfm_north_south[slice_index,:])  
                        cut_plot_utc = ma.masked_where(cut_mask_array,utc) 
                    
                    
                                  
                        emb169_tief_matlab_time = matlab_time[cuts_emb169_tief[0]:cuts_emb169_tief[1]]
                        #Data with the mean subtracted and holes filled
                        emb169_tief_centered_cuts_u[order_index] = dfm_west_east[slice_index,cuts_emb169_tief[0]:cuts_emb169_tief[1]]
                        emb169_tief_centered_cuts_v[order_index] = dfm_north_south[slice_index,cuts_emb169_tief[0]:cuts_emb169_tief[1]]
                        #original measurement data
                        emb169_tief_original_cuts_u[order_index] = west_east[slice_index,cuts_emb169_tief[0]:cuts_emb169_tief[1]]
                        emb169_tief_original_cuts_v[order_index] = north_south[slice_index,cuts_emb169_tief[0]:cuts_emb169_tief[1]]
                        #actual depth of the slice               
                        emb169_tief_depths_of_cuts[order_index] = depth[slice_index]
                        
                        emb169_tief_axarr1[order_index,0].plot(utc,dfm_west_east[slice_index,:])                                               
                        emb169_tief_axarr1[order_index,1].plot(utc,dfm_north_south[slice_index,:])
                        
                        emb169_tief_axarr1[order_index,0].plot(nan_plot_utc,nan_plot_west_east,"r")
                        emb169_tief_axarr1[order_index,1].plot(nan_plot_utc,nan_plot_north_south,"r")
                                                
                        emb169_tief_axarr2[0].axhline(y=depth[slice_index],color='g')
                        emb169_tief_axarr2[1].axhline(y=depth[slice_index],color='g')
                                                
                        
                    elif  flach_or_tief == "flach":
                    
                        
                        cut_mask_array[cuts_emb169_flach[0]:cuts_emb169_flach[1]] = True
                
                        cut_plot_west_east  = ma.masked_where(cut_mask_array,dfm_west_east[slice_index,:])
                        cut_plot_north_south = ma.masked_where(cut_mask_array,dfm_north_south[slice_index,:])  
                        cut_plot_utc = ma.masked_where(cut_mask_array,utc) 
                    
                        emb169_flach_matlab_time = matlab_time[cuts_emb169_flach[0]:cuts_emb169_flach[1]]
                        #Data with the mean subtracted and holes filled
                        emb169_flach_centered_cuts_u[order_index] = dfm_west_east[slice_index,:]
                        emb169_flach_centered_cuts_v[order_index] = dfm_north_south[slice_index,:]                        
                        #original measurement data
                        emb169_flach_original_cuts_u[order_index] = west_east[slice_index,:]                        
                        emb169_flach_original_cuts_v[order_index] = north_south[slice_index,:]
                        #actual depth of the slice
                        emb169_flach_depths_of_cuts[order_index] = depth[slice_index]
                        
                        emb169_flach_axarr1[order_index,0].plot(utc,dfm_west_east[slice_index,:])
                        emb169_flach_axarr1[order_index,1].plot(utc,dfm_north_south[slice_index,:])
                        
                        emb169_flach_axarr1[order_index,0].plot(nan_plot_utc,nan_plot_west_east,"r")
                        emb169_flach_axarr1[order_index,1].plot(nan_plot_utc,nan_plot_north_south,"r")
                        
                        
                        emb169_flach_axarr2[0].axhline(y=depth[slice_index],color='g')
                        emb169_flach_axarr2[1].axhline(y=depth[slice_index],color='g')
                        
                    else:
                        print(flach_or_tief)
                        assert(False)
                        
                if cruisename == "emb177":
                    if flach_or_tief == "tief":
                    
                        cut_mask_array[cuts_emb177_tief[0]:cuts_emb177_tief[1]] = True
                
                        cut_plot_west_east  = ma.masked_where(cut_mask_array,dfm_west_east[slice_index,:])
                        cut_plot_north_south = ma.masked_where(cut_mask_array,dfm_north_south[slice_index,:])  
                        cut_plot_utc = ma.masked_where(cut_mask_array,utc)                     
                    
                        emb177_tief_matlab_time = matlab_time[cuts_emb177_tief[0]:cuts_emb177_tief[1]]
                        #Data with the mean subtracted and holes filled
                        emb177_tief_centered_cuts_u[order_index] = dfm_west_east[slice_index,cuts_emb177_tief[0]:cuts_emb177_tief[1]]
                        emb177_tief_centered_cuts_v[order_index] = dfm_north_south[slice_index,cuts_emb177_tief[0]:cuts_emb177_tief[1]]
                        #original measurement data
                        emb177_tief_original_cuts_u[order_index] = west_east[slice_index,cuts_emb177_tief[0]:cuts_emb177_tief[1]]                        
                        emb177_tief_original_cuts_v[order_index] = north_south[slice_index,cuts_emb177_tief[0]:cuts_emb177_tief[1]]
                        #actual depth of the slice
                        emb177_tief_depths_of_cuts[order_index] = depth[slice_index]
                                               
                        emb177_tief_axarr1[order_index,0].plot(utc,dfm_west_east[slice_index,:])
                        emb177_tief_axarr1[order_index,1].plot(utc,dfm_north_south[slice_index,:])
                        
                        emb177_tief_axarr1[order_index,0].plot(nan_plot_utc,nan_plot_west_east,"r")
                        emb177_tief_axarr1[order_index,1].plot(nan_plot_utc,nan_plot_north_south,"r")                        
                        
                        for i in range(2):
                            emb177_tief_axarr1[order_index,i].axvline(x = utc[cuts_emb177_tief[0]],color='g')
                            emb177_tief_axarr1[order_index,i].axvline(x = utc[cuts_emb177_tief[1]],color='g')
                            emb177_tief_axarr1[order_index,i].set_ylim(-0.3,+0.3)
                                                
                        emb177_tief_axarr2[0].axhline(y=depth[slice_index],color='g')
                        emb177_tief_axarr2[1].axhline(y=depth[slice_index],color='g')
                        
                        
                    elif  flach_or_tief == "flach":
                    
                        cut_mask_array[cuts_emb177_flach[0]:cuts_emb177_flach[1]] = True
                        cut_plot_west_east  = ma.masked_where(cut_mask_array,dfm_west_east[slice_index,:])
                        cut_plot_north_south = ma.masked_where(cut_mask_array,dfm_north_south[slice_index,:])  
                        cut_plot_utc = ma.masked_where(cut_mask_array,utc)    
                    
                    
                    
                        emb177_flach_matlab_time = matlab_time[cuts_emb177_flach[0]:cuts_emb177_flach[1]]
                        #Data with the mean subtracted and holes filled
                        emb177_flach_centered_cuts_u[order_index] = dfm_west_east[slice_index,cuts_emb177_flach[0]:cuts_emb177_flach[1]]
                        emb177_flach_centered_cuts_v[order_index] = dfm_north_south[slice_index,cuts_emb177_flach[0]:cuts_emb177_flach[1]]
                        #original measurement data
                        emb177_flach_original_cuts_u[order_index] = west_east[slice_index,cuts_emb177_flach[0]:cuts_emb177_flach[1]]
                        emb177_flach_original_cuts_v[order_index] = north_south[slice_index,cuts_emb177_flach[0]:cuts_emb177_flach[1]]                        
                        #actual depth of the slice
                        emb177_flach_depths_of_cuts[order_index] = depth[slice_index]
                        
                        emb177_flach_axarr1[order_index,0].plot(utc,dfm_west_east[slice_index,:])
                        emb177_flach_axarr1[order_index,1].plot(utc,dfm_north_south[slice_index,:])
                                               
                        emb177_flach_axarr1[order_index,0].plot(nan_plot_utc,nan_plot_west_east,"r")
                        emb177_flach_axarr1[order_index,1].plot(nan_plot_utc,nan_plot_north_south,"r")

                        for i in range(2):
                            emb177_flach_axarr1[order_index,i].axvline(x = utc[cuts_emb177_flach[0]],color='g')
                            emb177_flach_axarr1[order_index,i].axvline(x = utc[cuts_emb177_flach[1]],color='g')
                            emb177_flach_axarr1[order_index,i].set_ylim(-0.3,+0.3)
                                                   
                        emb177_flach_axarr2[0].axhline(y=depth[slice_index],color='g')
                        emb177_flach_axarr2[1].axhline(y=depth[slice_index],color='g')
                        
                    else:
                        assert(False)
                        
                        
                if cruisename == "emb217":
                
                    if Is_amount_of_NaNs_in_west_east_tolerable == False:   
                        print(depth[slice_index])
                        print(Number_of_Nans_west_east)
                        print(acceptance_limit*dfm_west_east[i,:].size)
                
                
                
                
                    if flach_or_tief == "tief" or flach_or_tief == "Tief":
                    
                        #if order_index == 2:
                        #    cuts_emb217_flach = [0,-1]
                    
                        cut_mask_array[cuts_emb217_tief[0]:cuts_emb217_tief[1]] = True
                        cut_plot_west_east  = ma.masked_where(cut_mask_array,dfm_west_east[slice_index,:])
                        cut_plot_north_south = ma.masked_where(cut_mask_array,dfm_north_south[slice_index,:])  
                        cut_plot_utc = ma.masked_where(cut_mask_array,utc)    
                    
                    
                        emb217_tief_matlab_time = matlab_time[cuts_emb217_tief[0]:cuts_emb217_tief[1]]
                        #Data with the mean subtracted and holes filled
                        emb217_tief_centered_cuts_u[order_index] = dfm_west_east[slice_index,cuts_emb217_tief[0]:cuts_emb217_tief[1]]
                        emb217_tief_centered_cuts_v[order_index] = dfm_north_south[slice_index,cuts_emb217_tief[0]:cuts_emb217_tief[1]]
                        #original measurement data
                        emb217_tief_original_cuts_u[order_index] = west_east[slice_index,cuts_emb217_tief[0]:cuts_emb217_tief[1]]
                        emb217_tief_original_cuts_v[order_index] = north_south[slice_index,cuts_emb217_tief[0]:cuts_emb217_tief[1]]
                        #actual depth of the slice
                        emb217_tief_depths_of_cuts[order_index] = depth[slice_index]                                               
                        
                        emb217_tief_axarr1[order_index,0].plot(utc,dfm_west_east[slice_index,:])
                        emb217_tief_axarr1[order_index,1].plot(utc,dfm_north_south[slice_index,:])

                        emb217_tief_axarr1[order_index,0].plot(nan_plot_utc,nan_plot_west_east,"r")
                        emb217_tief_axarr1[order_index,1].plot(nan_plot_utc,nan_plot_north_south,"r")

                        for i in range(2):
                            emb217_tief_axarr1[order_index,i].axvline(x = utc[cuts_emb217_tief[0]],color='g')
                            emb217_tief_axarr1[order_index,i].axvline(x = utc[cuts_emb217_tief[1]],color='g')
                            emb217_tief_axarr1[order_index,i].set_ylim(-0.3,+0.3)
                                             
                        emb217_tief_axarr2[0].axhline(y=depth[slice_index],color='g')
                        emb217_tief_axarr2[1].axhline(y=depth[slice_index],color='g')                         
                        
                    elif  flach_or_tief == "flach" or flach_or_tief == "Flach":
                    
                        if order_index == 2:
                            cuts_emb217_flach = [0,-1]
                            
                        cut_mask_array[cuts_emb217_flach[0]:cuts_emb217_flach[1]] = True
                        cut_plot_west_east  = ma.masked_where(cut_mask_array,dfm_west_east[slice_index,:])
                        cut_plot_north_south = ma.masked_where(cut_mask_array,dfm_north_south[slice_index,:])  
                        cut_plot_utc = ma.masked_where(cut_mask_array,utc)    
                    
                    
                        emb217_flach_matlab_time = matlab_time[cuts_emb217_flach[0]:cuts_emb217_flach[1]]
                        #Data with the mean subtracted and holes filled
                        emb217_flach_centered_cuts_u[order_index] = dfm_west_east[slice_index,cuts_emb217_flach[0]:cuts_emb217_flach[1]]
                        emb217_flach_centered_cuts_v[order_index] = dfm_north_south[slice_index,cuts_emb217_flach[0]:cuts_emb217_flach[1]]
                        #original measurement data
                        emb217_flach_original_cuts_u[order_index] = west_east[slice_index,cuts_emb217_flach[0]:cuts_emb217_flach[1]]
                        emb217_flach_original_cuts_v[order_index] = north_south[slice_index,cuts_emb217_flach[0]:cuts_emb217_flach[1]]
                        #actual depth of the slice
                        emb217_flach_depths_of_cuts[order_index] = depth[slice_index]
                                             
                        emb217_flach_axarr1[order_index,0].plot(utc,dfm_west_east[slice_index,:])
                        emb217_flach_axarr1[order_index,1].plot(utc,dfm_north_south[slice_index,:])

                        emb217_flach_axarr1[order_index,0].plot(nan_plot_utc,nan_plot_west_east,"r")
                        emb217_flach_axarr1[order_index,1].plot(nan_plot_utc,nan_plot_north_south,"r")
                        
                        for i in range(2):
                            if order_index != 2:
                                emb217_flach_axarr1[order_index,i].axvline(x = utc[cuts_emb217_flach[0]],color='g')
                                emb217_flach_axarr1[order_index,i].axvline(x = utc[cuts_emb217_flach[1]],color='g')
                            emb217_flach_axarr1[order_index,i].set_ylim(-0.3,+0.3)
                      
                        emb217_flach_axarr2[0].axhline(y=depth[slice_index],color='g')
                        emb217_flach_axarr2[1].axhline(y=depth[slice_index],color='g')                        
                    
                    
                    
                        if order_index == 2:
                            cuts_emb217_flach = [0,8086]
                                
                    else:
                        print(flach_or_tief)
                        assert(False)



        
        
FILENAME_emb169_flach = "/home/ole/share-windows/adcp_slices/ADCP_horizontal_slice_emb169_flach.mat"
emb169_flach_matlab_dictionary = {}
emb169_flach_matlab_dictionary["matlab_time"] = emb169_flach_matlab_time
emb169_flach_matlab_dictionary["depths_of_cuts"] = emb169_flach_depths_of_cuts
emb169_flach_matlab_dictionary["u"] = emb169_flach_original_cuts_u
emb169_flach_matlab_dictionary["centered_u"] = emb169_flach_centered_cuts_u
emb169_flach_matlab_dictionary["v"] = emb169_flach_original_cuts_v
emb169_flach_matlab_dictionary["centered_v"] = emb169_flach_centered_cuts_v
#save data to a .mat to use in fancypsd.mat
sio.savemat(FILENAME_emb169_flach,emb169_flach_matlab_dictionary)
#picture_emb169_flach_cuts.savefig(plot1_name)        
        
FILENAME_emb169_tief = "/home/ole/share-windows/adcp_slices/ADCP_horizontal_slice_emb169_tief.mat"
emb169_tief_matlab_dictionary = {}
emb169_tief_matlab_dictionary["matlab_time"] = emb169_tief_matlab_time
emb169_tief_matlab_dictionary["depths_of_cuts"] = emb169_tief_depths_of_cuts
emb169_tief_matlab_dictionary["u"] = emb169_tief_original_cuts_u
emb169_tief_matlab_dictionary["centered_u"] = emb169_tief_centered_cuts_u
emb169_tief_matlab_dictionary["v"] = emb169_tief_original_cuts_v
emb169_tief_matlab_dictionary["centered_v"] = emb169_tief_centered_cuts_v
#save data to a .mat to use in fancypsd.mat
sio.savemat(FILENAME_emb169_tief,emb169_tief_matlab_dictionary)


FILENAME_emb177_flach = "/home/ole/share-windows/adcp_slices/ADCP_horizontal_slice_emb177_flach.mat"
emb177_flach_matlab_dictionary = {}
emb177_flach_matlab_dictionary["matlab_time"] = emb177_flach_matlab_time
emb177_flach_matlab_dictionary["depths_of_cuts"] = emb177_flach_depths_of_cuts
emb177_flach_matlab_dictionary["u"] = emb177_tief_original_cuts_u
emb177_flach_matlab_dictionary["centered_u"] = emb177_tief_centered_cuts_u
emb177_flach_matlab_dictionary["v"] = emb177_tief_original_cuts_v
emb177_flach_matlab_dictionary["centered_v"] = emb177_tief_centered_cuts_v
#save data to a .mat to use in fancypsd.mat
sio.savemat(FILENAME_emb177_flach,emb177_flach_matlab_dictionary)


FILENAME_emb177_tief = "/home/ole/share-windows/adcp_slices/ADCP_horizontal_slice_emb177_tief.mat"
emb177_tief_matlab_dictionary = {}
emb177_tief_matlab_dictionary["matlab_time"] = emb177_tief_matlab_time
emb177_tief_matlab_dictionary["depths_of_cuts"] = emb177_tief_depths_of_cuts
emb177_tief_matlab_dictionary["u"] = emb177_tief_original_cuts_u
emb177_tief_matlab_dictionary["centered_u"] = emb177_tief_centered_cuts_u
emb177_tief_matlab_dictionary["v"] = emb177_tief_original_cuts_v
emb177_tief_matlab_dictionary["centered_v"] = emb177_tief_centered_cuts_v
#save data to a .mat to use in fancypsd.mat
sio.savemat(FILENAME_emb177_tief,emb177_tief_matlab_dictionary)


FILENAME_emb217_flach = "/home/ole/share-windows/adcp_slices/ADCP_horizontal_slice_emb217_flach.mat"
emb217_flach_matlab_dictionary = {}
emb217_flach_matlab_dictionary["matlab_time"] = emb217_flach_matlab_time
emb217_flach_matlab_dictionary["depths_of_cuts"] = emb217_flach_depths_of_cuts
emb217_flach_matlab_dictionary["u"] = emb217_flach_original_cuts_u
emb217_flach_matlab_dictionary["centered_u"] = emb217_flach_centered_cuts_u
emb217_flach_matlab_dictionary["v"] = emb217_flach_original_cuts_v
emb217_flach_matlab_dictionary["centered_v"] = emb217_flach_centered_cuts_v
#save data to a .mat to use in fancypsd.mat
sio.savemat(FILENAME_emb217_flach,emb217_flach_matlab_dictionary)


FILENAME_emb217_tief = "/home/ole/share-windows/adcp_slices/ADCP_horizontal_slice_emb217_tief.mat"
emb217_tief_matlab_dictionary = {}
emb217_tief_matlab_dictionary["matlab_time"] = emb217_tief_matlab_time
emb217_tief_matlab_dictionary["depths_of_cuts"] = emb217_tief_depths_of_cuts
emb217_tief_matlab_dictionary["u"] = emb217_tief_original_cuts_u
emb217_tief_matlab_dictionary["centered_u"] = emb217_tief_centered_cuts_u
emb217_tief_matlab_dictionary["v"] = emb217_tief_original_cuts_v
emb217_tief_matlab_dictionary["centered_v"] = emb217_tief_centered_cuts_v
#save data to a .mat to use in fancypsd.mat
sio.savemat(FILENAME_emb217_tief,emb217_tief_matlab_dictionary)        




picture_emb169_flach_cuts.suptitle("emb169 flach")
picture_emb169_flach_cuts.tight_layout() #TODO top=0.934,bottom=0.081,left=0.041,right=0.979,hspace=0.181,wspace=0.06
picture_emb169_flach_cuts.set_size_inches(18,10.5)
picture_emb169_flach.suptitle("emb169 flach")
picture_emb169_flach.tight_layout()
picture_emb169_flach.set_size_inches(18,10.5)
picture_emb169_tief_cuts.suptitle("emb169 tief")
picture_emb169_tief_cuts.tight_layout()
picture_emb169_tief_cuts.set_size_inches(18,10.5)
picture_emb169_tief.suptitle("emb169 tief")
picture_emb169_tief.tight_layout()
picture_emb169_tief.set_size_inches(18,10.5)

picture_emb177_flach_cuts.suptitle("emb177 flach")
picture_emb177_flach_cuts.tight_layout()
picture_emb177_flach_cuts.set_size_inches(18,10.5)
picture_emb177_flach.suptitle("emb177 flach")
picture_emb177_flach.tight_layout()
picture_emb177_flach.set_size_inches(18,10.5)
picture_emb177_tief_cuts.suptitle("emb177 tief")
picture_emb177_tief_cuts.tight_layout()
picture_emb177_tief_cuts.set_size_inches(18,10.5)
picture_emb177_tief.suptitle("emb177 tief")
picture_emb177_tief.tight_layout()
picture_emb177_tief.set_size_inches(18,10.5)

picture_emb217_flach_cuts.suptitle("emb217 flach")
#picture_emb217_flach_cuts.tight_layout()
picture_emb217_flach_cuts.set_size_inches(18,10.5)
picture_emb217_flach.suptitle("emb217 flach")
picture_emb217_flach.tight_layout()
picture_emb217_flach.set_size_inches(18,10.5)
picture_emb217_tief_cuts.suptitle("emb217 tief")
picture_emb217_tief_cuts.tight_layout()
picture_emb217_tief_cuts.set_size_inches(18,10.5)
picture_emb217_tief.suptitle("emb217 tief")
picture_emb217_tief.tight_layout()
picture_emb217_tief.set_size_inches(18,10.5)

plt.show()        
    

