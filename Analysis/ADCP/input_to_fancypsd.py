##################################################################
#creates 3 horizonatal slices of the adcp measured velocity data per measuring station
#as an input for a matlab routine for computing the power spectral density

#Warning: This program is heavily hardcoded for emb169,emb177 and emb217
###################################################################
#TODO: Refactor the code?
##################################################################
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


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)

#defines which fraction of missing data is acceptable
acceptance_limit = 0.3
#if before the patching the mean gets subtracted
subtract_mean = False

#Preferences Plot:
vmin = -0.3 #lower y-limit of the velocity slices
vmax = +0.3 #upper y-limit of the velocity slices
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
desired_depths_emb217_flach = [20, 50, 75]

desired_depths_emb169_tief = [20,58,85]
desired_depths_emb177_tief = [40, 70, 85]
desired_depths_emb217_tief = [20, 50, 85]

#trim to avoid edge effects (manually picked) 
trim_emb169_flach = [0,-1]
trim_emb177_flach = [87,12778]
trim_emb217_flach = [0,8086]

trim_emb169_tief = [0,-1]
trim_emb177_tief = [60,5491]
trim_emb217_tief = [100,8165]

mooring_names = ["emb217_flach"]#["emb169_flach","emb169_tief","emb177_flach","emb177_tief","emb217_flach","emb217_tief"]

for mooring in  mooring_names: 

    f1, axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
    f2, axarr2 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True) 
    
    original_cuts_u = [[],[],[],]
    original_cuts_v = [[],[],[],]
    patched_cuts_u = [[],[],[],]
    patched_cuts_v = [[],[],[],]
    depths_of_cuts = [[],[],[],]
    nan_mask = [[],[],[],]
    matlab_time = [[],[],[],]     
     
    print(mooring)
         
    for FOLDERNAME in LIST_OF_FOLDERS:
        path = pathlib.Path(FOLDERNAME)
        DATAFILENAMES = []

        splitted_foldername = FOLDERNAME.split("/")
        cruisename = splitted_foldername[5]
        flach_or_tief = splitted_foldername[8][3:].lower()
            
        #Exception to account for a different naming convention
        if len(flach_or_tief) > 6:
            flach_or_tief = flach_or_tief[6:]
        
        current_mooring = cruisename+"_"+flach_or_tief
        #print(current_mooring,mooring)
        if current_mooring != mooring:
            continue
            
        
        #go through all files of specified folder and select only files ending with .mat
        for p in path.iterdir():
            all_files_name = str(p.parts[-1])

            if all_files_name[-4:] == ".mat":
                DATAFILENAMES.append(str(p.parts[-1]))

        DATAFILENAMES = sorted(DATAFILENAMES) 
        

        for DATAFILENAME in DATAFILENAMES:
            datafile_path = FOLDERNAME+"/"+DATAFILENAME
            
            data = sio.loadmat(datafile_path)

            data = data["adcpavg"]
            substructure = data.dtype            
            
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
            #print("timestep in minutes",timestep*24*60) #result circa 900 s = 15min
            #print("timestep in seconds",timestep*24*60*60)

            
            if subtract_mean == True:    
            #dfm = difference from mean 
                try:
                    dfm_west_east = west_east[:,:]-np.reshape(np.nanmean(west_east[:,:],axis=1),(depth.size,1))
                    dfm_north_south = north_south[:,:]-np.reshape(np.nanmean(north_south[:,:],axis=1),(depth.size,1))
                except RuntimeWarning:
                    print("depth with no data")
                #TODO Why has this no effect? How to catch the warning
                #https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing  
                
            else:
                dfm_west_east = np.copy(west_east)   
                dfm_north_south = np.copy(north_south)    
                



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
                    cut_points = trim_emb169_tief
                elif  flach_or_tief == "flach":
                    horizontal_cut = desired_depths_emb169_flach   
                    cut_points = trim_emb169_flach                    
            elif cruisename == "emb177":
                if flach_or_tief == "tief":
                    horizontal_cut = desired_depths_emb177_tief       
                    cut_points = trim_emb177_tief                     
                elif  flach_or_tief == "flach":
                    horizontal_cut = desired_depths_emb177_flach
                    cut_points = trim_emb177_flach 
                                                    
            elif cruisename == "emb217":
                if flach_or_tief == "tief" or flach_or_tief == "Tief":
                    horizontal_cut = desired_depths_emb217_tief
                    cut_points = trim_emb217_tief
                elif  flach_or_tief == "flach" or flach_or_tief == "Flach":
                    horizontal_cut = desired_depths_emb217_flach
                    cut_points = trim_emb217_flach
            
            
            img2_1 = axarr2[0].pcolormesh(utc,depth,west_east, vmin = vmin, vmax = vmax, cmap = cmap)
            img2_2 = axarr2[1].pcolormesh(utc,depth,north_south, vmin = vmin, vmax = vmax, cmap = cmap)
            colorbar(img2_1).set_label("velocity [m/s]")
            colorbar(img2_2).set_label("velocity [m/s]")
            
            axarr2[0].xaxis.set_major_locator(mdates.DayLocator())
            axarr2[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
            axarr2[0].xaxis.set_major_formatter(hfmt)

            axarr2[0].set_ylim(bottom = 0, top = 90)
            axarr2[0].invert_yaxis()
        
            axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
            axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
            axarr2[1].xaxis.set_major_formatter(hfmt)

            axarr2[1].set_ylim(bottom = 0, top = 90)
            axarr2[1].invert_yaxis()      
            
            #axarr2[0].set_xlabel("cycles per hour")   
            axarr2[1].set_xlabel("2019")   
            axarr2[0].set_ylabel("pressure [dbar]")  
            axarr2[1].set_ylabel("pressure [dbar]")  
                                       
            print("\n",datafile_path)
            print("from ",depth[0]," to ",depth[-1])
            
            
            for desired_depth in horizontal_cut:
            
                #finds the index of the depth closely to the prior determined desired depth value
                slice_index = np.argmin(np.absolute(depth-desired_depth))
                
                #print("Assert:",np.count_nonzero(np.isnan(dfm_north_south[slice_index,:])) == 0)
                assert(np.count_nonzero(np.isnan(dfm_north_south[slice_index,:])) == 0)
                assert(np.count_nonzero(np.isnan(dfm_west_east[slice_index,:])) == 0)           
                
                #prevents that the nearest edge gets chosen
                if (abs(depth[slice_index]-desired_depth)) < 1:

                    
                    #which of the 3 slices the current is (used for plotting)
                    order_index = np.argmin(np.absolute(np.asarray(horizontal_cut)-desired_depth))
                    
                    #TODO
                    nan_plot_west_east = ma.masked_where(~np.isnan(west_east[slice_index,:]),dfm_west_east[slice_index,:])
                    nan_plot_north_south = ma.masked_where(~np.isnan(west_east[slice_index,:]),dfm_north_south[slice_index,:])
                    nan_plot_utc = ma.masked_where(np.isnan(west_east[slice_index,:]),utc)
                                    
                    cut_mask_array = np.zeros(west_east[slice_index,:].size, dtype=bool)
                    
                    print("\t",depth[slice_index],desired_depth,order_index)
                    #print(cruisename,flach_or_tief,depth[slice_index],desired_depth,order_index)


                    if cruisename == "emb217" and flach_or_tief == "flach" and order_index == 2:
                        temp = cut_points
                        cut_points = [0,-1]
                                                
                    cut_mask_array[cut_points[0]:cut_points[1]] = True
                    cut_plot_west_east  = ma.masked_where(cut_mask_array,dfm_west_east[slice_index,:])
                    cut_plot_north_south = ma.masked_where(cut_mask_array,dfm_north_south[slice_index,:])  
                    cut_plot_utc = ma.masked_where(cut_mask_array,utc)    
                
                
                    matlab_time[order_index] = (rtc[cut_points[0]:cut_points[1]]).flatten()
                    #Data with the mean subtracted and holes filled
                    patched_cuts_u[order_index] = dfm_west_east[slice_index,cut_points[0]:cut_points[1]]
                    patched_cuts_v[order_index] = dfm_north_south[slice_index,cut_points[0]:cut_points[1]]
                    #original measurement data
                    original_cuts_u[order_index] = west_east[slice_index,cut_points[0]:cut_points[1]]
                    original_cuts_v[order_index] = north_south[slice_index,cut_points[0]:cut_points[1]]
                    #actual depth of the slice
                    depths_of_cuts[order_index] = depth[slice_index]
                    #mask where dfm was patched
                    nan_mask[order_index] = (~np.isnan(west_east[slice_index,:]))[cut_points[0]:cut_points[1]]
                    
                                                                 
                    axarr1[order_index,0].plot(utc,dfm_west_east[slice_index,:])
                    axarr1[order_index,1].plot(utc,dfm_north_south[slice_index,:])

                    axarr1[order_index,0].plot(nan_plot_utc,nan_plot_west_east,"r")
                    axarr1[order_index,1].plot(nan_plot_utc,nan_plot_north_south,"r")
                    
                    for i in range(2):
                        if order_index != 2:
                            axarr1[order_index,i].axvline(x = utc[cut_points[0]],color='g')
                            axarr1[order_index,i].axvline(x = utc[cut_points[1]],color='g')
                        axarr1[order_index,i].set_ylim(-0.3,+0.3)
                  
                    axarr2[0].axhline(y=depth[slice_index],color='g')
                    axarr2[1].axhline(y=depth[slice_index],color='g')                        
                
                
                
                    if cruisename == "emb217" and flach_or_tief == "flach" and order_index == 2:
                        cut_points = temp
                            


    print(matlab_time)
    print(np.shape(matlab_time))
    print(np.shape(matlab_time[0]))
    print(type(matlab_time))
    print(type(matlab_time[0]))
    f1.suptitle(mooring)
    #f1.tight_layout()
    f1.set_size_inches(18,10.5)
    f2.suptitle(mooring)
    f2.tight_layout()
    #f2.set_size_inches(18,10.5)
    

    
    #f2.subplots_adjust(top=0.95,bottom=0.056,left=0.05,right=0.961,hspace=0.2,wspace=0.2)
    f2.set_size_inches(1.618*7.2,7.2)
    plot2_name = plot_name = "/home/ole/share-windows/adcp_slices/"+mooring+"_cut_positions"
    f2.savefig(plot2_name)
    
    SAVEFILENAME = "/home/ole/share-windows/adcp_slices/ADCP_horizontal_slice_"+mooring+".mat" 
    matlab_dictionary = {}
    matlab_dictionary["matlab_time"] = matlab_time
    matlab_dictionary["depths_of_cuts"] = depths_of_cuts
    matlab_dictionary["u"] = original_cuts_u
    matlab_dictionary["patched_u"] = patched_cuts_u
    matlab_dictionary["v"] = original_cuts_v
    matlab_dictionary["patched_v"] = patched_cuts_v
    matlab_dictionary["nan_mask"] = nan_mask
    #save data to a .mat to use in fancypsd.mat
    sio.savemat(SAVEFILENAME,matlab_dictionary)     
    #print("!!!!!!!!!!!!!!!!!!!!!!!!") 
    #print("Saved as"+SAVEFILENAME)
    #print("!!!!!!!!!!!!!!!!!!!!!!!!")     
plt.show()        
  

