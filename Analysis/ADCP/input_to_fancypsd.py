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

#defines which fraction of missing data is acceptable
acceptance_limit = 0.2

#for the UBUNTU Laptop
#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/adcp/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/adcp/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP1200/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/adcp/data"]

#for the VM
LIST_OF_FOLDERS = ["/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_flach/adcp/data","/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_tief/adcp/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-flach/adcp/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-tief/adcp/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/ADCP1200/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/adcp/data"]
#Handpicked from the complete ADCP plots
desired_depths_emb169 = [20, 50, 60]
desired_depths_emb177 = [20, 50, 60]
desired_depths_emb217 = [20, 50, 60]
    
for FOLDERNAME in LIST_OF_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    cruisename = splitted_foldername[5]
    flach_or_tief = splitted_foldername[8]
    
    #TODO Explain why?
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
            horizontal_cut = desired_depths_emb169
            
        if cruisename == "emb177":
            horizontal_cut = desired_depths_emb177
            
        if cruisename == "emb217":
            horizontal_cut = desired_depths_emb217
        
        print("\n",datafile_path)
        print("from ",depth[0]," to ",depth[-1])
        for desired_depth in horizontal_cut:
            index = np.argmin(np.absolute(depth-desired_depth))
            if (abs(depth[index]-desired_depth)) < 5:

                
                print(cruisename,flach_or_tief,depth[index],desired_depth)

                #Data with the mean subtracted and holes filled
                dfm_north_south[i,:]
            
                #original measurement data
                north_south[i,:]
        
        #sio.savemat(FILENAME,
        
        
        
        #create ouput pictures, which will be filled later in the code
        #figure 1 for overview
        #f1, axarr1 = plt.subplots(nrows = 1, ncols = 1)

        #figure 2 for velocties and spectra
        #f2, axarr2 = plt.subplots(nrows = 3, ncols = 2, sharex = "row", sharey=True)

        

    
