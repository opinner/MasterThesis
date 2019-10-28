#---------------------------------------------------------------------------#
#Fehler berechnen
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



#Constants
#---------------------------------------------------------------------------


#Set which fraction of missing data per depth is acceptable
acceptance_limit = 0.2

#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/adcp/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/adcp/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Lars_TC/ADCP/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/adcp/data"]

LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data"]

for FOLDERNAME in LIST_OF_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    cruisename = splitted_foldername[5]
    flach_or_tief = splitted_foldername[8]
    
    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".mat":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    
    #create ouput pictures, which will be filled later in the code
    
    #figure 1 for the measurements
    f1, axarr1 = plt.subplots(1)

    #figure 2 for the test
    f2, axarr2 = plt.subplots(2, sharex=True)#, sharey = True)
    
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        print("--------------------------------------\nfile:")
        print(datafile_path[25:])
        
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

        

        """
        #dfm = difference from mean 
        try:
            dfm_west_east = west_east[:,:]-np.reshape(np.nanmean(west_east[:,:],axis=1),(number_of_depth_bins,1))
            dfm_north_south = north_south[:,:]-np.reshape(np.nanmean(north_south[:,:],axis=1),(number_of_depth_bins,1))
        except RuntimeWarning:
            print("depth with no data")
        #TODO Why has this no effect? How to catch the warning
        #https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing  

        print("Interpolation at following depths:")
        print("row index\t depth\t amount of NaN\t ")
        for i in np.arange(number_of_depth_bins):
            Number_of_Nans_west_east = np.count_nonzero(np.isnan(dfm_west_east[i,:]))  
            Number_of_Nans_north_south = np.count_nonzero(np.isnan(dfm_north_south[i,:])) 

            

            #Only accepted if missing data is less than the set fraction "acceptance limit" of the total data
            Is_amount_of_NaNs_in_west_east_tolerable = Number_of_Nans_west_east < acceptance_limit*dfm_west_east[i,:].size
            Is_amount_of_NaNs_in_north_south_tolerable = Number_of_Nans_north_south < acceptance_limit*dfm_north_south[i,:].size

            if (Is_amount_of_NaNs_in_west_east_tolerable and Number_of_Nans_west_east != 0) :
            
                print("\t",i,"\t",np.around(depth[i],3),"\t",Number_of_Nans_west_east) #Depths should be the same for west_east and north_south
            
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
        """
        
          
        # Dissipation rate calculation
        #-------------------------------------------------------------------------------------------------
        averaging_period_in_minutes = 16 #in minutes
        
        print((utc[-1]-utc[0]))
        print(int((utc[-1]-utc[0]).total_seconds()/60))
        
        total_measurement_time_in_minutes = int((utc[-1]-utc[0]).total_seconds()/60)
        number_of_subarrays = round(total_measurement_time_in_minutes/averaging_period_in_minutes)
        
        print("Number",number_of_subarrays)
        
        #np.split(vertical_v,number_of_subarrays)  #has to split into sub-arrays of equal length
        
        print("shape:",np.shape(vertical_v))
        sub_arrays = np.array_split(vertical_v,number_of_subarrays,axis = 1)  #creates "almost equal" sub-arrays
        sub_time_arrays = np.array_split(utc,number_of_subarrays) 
        print("length:",len(sub_arrays))
        
        #for i in np.arange(number_of_depth_bins):
        #    print(i,depth[i])
        

        """
        #all depths at the same time
        #-----------------------------------------------------------------------
        
        
        #loop over all subarrays
        for sub_array in sub_arrays:
        
            summation = np.zeros(number_of_depth_bins,number_of_depth_bins)
        
            #loop over the timesteps in a sub array
            for i in np.arange(np.shape(sub_array)[1]):
                summation += np.reshape(sub_array[:,i],(-1,1)) - np.reshape(sub_array[:,i],(-1,1)).T    
        
            
            summation/np.shape(sub_array)[1]
                
        
        #for timestep in np.arange(sub_arrays[0]
        #fluctuation_at_z = sub_arrays[0] - np.reshape(np.nanmean(sub_arrays[0],axis = 1),(60,1))
   
        print("fluctuation",np.shape(fluctuation_at_z))
        """

        #more easy, only one depth 
        # TODO Assumes a depth bin Size of 1m
        #-----------------------------------------------------------------------
        current_depth_index = 40
        maximum_distance_r = 12        
        velocity_difference_up = [] #np.zeros((number_of_subarrays,number_of_depth_bins,maximum_distance_r))
        velocity_difference_down = []
        
        #------------------------------------------------------------------------------------------------------------
        current_time_index = 500
        current_depth_index = 40
        sub_array = sub_arrays[current_time_index]
        point_in_time = sub_time_arrays[current_time_index][-1]
        
        mean_velocity = np.nanmean(sub_array[current_depth_index,:])
        summe_up = np.zeros(min(maximum_distance_r,number_of_depth_bins-1-current_depth_index))
        summe_down = np.zeros(min(maximum_distance_r,current_depth_index-1))
    
        def func(r,N,e):
            C_v = 2.1 #from Wiles, Rippeth et al. 2006
            A = C_v**2 * e**(2/3)
            return N+A*r**(2/3)
            
            
        #loop over the timesteps in a sub array (small timesteps)
        for i in np.arange(np.shape(sub_array)[1]):
            
            upwards = np.zeros(min(maximum_distance_r,number_of_depth_bins-1-current_depth_index))-1
            downwards = np.zeros(min(maximum_distance_r,current_depth_index-1))-1
        
            fluctuation_at_z = sub_array[current_depth_index,i] - mean_velocity
            
            #print("upwards shape",np.shape(upwards))
                   
            #TODO time average
            #Loops over the distances r
        
            #velocity difference upwards
            for j in np.arange(0,len(upwards)):
                fluctuation_at_z_plus_r = sub_array[current_depth_index+j+1,i] - mean_velocity
        
                upwards[j] = (fluctuation_at_z - fluctuation_at_z_plus_r)**2
        
            summe_up += upwards  
        
            #velocity difference upwards
            for j in np.arange(0,len(downwards)):
                fluctuation_at_z_plus_r = sub_array[current_depth_index-j-1,i] - mean_velocity
        
                downwards[j] = (fluctuation_at_z - fluctuation_at_z_plus_r)**2
        
            summe_down += downwards
            
    
        structure_function_up = np.asarray(summe_up/np.shape(sub_array)[1])
        
        structure_function_down = np.asarray(summe_down/np.shape(sub_array)[1])
    
        #Dissipation up
        #array [1,..,12] for a maximum distance of 12
        distance_up = np.arange(1,min(maximum_distance_r+1,number_of_depth_bins-current_depth_index))
         
        mask = ~np.isnan(structure_function_up)
        adjusted_distance_up = distance_up[mask]
        adjusted_difference_up = structure_function_up[mask]
        number_of_fitpoints_up = adjusted_difference_up.size
        
        try:    
            popt_up, pcov = curve_fit(func, adjusted_distance_up, adjusted_difference_up)
            rate_of_dissipation_up = popt_up[1]
        except TypeError:
            rate_of_dissipation_up = np.nan
        
        
        
        #Dissipation down
        #array [1,..,12] for a maximum distance of 12
        distance_down = np.arange(1,min(maximum_distance_r+1,current_depth_index))
         
        #print(distance, np.shape(distance))
        mask = ~np.isnan(structure_function_down)
        adjusted_distance_down = distance_down[mask]
        adjusted_difference_down = structure_function_down[mask]
        number_of_fitpoints_down = adjusted_difference_down.size
        
        
        try:    
            popt_down, pcov_down = curve_fit(func, adjusted_distance_down, adjusted_difference_down)
            rate_of_dissipation_down = popt_down[1]
        except TypeError:
            rate_of_dissipation_down = np.nan
        

        
        #weighted average
        total_number_of_points = number_of_fitpoints_up + number_of_fitpoints_down
        if total_number_of_points != 0:
            a = number_of_fitpoints_up/total_number_of_points
            b = number_of_fitpoints_down/total_number_of_points
            avrg_dissipation = a * rate_of_dissipation_up + b * rate_of_dissipation_down 
        else:
            avrg_dissipation = np.nan
        

        print("dissipation up:",rate_of_dissipation_up)
        print("dissipation down:",rate_of_dissipation_down)
        print("averaged dissipation:",avrg_dissipation)
        print("at depth",depth[current_depth_index])
        print("at time",sub_time_arrays[current_time_index][0],"\n")
        
        
 
        
        #Plot
        #---------------------------------------------------------------
        axarr1.plot(distance_up,structure_function_up,"r")
        axarr1.plot(distance_down,structure_function_down,"b")
        axarr1.plot(adjusted_distance_up, func(adjusted_distance_up, *popt_up), 'r--',label='fit: N=%5.3f, e=%5.3f' % tuple(popt_up))
        axarr1.plot(adjusted_distance_down, func(adjusted_distance_down, *popt_down), 'b--',label='fit: N=%5.3f, e=%5.3f' % tuple(popt_down)) 
        
        #axarr2[1].plot(distance,velocity_difference_up[600,:]*10**3)
        #axarr2[1].plot(distance,velocity_difference_down[600,:]*10**3)
        #axarr2[1].plot(distance,velocity_difference_up[500,:])

        """ 
        #Crosscheck 
        check_array = sub_arrays[500] 
        z_array = check_array[current_depth_index]
        r_array = check_array[current_depth_index+7]
        
        
        cross_check = np.mean((z_array - r_array)**2)
        
        print("check array",np.shape(check_array))
        print("chrosscheck",cross_check == velocity_difference_up[500,6])
        print(velocity_difference_up[500,6])
        print(cross_check)
                  
        #https://en.wikipedia.org/wiki/Bilinear_interpolation
        #https://stackoverflow.com/questions/32800623/how-to-get-the-fft-of-a-numpy-array-to-work
        #https://docs.scipy.org/doc/numpy/reference/generated/numpy.interp.html   
        # Is interpolation necessary? Data is equidistant, so a mean is sufficient?   
        """

        #set plot limit to either zero or the minimum depth
        bottom_limit = max(0,min(depth))
    

        #Einstellungen Plot:
        vmin = -0.3
        vmax = +0.3     

        #fill figure 1 with data
        #img1_1 = axarr1[0].pcolormesh(utc,depth,vertical_v, vmin = -0.0075, vmax = 0.0075, cmap = plt.cm.RdYlBu_r)
        #img1_2 = axarr1[1].pcolormesh(utc,depth,vertical_v, vmin = -0.0075, vmax = 0.0075, cmap = plt.cm.RdYlBu_r)

    
    #Preferences of the plots

    #figure 1
    #set the xticks (dateformat and tick location)
    """
    hfmt = mdates.DateFormatter('%d %b')
    axarr1[1].xaxis.set_major_locator(mdates.DayLocator())
    axarr1[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
    axarr1[1].xaxis.set_major_formatter(hfmt)
    
    axarr1[1].set_ylim(bottom = bottom_limit)
    axarr1[1].invert_yaxis()
    """
    
    title_fig1_1 = "adcp "+cruisename+" "+flach_or_tief
    title_fig1_2 = "adcp "+cruisename+" "+flach_or_tief+" north component"

    axarr1.set_title(title_fig1_1)

    #colorbar(img1_1).set_label('velocity [m/s]')
    #colorbar(img1_2).set_label('velocity [m/s]')

    f1.set_size_inches(12,7)
    
    #Save the plot as png
    plot1_name = "./pictures/"+"adcp_"+cruisename+"_"+flach_or_tief 
    #f1.savefig(plot1_name)
    
    
plt.show()
    
    
