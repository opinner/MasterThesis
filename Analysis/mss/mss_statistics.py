############################################################
#this program loads all profiles from a cruise, removes outliers
#and compute mean and standard deviation of profiles
#in choosable longitude intervals

##############################################################
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
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
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw.conversions as gswc
import gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import scipy.stats as ss
import warnings
warnings.filterwarnings('ignore') #ignoring warnings, specially numpy warnings

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)
    
    
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169","/home/ole/share-windows/processed_mss/emb177"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb217"]
averaging_intervals_borders = [20.55,20.62]
#averaging_intervals_borders = np.linspace(20.48,20.7,9)
number_of_intervals = len(averaging_intervals_borders)+1
first_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)

outlier_factor = 1.5
#TODO change dependent on the cruise?

for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print(cruisename)    
    #2 borders = 3 intervals
    print(averaging_intervals_borders)

    pot_density_statistic = [None] * number_of_intervals
    salinity_statistic = [None] * number_of_intervals
    temperature_statistic = [None] * number_of_intervals
    dissipation_statistic = [None] * number_of_intervals
    oxygen_sat_statistic = [None] * number_of_intervals
    oxygen_flux_statistic = [None] * number_of_intervals
    all_dissipation_statistic = np.asarray([])
        
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".npz":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    #print(DATAFILENAMES)
    
    
    for DATAFILENAME in DATAFILENAMES:
    
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        transect_name = DATAFILENAME[:-4]
    
        #skip the short "S206" transects
        if transect_name[0:4] == "S106":
            print(transect_name,"skipped")
            continue
            
        print("\n",transect_name)
    
        #define the pictures
        #f1, axarr1 = plt.subplots(nrows = 1, ncols = 7, sharey = True)
        #f2, axarr2 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
    
    

        
        
        data = np.load(datafile_path)
        
        try:
            number_of_profiles = data["number_of_profiles"] #
            lat = data["lat"] #Latitude of the profiles
            lon = data["lon"] #Longitude of the profiles
            distance = data["distance"] #distance from the starting profile (monotonically increasing)
            
            bathymetrie = data["bathymetrie"]
            list_of_bathymetrie_indices = data["list_of_bathymetrie_indices"]
            BBL = data["BBL"]
            list_of_BBL_indices = data["list_of_BBL_indices"]
            BBL_range = data["BBL_range"]
            list_of_BBL_range_indices = data["list_of_BBL_range_indices"]
            
            #interp_pressure = data["interp_pressure"]
            #oxygen_grid = data["oxygen_grid"]
            #oxygen_sat_grid = data["oxygen_sat_grid"]
            #salinity_grid = data["salinity_grid"]
            #consv_temperature_grid = data["consv_temperature_grid"]
            #pot_density_grid = data["pot_density_grid"]
            
            eps_pressure = data["eps_pressure"]
            eps_grid = data["eps_grid"]
            corrected_eps_grid = data["corrected_eps_grid"]
            corrected_eps_wiki_grid = data["corrected_eps_wiki_grid"]
            eps_consv_temperature_grid = data["eps_consv_temperature_grid"]
            eps_oxygen_grid = data["eps_oxygen_grid"] 
            eps_oxygen_sat_grid = data["eps_oxygen_sat_grid"]         
            eps_salinity_grid = data["eps_salinity_grid"]
            
            eps_N_squared_grid = data["eps_N_squared_grid"]
            eps_density_grid = data["eps_density_grid"]
            eps_pot_density_grid = data["eps_pot_density_grid"]
            #eps_viscosity_grid = data["eps_viscosity_grid"]
            eps_Reynolds_bouyancy_grid = data["eps_Reynolds_bouyancy_grid"]
            corrected_eps_Reynolds_bouyancy_grid = data["corrected_eps_Reynolds_bouyancy_grid"]
            eps_wiki_Reynolds_bouyancy_grid = data["eps_wiki_Reynolds_bouyancy_grid"]
            corrected_eps_wiki_Reynolds_bouyancy_grid = data["corrected_eps_wiki_Reynolds_bouyancy_grid"]
        
        except KeyError:
            print(transect_name," is skipped, Error during loading data")
            continue
            
            
        """
        number_of_profiles              number of profiles/casts in the transect
        lat                             latitude in degrees (as a float) of the casts
        lon                             longitude in degrees (as a float) of the casts
        distance                        distance in km from the starting point of the transect
        
        interp_pressure                 equidistant 1D pressure array between the highest and the lowest measured pressure value
        oxygen_grid                     oxygen concentration in in micromol per kg as a grid (number_of_profiles x len(interp_pressure))
        salinity_grid                   salinity in g/kg as a grid (number_of_profiles x len(interp_pressure)) 
        consv_temperature_grid          conservative temperature in degrees Celsius as a grid (number_of_profiles x len(interp_pressure))
        density_grid                    density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
        
        eps_pressure                    pressure values to the dissipation rate values (the pressure distance between points is bigger than in interp_pressure) 
        eps_grid                        measured dissipation rate values (number_of_profiles x len(eps_pressure))
        eps_consv_temperature_grid      conservative temperature as a grid (number_of_profiles x len(eps_pressure))
        eps_salinity_grid               absolute salinity in g/kg as a grid (number_of_profiles x len(eps_pressure)) 
        eps_oxygen_grid                 oxygen concentration in micromol per kg as a grid (number_of_profiles x len(eps_pressure))
        eps_N_squared_grid              N^2, the Brunt-Vaisala frequency in 1/s^2 as a grid (number_of_profiles x len(eps_pressure))
        eps_pot_density_grid            potential density in kg/m^3 as a grid (number_of_profiles x len(eps_pressure))
        
        eps_viscosity_grid
        eps_Reynolds_bouyancy_grid
        corrected_eps_Reynolds_bouyancy_grid 
        eps_wiki_Reynolds_bouyancy_grid
        corrected_eps_wiki_Reynolds_bouyancy_grid 
        
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
        """
        
        print("Number of profiles:",number_of_profiles)
        
        print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
        
        eps_N_grid = np.sqrt(eps_N_squared_grid)
        #ozmidov scale
        ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))
        
        #conversion from pressure coordinates to depth
        eps_depth = gswc.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
        #bathymetrie_in_m = gswc.z_from_p(bathymetrie,np.mean(lat))
        
        #eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
        
        #distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetrie_in_m,(-1,1))
        #boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)
            
       
        #compute the 3 oxygen fluxes with different parametrizations
        oxygen_flux_osborn_grid = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_Skif_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmean(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        outlier_count = 0
        
        for profile in range(number_of_profiles):
            for index,border in enumerate(averaging_intervals_borders):
                #print("index",index)
                if lon[profile]  <= border:
                    interval_number = index
                    break
                    
                #longitude of the profile was greater as the last border --> last interval
                elif lon[profile] > averaging_intervals_borders[-1]:
                    interval_number =  len(averaging_intervals_borders)             
                
            #print(np.round(lon[profile],3),interval_number),np.round(np.nanmedian(np.log10(eps_grid[profile,30:-30])),3),np.shape(np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))))
                       

            if np.nanmean(eps_oxygen_sat_grid[profile]) < 0:
                print(cruisename,transect_name,"negative oxygen values")
                continue

            #if the profile contains only nan values, profile is skipped
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,:])):
                print("NaN profile")
                continue
                
            #right now the criterion is only valid for emb217
            if cruisename == "emb217":
            #check for an outlier profile, ergo too high dissipation rates compared with the surrounding
                if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+2*spread_of_profile_medians):      
                    #print("\toutlier")
                    outlier_count += 1
                    continue
           
            if cruisename == "emb177":
                #index for a depth of 50db
                test_index = np.nanargmin(np.abs(eps_pressure-50))
                #print(eps_pressure[test_index],eps_oxygen_sat_grid[profile,test_index])
                
                #test if the saturation at that depth is under a certain level
                if eps_oxygen_sat_grid[profile,test_index] < 50:
                    print("Halocline is too high!")
                    outlier_count += 1
                    continue
                    
            all_dissipation_statistic = np.append(all_dissipation_statistic,eps_grid[profile])
                
            #if the list is empty
            if np.any(oxygen_flux_statistic[interval_number]) == None:
                #fill it with a reshaped profile
                oxygen_flux_statistic[interval_number] = np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))
                salinity_statistic[interval_number] = np.reshape(eps_salinity_grid[profile,:],(1,-1))
                temperature_statistic[interval_number] = np.reshape(eps_consv_temperature_grid[profile,:],(1,-1)) 
                dissipation_statistic[interval_number] = np.reshape(eps_grid[profile,:],(1,-1)) 
                oxygen_sat_statistic[interval_number] = np.reshape(eps_oxygen_sat_grid[profile,:],(1,-1)) 
                pot_density_statistic[interval_number] = np.reshape(eps_pot_density_grid[profile,:],(1,-1)) 
            else:
                #concatenate all further profiles to the ones already in the array
                oxygen_flux_statistic[interval_number] = np.concatenate((oxygen_flux_statistic[interval_number],np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))),axis=0)
                salinity_statistic[interval_number] = np.concatenate((salinity_statistic[interval_number],np.reshape(eps_salinity_grid[profile,:],(1,-1))),axis=0)
                temperature_statistic[interval_number] = np.concatenate((temperature_statistic[interval_number],np.reshape(eps_consv_temperature_grid[profile,:],(1,-1))),axis=0)
                dissipation_statistic[interval_number] = np.concatenate((dissipation_statistic[interval_number],np.reshape(eps_grid[profile,:],(1,-1))),axis=0)
                oxygen_sat_statistic[interval_number] = np.concatenate((oxygen_sat_statistic[interval_number],np.reshape(eps_oxygen_sat_grid[profile,:],(1,-1))),axis=0) 
                pot_density_statistic[interval_number] = np.concatenate((pot_density_statistic[interval_number],np.reshape(eps_pot_density_grid[profile,:],(1,-1))),axis=0)
        print("removed",outlier_count,"profiles as outliers")
        
        
        ###############################################################################################
        #Plotting of the outlier check pictures
        ###############################################################################################

        cmap_hot = plt.get_cmap('hot_r')
        cmap_hot.set_bad(color = 'lightgrey')
        
        f1, axarr1 = plt.subplots(nrows = 2, ncols = 1, sharex = True)
        axarr1[0].plot(lon,np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        axarr1[0].axhline(transect_median,)
        axarr1[0].axhline(transect_median+outlier_factor*spread_of_profile_medians, ls = "--")
        axarr1[0].axhline(transect_median-outlier_factor*spread_of_profile_medians, ls = "--")
        axarr1[0].set_xlabel("longitude")        
        axarr1[0].set_ylabel(r"median(log10($\epsilon$)")  
               
        plotmesh_longitude = np.append(lon,2*lon[-1]-lon[-2])
        img1_1 = axarr1[1].pcolormesh(plotmesh_longitude,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -9.5, cmap = cmap_hot)
        
        #Draw interval borders
        #for border in averaging_intervals_borders:
        #    axarr1[1].axvline(border)
            
        axarr1[1].invert_yaxis()
        axarr1[1].set_ylabel("pressure [dbar]")   
        axarr1[1].set_xlabel("longitude") 
        colorbar(img1_1).set_label(r"log10(dissipation) $[m^2 s^{-3}]$")     
        f1.set_size_inches(9,5)
        f1.tight_layout()
         
        f1.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/outlier_visual_check/"+cruisename+"_"+transect_name+"_outlier", DPI = 300)         
        #plt.show()      
        plt.close(fig = "all")

        
    ###########################################################################################################################
    
    for interval in oxygen_flux_statistic:
        print(np.shape(interval))
            
    #compute mean and std over the saved intervals
    mean_oxygen_flux = [None] * number_of_intervals
    median_oxygen_flux = [None] * number_of_intervals
    upper_percentile_oxygen_flux = [None] * number_of_intervals
    lower_percentile_oxygen_flux = [None] * number_of_intervals
    max_oxygen_flux = [None] * number_of_intervals
    min_oxygen_flux = [None] * number_of_intervals
        
    mean_salinity = [None] * number_of_intervals
    median_salinity = [None] * number_of_intervals
    upper_percentile_salinity = [None] * number_of_intervals
    lower_percentile_salinity = [None] * number_of_intervals
    
    mean_temperature = [None] * number_of_intervals
    median_temperature = [None] * number_of_intervals
    upper_percentile_temperature = [None] * number_of_intervals
    lower_percentile_temperature = [None] * number_of_intervals
    
    log_mean_dissipation = [None] * number_of_intervals
    log_median_dissipation = [None] * number_of_intervals
    
    mean_dissipation = [None] * number_of_intervals
    median_dissipation = [None] * number_of_intervals
    upper_percentile_dissipation = [None] * number_of_intervals
    lower_percentile_dissipation = [None] * number_of_intervals
    second_upper_percentile_dissipation = [None] * number_of_intervals
    second_lower_percentile_dissipation = [None] * number_of_intervals
    
    skew_dissipation = [None] * number_of_intervals
    kurtosis_dissipation = [None] * number_of_intervals
            
    mean_oxygen_sat = [None] * number_of_intervals
    median_oxygen_sat = [None] * number_of_intervals
    upper_percentile_oxygen_sat = [None] * number_of_intervals
    lower_percentile_oxygen_sat = [None] * number_of_intervals
    
    mean_pot_density = [None] * number_of_intervals
    upper_percentile_pot_density = [None] * number_of_intervals
    lower_percentile_pot_density = [None] * number_of_intervals
            
    for index in range(number_of_intervals):
        mean_oxygen_flux[index] = np.nanmean(oxygen_flux_statistic[index],axis=0)
        median_oxygen_flux[index] = np.nanmedian(oxygen_flux_statistic[index],axis=0)
        upper_percentile_oxygen_flux[index] = np.nanpercentile(oxygen_flux_statistic[index], first_percentile, axis = 0)
        lower_percentile_oxygen_flux[index] = np.nanpercentile(oxygen_flux_statistic[index], 100 - first_percentile, axis = 0)
        max_oxygen_flux[index] = np.nanpercentile(oxygen_flux_statistic[index], 97.72, axis = 0)
        min_oxygen_flux[index] = np.nanpercentile(oxygen_flux_statistic[index], 100-97.72, axis = 0)
            
        mean_salinity[index] = np.nanmean(salinity_statistic[index],axis=0)
        median_salinity[index] = np.nanmedian(salinity_statistic[index],axis=0)
        upper_percentile_salinity[index] = np.nanpercentile(salinity_statistic[index], first_percentile, axis = 0)
        lower_percentile_salinity[index] = np.nanpercentile(salinity_statistic[index], 100-first_percentile, axis = 0)
        
        mean_temperature[index] = np.nanmean(temperature_statistic[index],axis=0)
        median_temperature[index] = np.nanmedian(temperature_statistic[index],axis=0)
        upper_percentile_temperature[index] = np.nanpercentile(temperature_statistic[index], first_percentile, axis = 0)
        lower_percentile_temperature[index] = np.nanpercentile(temperature_statistic[index], 100-first_percentile, axis = 0)
        
        log_mean_dissipation[index] = np.nanmean(np.log10(dissipation_statistic[index]),axis=0)
        log_median_dissipation[index] = np.nanmedian(np.log10(dissipation_statistic[index]),axis=0) 
        mean_dissipation[index] = np.log10(np.nanmean(dissipation_statistic[index],axis=0))
        median_dissipation[index] = np.log10(np.nanmedian(dissipation_statistic[index],axis=0)) 
    
        upper_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], first_percentile, axis = 0))
        lower_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], 100-first_percentile, axis = 0))

        second_upper_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], 97.72, axis = 0))
        second_lower_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], 100-97.72, axis = 0))
        
        skew_dissipation[index] = ss.skew(np.log10(dissipation_statistic[index]),axis=0, nan_policy = "omit")
        kurtosis_dissipation[index] = ss.kurtosis(np.log10(dissipation_statistic[index]),axis=0, nan_policy = "omit")
                
        mean_oxygen_sat[index] = np.nanmean(oxygen_sat_statistic[index],axis=0) 
        median_oxygen_sat[index] = np.nanmedian(oxygen_sat_statistic[index],axis=0)
        upper_percentile_oxygen_sat[index] = np.nanpercentile(oxygen_sat_statistic[index], first_percentile, axis = 0)
        lower_percentile_oxygen_sat[index] = np.nanpercentile(oxygen_sat_statistic[index], 100-first_percentile, axis = 0)

        mean_pot_density[index] = np.nanmean(pot_density_statistic[index],axis=0) 
        upper_percentile_pot_density[index] = np.nanpercentile(pot_density_statistic[index], first_percentile, axis = 0)
        lower_percentile_pot_density[index] = np.nanpercentile(pot_density_statistic[index], 100-first_percentile, axis = 0)
                  
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
        
    #TODO f_all,all_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True)       
    f_flux,flux_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_salinity,sal_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_temperature,temp_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_dissipation,dissipation_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_skew,skew_axarr = plt.subplots(nrows = 1, ncols = 2*number_of_intervals, sharey = True) 
    f_oxygen,oxygen_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_density,density_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
                    
    for index in range(number_of_intervals):
    
        if index == 0:
            title = "[-"+str(np.round(averaging_intervals_borders[0],2))+"]"
        elif index == number_of_intervals-1:
            title = "["+str(np.round(averaging_intervals_borders[-1],2))+"-]"
        else:
            title = "["+str(np.round(averaging_intervals_borders[index-1],2))+"-"+str(np.round(averaging_intervals_borders[index],2))+"]"  
            
        flux_axarr[index].plot(mean_oxygen_flux[index],eps_pressure,"k", label = "mean")
        #flux_axarr[index].plot(median_oxygen_flux[index],eps_pressure, "--", color = "tab:blue")
        #flux_axarr[index].fill_betweenx(eps_pressure,mean_oxygen_flux[index]-std_oxygen_flux[index],mean_oxygen_flux[index]+std_oxygen_flux[index], alpha = 0.5)
        flux_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        flux_axarr[index].fill_betweenx(eps_pressure,upper_percentile_oxygen_flux[index],lower_percentile_oxygen_flux[index], alpha = 1, color = "tab:blue", label = "84.13% percentile")
        flux_axarr[index].fill_betweenx(eps_pressure,upper_percentile_oxygen_flux[index],max_oxygen_flux[index], alpha = 0.4, color = "tab:blue",label = "97.72% percentile")
        flux_axarr[index].fill_betweenx(eps_pressure,min_oxygen_flux[index],lower_percentile_oxygen_flux[index], alpha = 0.4, color = "tab:blue")
                        
        sal_axarr[index].plot(mean_salinity[index],eps_pressure, "k", label = "mean")
        #sal_axarr[index].fill_betweenx(eps_pressure,mean_salinity[index]-std_salinity[index],mean_salinity[index]+std_salinity[index], alpha = 0.5)
        sal_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        sal_axarr[index].fill_betweenx(eps_pressure,upper_percentile_salinity[index],lower_percentile_salinity[index], alpha = 0.6, label = "84.13% percentile")
        
        temp_axarr[index].plot(mean_temperature[index],eps_pressure, "k", label = "mean")
        #temp_axarr[index].fill_betweenx(eps_pressure,mean_temperature[index]-std_temperature[index],mean_temperature[index]+std_temperature[index], alpha = 0.5)
        temp_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        temp_axarr[index].fill_betweenx(eps_pressure,upper_percentile_temperature[index],lower_percentile_temperature[index], alpha = 0.6, label = "84.13% percentile")
        
        #dissipation_axarr[index].plot(mean_dissipation[index],eps_pressure, label = "log of means")
        #dissipation_axarr[index].fill_betweenx(eps_pressure,mean_dissipation[index]-std_dissipation[index],mean_dissipation[index]+std_dissipation[index], alpha = 0.5)
        dissipation_axarr[index].plot(log_mean_dissipation[index],eps_pressure, "k", label = "mean of logs")
        dissipation_axarr[index].fill_betweenx(eps_pressure,upper_percentile_dissipation[index],lower_percentile_dissipation[index], alpha = 0.6, label = "84.13% percentile")
        dissipation_axarr[index].fill_betweenx(eps_pressure,second_upper_percentile_dissipation[index],upper_percentile_dissipation[index], alpha = 0.4, color = "tab:blue",label = "97.72% percentile")
        dissipation_axarr[index].fill_betweenx(eps_pressure,second_lower_percentile_dissipation[index],lower_percentile_dissipation[index], alpha = 0.4, color = "tab:blue")
                        

        #dissipation_axarr[index].fill_betweenx(eps_pressure,log_mean_dissipation[index]-log_std_dissipation[index],log_mean_dissipation[index]+log_std_dissipation[index], alpha = 0.4, label = "std")
        dissipation_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
         
        skew_axarr[2*index].plot(skew_dissipation[index],eps_pressure, "k-",label = "skew")
        skew_axarr[2*index].plot(kurtosis_dissipation[index],eps_pressure, "-", color = "tab:green", label = "kurtosis") 
        skew_axarr[2*index+1].fill_betweenx(eps_pressure,upper_percentile_dissipation[index],lower_percentile_dissipation[index], alpha = 0.6, label = "84.13% percentile")
        skew_axarr[2*index+1].fill_betweenx(eps_pressure,second_upper_percentile_dissipation[index],upper_percentile_dissipation[index], alpha = 0.4, color = "tab:blue",label = "97.72% percentile")
        skew_axarr[2*index+1].fill_betweenx(eps_pressure,second_lower_percentile_dissipation[index],lower_percentile_dissipation[index], alpha = 0.4, color = "tab:blue")
                        
        skew_axarr[2*index+1].plot(log_mean_dissipation[index],eps_pressure, "k", label = "mean of logs")
        skew_axarr[2*index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        skew_axarr[2*index+1].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
                
        oxygen_axarr[index].plot(mean_oxygen_sat[index],eps_pressure,"k", label = "mean")
        #oxygen_axarr[index].fill_betweenx(eps_pressure,mean_oxygen_sat[index]-std_oxygen_sat[index],mean_oxygen_sat[index]+std_oxygen_sat[index], alpha = 0.5)
        oxygen_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        oxygen_axarr[index].fill_betweenx(eps_pressure,upper_percentile_oxygen_sat[index],lower_percentile_oxygen_sat[index], alpha = 0.6, label = "84.13% percentile")

        density_axarr[index].plot(mean_pot_density[index],eps_pressure,"k", label = "mean")
        density_axarr[index].set_title(str(np.shape(pot_density_statistic[index])[0])+" profiles\n"+title)
        density_axarr[index].fill_betweenx(eps_pressure,upper_percentile_pot_density[index],lower_percentile_pot_density[index], alpha = 0.6, label = "84.13% percentile")
                        
        flux_axarr[index].set_xlabel(r"FO [mmol/(m$^2$d]")
        oxygen_axarr[index].set_xlabel(r"O$_2$ saturation [$\%$]") 
        sal_axarr[index].set_xlabel("salinity [SA]") 
        temp_axarr[index].set_xlabel("consv_temperature [C]")
        dissipation_axarr[index].set_xlabel(r"log10($\epsilon$) [m$^2$ s$^{-3}$]")    
        skew_axarr[2*index+1].set_xlabel(r"log10($\epsilon$) [m$^2$ s$^{-3}$]")    
        density_axarr[index].set_xlabel(r"potential density $\sigma$ [kg/m$^3$]")    
        
        flux_axarr[index].legend(loc = "upper center")
        oxygen_axarr[index].legend()
        sal_axarr[index].legend()
        temp_axarr[index].legend()
        dissipation_axarr[index].legend()
        skew_axarr[index].legend()
        skew_axarr[index+1].legend()
        
        
        
    import matplotlib.ticker as ticker
    tick_spacing = 10
        
    flux_axarr[0].invert_yaxis()   
    flux_axarr[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
    flux_axarr[0].set_xlim((-10,10))
    if cruisename == "emb177":
        flux_axarr[0].set_xlim((-70,40))   
    flux_axarr[0].set_ylabel("pressure [dbar]")
    f_flux.suptitle("Oxygen flux")
       
    f_salinity.suptitle(cruisename+" Salinity")
    sal_axarr[0].invert_yaxis()
    sal_axarr[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
    sal_axarr[0].set_ylabel("pressure [dbar]")
        
    f_temperature.suptitle(cruisename+" temperature")
    temp_axarr[0].invert_yaxis()
    temp_axarr[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
    temp_axarr[0].set_ylabel("pressure [dbar]")
    
    f_dissipation.suptitle(cruisename+" dissipation")
    #dissipation_axarr[0].set_xlim((-10,-5))
    dissipation_axarr[0].invert_yaxis()
    dissipation_axarr[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
    dissipation_axarr[0].set_ylabel("pressure [dbar]")
    
    f_oxygen.suptitle(cruisename+" Oxygen saturation")
    oxygen_axarr[0].invert_yaxis()
    oxygen_axarr[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
    oxygen_axarr[0].set_ylabel("pressure [dbar]")    

    f_density.suptitle(cruisename+" density")
    density_axarr[0].invert_yaxis()
    density_axarr[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
    density_axarr[0].set_ylabel("pressure [dbar]")    
    
    f_skew.suptitle(cruisename+" Properties of the disspation distribution")
    skew_axarr[0].invert_yaxis()
    skew_axarr[0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
    skew_axarr[0].set_ylabel("pressure [dbar]")
    
    f_flux.set_size_inches(18,10.5)
    f_flux.tight_layout() 
    f_flux.subplots_adjust(top=0.9)
    f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_oxygen_flux",dpi=300)
    
    f_temperature.set_size_inches(18,10.5)
    f_temperature.tight_layout() 
    f_temperature.subplots_adjust(top=0.9)
    f_temperature.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_temperature",dpi=300)
    
    f_dissipation.set_size_inches(18,10.5)
    f_dissipation.tight_layout() 
    f_dissipation.subplots_adjust(top=0.9)
    f_dissipation.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_dissipation",dpi=300)
    
    f_oxygen.set_size_inches(18,10.5)
    f_oxygen.tight_layout() 
    f_oxygen.subplots_adjust(top=0.9)
    f_oxygen.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_oxygen",dpi=300)

    f_salinity.set_size_inches(18,10.5)
    f_salinity.tight_layout() 
    f_salinity.subplots_adjust(top=0.9)
    f_salinity.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_salinity",dpi=300)

    f_density.set_size_inches(18,10.5)
    f_density.tight_layout() 
    f_density.subplots_adjust(top=0.9)
    f_density.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_density",dpi=300)
    
    f_skew.set_size_inches(18,10.5)
    f_skew.tight_layout() 
    f_skew.subplots_adjust(top=0.9)
    f_skew.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_skew_dissipation",dpi=300)
    
    
    """
    color = ["tab:blue","tab:red","tab:green"]
    scatter,scatter_axarr = plt.subplots(nrows = 1, ncols = 2, sharey = True) 
    for index in range(number_of_intervals):
        scatter_axarr[0].plot(mean_oxygen_flux[index],mean_density[index],".", color = color[index])
        scatter_axarr[1].plot(log_mean_dissipation[index],mean_density[index],".", color = color[index])
    scatter_axarr[0].invert_yaxis()
    scatter_axarr[0].set_xlim(-45,20)
    scatter_axarr[0].set_title("oxygen flux")
    scatter_axarr[1].set_title("dissipation")
    scatter.set_size_inches(18,10.5)
    scatter.tight_layout()    
    """
    """    
    hist_dissip,hist_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
    
    #filtered_list = list(filter(None, dissipation_statistic))
    #flat_list = [item for sublist in filtered_list for item in sublist]
    #flattened_array = np.asarray(flat_list).flatten()
    
    #print("profile_count",profile_count)
    print(np.shape(all_dissipation_statistic))
    hist_axarr.hist(np.log10(all_dissipation_statistic),bins = 40,edgecolor='black', linewidth=1.2, log = False)
    hist_axarr.axvline(np.nanmean(np.log10(all_dissipation_statistic)), linewidth=2.4, c = "tab:red", label = "mean of logs")
    hist_axarr.axvline(np.log10(np.nanmean(all_dissipation_statistic)), linewidth=2.4, c = "tab:green", label = "log of means")
    hist_axarr.axvline(np.nanmedian(np.log10(all_dissipation_statistic)), linewidth=2.4, ls = "-", c = "tab:orange", label = "median of logs")
    hist_axarr.axvline(np.log10(np.nanmedian(all_dissipation_statistic)), linewidth=2.4, ls = ":", c = "k", label = "log of median")
    
    hist_axarr.legend()
    
    
    textstr = "Properties of the raw distribution:\nSkew = "+str(np.round(ss.skew(all_dissipation_statistic,nan_policy = "omit"),3))+"\n"+"Excess Kurtosis = "+str(np.round(ss.kurtosis(all_dissipation_statistic, nan_policy = "omit"),3))+"\n\nProperties of the logarithmic distribution:\nSkew = "+str(np.round(ss.skew(np.log10(all_dissipation_statistic),nan_policy = "omit"),3))+"\n"+"Excess Kurtosis = "+str(np.round(ss.kurtosis(np.log10(all_dissipation_statistic), nan_policy = "omit"),3))       
    
    props = dict(boxstyle='square', facecolor = "white")
    hist_axarr.text(0.3, 0.97, textstr, transform=hist_axarr.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
         
    hist_axarr.set_xlabel(r"log10($\epsilon$) $[m^2 s^{-3}]$") 
    hist_axarr.set_ylabel("frequency [#]")
                
    hist_dissip.suptitle(cruisename+": dissipation distribution")
    hist_dissip.set_size_inches(18,10.5)
    hist_dissip.tight_layout() 
    hist_dissip.subplots_adjust(top=0.94)
    hist_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_dissipation_hist",dpi=300)
    """
    plt.show()
    
    
    
    
    
