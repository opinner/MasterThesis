############################################################
#this program loads all profiles from a cruise, removes outliers
#and retrieves the maximum oxyen flux in the lowermost meters 
#of the water column in choosable longitude intervals

#TODO plot mean dissipation per transect
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

import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
    
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217"]#,"/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb169"]#,"/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb217"]

rolling_window_size = 12

flux_through_halocline = True #set to True if the flux trough the Halocline instead of the BBL should be computed 
density_interval = True

#only important for a pressure interval 
upper_bound_halocline_in_db = 57 #67
lower_bound_halocline_in_db = 72 #77

#only important for a density interval
#upper_bound_halocline_as_density = 1007.25 #1006.24 
#lower_bound_halocline_as_density = 1007.75 #1008.90 

height_above_ground = 5 #Size of the averaging interval above ground for the BBL, has no meaning if (flux_through_halocline == True)
maximum_reasonable_flux = 500 #float('Inf') #200 #Fluxes above this value will be discarded
acceptable_slope = 2 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 

flux_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_flux_percentile = 97.72
dissip_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_dissip_percentile = 97.72

number_of_dissipation_subplots = 1 #Decide if both the mean and the median subplots is shown or only the mean

 
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    number_of_fluxes_over_the_threshold = 0
    total_number_of_fluxes = 0
    number_of_zero_flux = 0
    amount_of_missing_values = 0
    total_number_of_chosen_profiles = 0    
    total_number_of_valid_profiles = 0
    
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print(cruisename)
    
    if cruisename == "emb217":
        upper_bound_halocline_as_density = 1006.4 #1005.75
        lower_bound_halocline_as_density = 1008.5 #1006.25
    elif cruisename == "emb177":
        upper_bound_halocline_as_density = 1006.9 #1007.4
        lower_bound_halocline_as_density = 1008.2 #1007.9   
    elif cruisename == "emb169":
        upper_bound_halocline_as_density = 1006.5 
        lower_bound_halocline_as_density = 1008.6     

    dissipation_list = []
    BB_flux_list = []
    Shih_flux_list = []
    Osborn_flux_list = []
    longitude_list = []
    bathymetry_list = []
    bathymetry_longitude_list = []
    interval_list = []
    
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
        
        #something is not correct with this measurement
        if cruisename == "emb169" and transect_name[0:4] == "TS13":
            print(transect_name,"skipped, measurement looks wrong")
            continue
        
        #if cruisename == "emb169" and transect_name[0:4] in ["TS111","TS113","TS113","TS114","TS115",]:
        #    print(transect_name,"skipped, measurement looks wrong")
        #    continue
                        
        print("\n",transect_name)
            
        
        data = np.load(datafile_path)
        
        try:
            number_of_profiles = data["number_of_profiles"] #
            lat = data["lat"] #Latitude of the profiles
            lon = data["lon"] #Longitude of the profiles
            distance = data["distance"] #distance from the starting profile (monotonically increasing)
            
            interp_pressure = data["interp_pressure"]
            oxygen_grid = data["oxygen_grid"]
            salinity_grid = data["salinity_grid"]
            consv_temperature_grid = data["consv_temperature_grid"]
            density_grid = data["density_grid"]
            
            eps_pressure = data["eps_pressure"]
            eps_grid = data["eps_grid"]
            corrected_eps_grid = data["corrected_eps_grid"]
            corrected_eps_wiki_grid = data["corrected_eps_wiki_grid"]
            eps_consv_temperature_grid = data["eps_consv_temperature_grid"]
            eps_oxygen_sat_grid = data["eps_oxygen_sat_grid"]
            eps_oxygen_grid = data["eps_oxygen_grid"] 
            
            eps_N_squared_grid = data["eps_N_squared_grid"]
            eps_density_grid = data["eps_density_grid"]
            eps_pot_density_grid = data["eps_pot_density_grid"]
            #eps_viscosity_grid = data["eps_viscosity_grid"]
            eps_Reynolds_bouyancy_grid = data["eps_Reynolds_bouyancy_grid"]
            corrected_eps_Reynolds_bouyancy_grid = data["corrected_eps_Reynolds_bouyancy_grid"]
            eps_wiki_Reynolds_bouyancy_grid = data["eps_wiki_Reynolds_bouyancy_grid"]
            corrected_eps_wiki_Reynolds_bouyancy_grid = data["corrected_eps_wiki_Reynolds_bouyancy_grid"]
            
            """
            number_of_profiles              number of profiles/casts in the transect
            lat                             latitude in degrees (as a float) of the casts
            lon                             longitude in degrees (as a float) of the casts
            distance                        distance in km from the starting point of the transect
            
            interp_pressure                 equidistant 1D pressure array between the highest and the lowest measured pressure value
            oxygen_grid                     oxygen concentration in in micromol per litre as a grid (number_of_profiles x len(interp_pressure))
            salinity_grid                   salinity in g/kg as a grid (number_of_profiles x len(interp_pressure)) 
            consv_temperature_grid          conservative temperature in degrees Celsius as a grid (number_of_profiles x len(interp_pressure))
            density_grid                    density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
            
            eps_pressure                    pressure values to the dissipation rate values (the pressure distance between points is bigger than in interp_pressure) 
            eps_grid                        measured dissipation rate values (number_of_profiles x len(eps_pressure))
            eps_consv_temperature_grid      conservative temperature as a grid (number_of_profiles x len(eps_pressure))
            eps_oxygen_grid                 oxygen concentration in micromol per litre as a grid (number_of_profiles x len(eps_pressure))
            eps_N_squared_grid              N^2, the Brunt-Vaisala frequency in 1/s^2 as a grid (number_of_profiles x len(eps_pressure))
            eps_density_grid                density in kg/m^3 as a grid (number_of_profiles x len(eps_pressure))
            
            eps_viscosity_grid
            eps_Reynolds_bouyancy_grid
            corrected_eps_Reynolds_bouyancy_grid 
            eps_wiki_Reynolds_bouyancy_grid
            corrected_eps_wiki_Reynolds_bouyancy_grid 
            """
        
        except KeyError:
            print(transect_name," is skipped, Error during loading data")
            continue
            
        print("Number of profiles:",number_of_profiles)
        
        print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
        
        #calculate the idices of the bottom and some meters above that
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid,height_above_ground = height_above_ground)
        """
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
        """
        bathymetrie,list_of_bathymetrie_indices = results[0]
        #BBL,list_of_BBL_indices = results[1] #not needed here
        BBL_range,list_of_BBL_range_indices = results[2]
        
        eps_N_grid = np.sqrt(eps_N_squared_grid)
        #ozmidov scale
        ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))
        
        #conversion from pressure coordinates to depth
        eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
        bathymetrie_in_m = gsw.z_from_p(bathymetrie,np.mean(lat))
        
        eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
        
        distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetrie_in_m,(-1,1))
        boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)

        
        oxygen_flux_osborn_grid = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_Skif_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0
        lon_without_outliers = []
        bathymetrie_without_outliers = []
        
        transect_oxygen_flux_statistic = []
        

        #compare the pressure value of the lowest valid data point per profile th the neighbouring profiles to determine outliers
        list_of_short_profiles = thesis.get_list_of_short_profiles(number_of_profiles,bathymetrie,acceptable_slope)
        
        #hard coded to accout for one short profile, that my condtions doesn't recognize
        #if the boolean array contains at least one true
        if np.any(lon == 20.466176666666666):
            list_of_short_profiles.append(np.argmax(lon == 20.466176666666666))
        
        for profile in range(number_of_profiles):
        
            #if the current profile is too short, skip it
            if profile in list_of_short_profiles:
                print(str(lon[profile])+": short profile")
                continue
                
            
            if np.any(eps_oxygen_sat_grid[profile,:] < -0.01):
                temp = eps_oxygen_sat_grid[profile,:]
                print(cruisename,transect_name," is skipped due to ",len(temp[temp<0])," negative oxygen values")
                if len(temp[temp<0]) < 10:
                    print(temp[temp<0])
                continue
                
            #Determine the averaging intervall
            if flux_through_halocline == True:
                if density_interval == True:
                    #print("\n\n",np.shape(eps_pot_density_grid),np.shape(eps_pot_density_grid[profile]),np.nanmin(eps_pot_density_grid[profile]), np.nanmax(eps_pot_density_grid[profile]))
                    
                    from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-upper_bound_halocline_as_density))     
                    to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-lower_bound_halocline_as_density))
                     
                    #print(from_index,to_index,"\n\n")
                                    
                else: #use pressure intervall
                    from_index =  np.argmin(abs(eps_pressure-upper_bound_halocline_in_db))     
                    to_index = np.argmin(abs(eps_pressure-lower_bound_halocline_in_db))
            
            #use the indices that are {height_above_ground} meters above the bathymetry
            else:     
                from_index = int(list_of_BBL_range_indices[profile]) 
                to_index = int(list_of_bathymetrie_indices[profile])
                
                
            #if the profile contains only nan values, profile is skipped
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,:])): #from_index:to_index
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
            
            total_number_of_valid_profiles+=1
                                    
            #if the water colum portion contains only nan values, save only the bathymetrie then skip it
            #useful if the interval to average over is deeper than the current bathymetrie
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,from_index:to_index])): # or  to_index == list_of_bathymetrie_indices[profile]-1:

                #find the correct position in the sorted list
                for index,value in enumerate(bathymetry_longitude_list):
                    if value > lon[profile]:
                        list_position = index
                        break
                    elif index == len(bathymetry_longitude_list)-1:
                        list_position = index+1
                        break
                
                #if the list is empty                   
                if len(bathymetry_longitude_list) == 0:   
                    bathymetry_list.append(bathymetrie[profile])
                    bathymetry_longitude_list.append(lon[profile])
                    interval_list.append([np.nan,np.nan])
                else:
                    bathymetry_list.insert(list_position,bathymetrie[profile])
                    bathymetry_longitude_list.insert(list_position,lon[profile])
                    interval_list.insert(list_position,[np.nan,np.nan])
                continue
   
                
            #find the correct position in the sorted list
            for index,value in enumerate(bathymetry_longitude_list):
                if value > lon[profile]:
                    list_position = index
                    break
                elif index == len(bathymetry_longitude_list)-1:
                    list_position = len(bathymetry_longitude_list)
                    break
                                
            if len(bathymetry_longitude_list) == 0:   
                bathymetry_list.append(bathymetrie[profile])
                bathymetry_longitude_list.append(lon[profile])
                interval_list.append([eps_pressure[from_index],eps_pressure[to_index]])
            else:
                bathymetry_list.insert(list_position,bathymetrie[profile])
                bathymetry_longitude_list.insert(list_position,lon[profile])
                interval_list.insert(list_position,[eps_pressure[from_index],eps_pressure[to_index]])
                                            
            #find the correct position in the sorted list
            for index,value in enumerate(longitude_list):
                if value > lon[profile]:
                    list_position = index
                    break
                elif index == len(longitude_list)-1:
                    list_position = len(longitude_list)
                    break
            
                 
            if len(longitude_list) == 0:    
                dissipation_list.append(eps_grid[profile,from_index:to_index])
                BB_flux_list.append(oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.append(oxygen_flux_Skif_grid[profile,from_index:to_index])
                #Osborn_flux_list.insert(list_position,)
                longitude_list.append(lon[profile])
            
            
            else:
                
                #Sort the current profile into the list            
                dissipation_list.insert(list_position,eps_grid[profile,from_index:to_index])
                BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.insert(list_position,oxygen_flux_Skif_grid[profile,from_index:to_index])
                #Osborn_flux_list.insert(list_position,)
                longitude_list.insert(list_position,lon[profile])


            #print(longitude_list)        
            assert(np.all(longitude_list == sorted(longitude_list)))

            total_number_of_chosen_profiles+=1

        """
        print("removed",outlier_count,"profiles as outliers")
        print("removed",count_of_short_profiles,"profiles as they did not reach the sea floor")

        ###############################################################################################
        #Plotting of the maximum flux values per transect
        ###############################################################################################
        f1, axarr1 = plt.subplots(nrows = 1, ncols = 1, sharex = True)
        
        bathymetrie_axes = axarr1.twinx()
        bathymetrie_axes.set_ylim((min(bathymetrie_without_outliers)-5,max(bathymetrie_without_outliers)))
        bathymetrie_axes.invert_yaxis()
        bathymetrie_axes.fill_between(lon_without_outliers,bathymetrie_without_outliers, np.ones(len(lon_without_outliers))*max(bathymetrie_without_outliers),color = "lightgrey", alpha = 0.8)
    
        axarr1.set_xlabel("longitude [degree]")    
        axarr1.set_ylabel(r"BB oxygen flux [mmol/(m$^2$*d]")
        bathymetrie_axes.set_ylabel("pressure [dbar]")
        #plot maximum flux
        axarr1.plot(lon_without_outliers,transect_oxygen_flux_statistic[:,1])
        
        #plot minimum flux
        axarr1.plot(lon_without_outliers,transect_oxygen_flux_statistic[:,0])
        axarr1.set_xlabel("longitude")        
                   
        f1.set_size_inches(9,5)
        f1.tight_layout()
        f1.subplots_adjust(top=0.95)
        f1.suptitle("max oxygen flux "+cruisename+" "+transect_name+" "+str(len(lon_without_outliers))+" profiles")
        
        if flux_through_halocline == True:
            f1.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/max_flux_transects/"+cruisename+"_"+transect_name+"_halocline", DPI = 300)   
        else: 
            f1.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/max_flux_transects/"+cruisename+"_"+transect_name+"_BBL", DPI = 300)        
        #plt.show()      
        plt.close(fig = "all")
        """
        
    ###########################################################################################################################
    assert(len(longitude_list) != 0)
    assert(np.all(longitude_list == sorted(longitude_list)))
    assert(np.all(bathymetry_longitude_list == sorted(bathymetry_longitude_list)))        
    interval_list = np.asarray(interval_list)
    
    #compute mean and std over the saved intervals
    mean_flux = [None] * total_number_of_chosen_profiles
    median_flux = [None] * total_number_of_chosen_profiles
    upper_percentile_flux = [None] * total_number_of_chosen_profiles
    lower_percentile_flux = [None] * total_number_of_chosen_profiles
    second_upper_percentile_flux = [None] * total_number_of_chosen_profiles
    second_lower_percentile_flux = [None] * total_number_of_chosen_profiles
  
    #bathymetrie_mean = [None] * number_of_intervals
    #bathymetrie_percentile = [None] * number_of_intervals
    
    log_mean_dissipation = [None] * total_number_of_chosen_profiles
    arith_mean_dissipation = [None] * total_number_of_chosen_profiles
    median_dissipation = [None] * total_number_of_chosen_profiles
    lower_percentile_dissip = [None] * total_number_of_chosen_profiles
    upper_percentile_dissip = [None] * total_number_of_chosen_profiles

    second_lower_percentile_dissip = [None] * total_number_of_chosen_profiles
    second_upper_percentile_dissip = [None] * total_number_of_chosen_profiles

    """
    mean_dissipation_med = [None] * number_of_profiles
    median_dissipation_med = [None] * number_of_profiles
    #std_dissipation_med = [None] * number_of_profiles
    lower_percentile_dissip_med = [None] * number_of_profiles
    upper_percentile_dissip_med = [None] * number_of_profiles

    second_lower_percentile_dissip_med = [None] * number_of_profiles
    second_upper_percentile_dissip_med = [None] * number_of_profiles
    """
    
    #compute statistical properties of the saved values
    for index in range(total_number_of_chosen_profiles):
        temp_flux = Shih_flux_list[index]
        number_of_fluxes_over_the_threshold += np.sum(np.abs(temp_flux)>maximum_reasonable_flux)
        number_of_zero_flux += np.sum(np.abs(temp_flux)==0)
        amount_of_missing_values += np.sum(np.isnan(temp_flux))
        #count the number of flux data points        
        total_number_of_fluxes += temp_flux.size
        
        temp_flux[np.abs(temp_flux)>maximum_reasonable_flux] = np.nan
    
        mean_flux[index] = np.nanmean(temp_flux)
        median_flux[index] = np.nanmedian(temp_flux)
        
        upper_percentile_flux[index] = np.nanpercentile(temp_flux, flux_percentile)
        lower_percentile_flux[index] = np.nanpercentile(temp_flux, 100-flux_percentile)
        second_upper_percentile_flux[index] = np.nanpercentile(temp_flux, second_flux_percentile)
        second_lower_percentile_flux[index] = np.nanpercentile(temp_flux, 100-second_flux_percentile)
                        
        """        
        mean_min_flux[index] = np.nanmean(oxygen_flux_statistic[index][:,0],axis=0)
        median_min_flux[index] = np.nanmedian(oxygen_flux_statistic[index][:,0],axis=0)

        upper_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], flux_percentile)
        lower_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], 100-flux_percentile)
        """        
        #bathymetrie_mean[index] = np.nanmean(bathymetrie_statistic[index])
        
        log_mean_dissipation[index] = np.nanmean(np.log10(dissipation_list[index]))
        arith_mean_dissipation[index] = np.log10(np.nanmean(dissipation_list[index]))
        median_dissipation[index] = np.log10(np.nanmedian(dissipation_list[index]))
        upper_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], dissip_percentile))
        lower_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], 100-dissip_percentile))
        second_upper_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], second_dissip_percentile))
        second_lower_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], 100-second_dissip_percentile))
                    
                       
                           
    rolling_mean_flux = [None] * total_number_of_chosen_profiles
    rolling_median_flux = [None] * total_number_of_chosen_profiles
    rolling_upper_percentile_flux = [None] * total_number_of_chosen_profiles
    rolling_lower_percentile_flux = [None] * total_number_of_chosen_profiles
    rolling_second_upper_percentile_flux = [None] * total_number_of_chosen_profiles
    rolling_second_lower_percentile_flux = [None] * total_number_of_chosen_profiles
    
    rolling_log_mean_dissipation = [None] * total_number_of_chosen_profiles
    rolling_arith_mean_dissipation = [None] * total_number_of_chosen_profiles
    rolling_median_dissipation = [None] * total_number_of_chosen_profiles
    rolling_lower_percentile_dissip = [None] * total_number_of_chosen_profiles
    rolling_upper_percentile_dissip = [None] * total_number_of_chosen_profiles
    rolling_second_upper_percentile_dissip = [None] * total_number_of_chosen_profiles
    rolling_second_lower_percentile_dissip = [None] * total_number_of_chosen_profiles
            
    max_longitude_gap_index = np.argmax(np.diff(longitude_list))
    max_longitude_gap = np.diff(longitude_list)[max_longitude_gap_index]
    print("LONGITUDE GAP",max_longitude_gap)
    #compute rolling average
    for index in range(total_number_of_chosen_profiles):
    
        #controls that the mean is not computed over too distant points
        if max_longitude_gap > 0.01:
            if ((index+rolling_window_size//2) >= max_longitude_gap_index+1) and ((index-rolling_window_size//2)<=max_longitude_gap_index+1):
                
                rolling_mean_flux[index] = np.nan
                rolling_median_flux[index] = np.nan  
                rolling_upper_percentile_flux[index] = np.nan
                rolling_lower_percentile_flux[index] = np.nan
                rolling_second_upper_percentile_flux[index] = np.nan
                rolling_second_lower_percentile_flux[index] = np.nan
                                
                rolling_log_mean_dissipation[index] = np.nan
                rolling_arith_mean_dissipation[index] = np.nan
                rolling_median_dissipation[index] = np.nan   
                rolling_lower_percentile_dissip[index] =  np.nan
                rolling_upper_percentile_dissip[index] =  np.nan
                rolling_second_upper_percentile_dissip[index] =  np.nan
                rolling_second_lower_percentile_dissip[index] =  np.nan
                continue
    
        if rolling_window_size ==1:
            rolling_mean_flux[index] = mean_flux[index]
            rolling_median_flux[index] = mean_flux[index]



            #print(index,longitude_list[index],np.round(mean_flux[index],3))
            
        else:
            try:
                """
                temp_array = mean_flux[index-(rolling_window_size//2):index+rolling_window_size//2]
                weights = (rolling_window_size - np.abs(np.arange(-(rolling_window_size//2),+rolling_window_size//2))) / np.sum((rolling_window_size - np.abs(np.arange(-(rolling_window_size//2),+rolling_window_size//2))))
                
                
                #print(index,rolling_window_size//2)
                #print(index-(rolling_window_size//2),index+rolling_window_size//2)
                #print(np.shape(temp_array),np.shape(weights))
                #print(np.shape(np.abs(np.arange(-(rolling_window_size//2,)+rolling_window_size//2))))
                #print(np.sum(weights))
                
                
                assert((np.sum(weights)-1)<0.0001)
                rolling_mean_flux[index] = np.sum(np.asarray(weights) * np.asarray(temp_array))
                """
                
                rolling_mean_flux[index] = np.nanmean(mean_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_median_flux[index] = np.nanmean(median_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_upper_percentile_flux[index] = np.nanmean(upper_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_lower_percentile_flux[index] = np.nanmean(lower_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_second_upper_percentile_flux[index] = np.nanmean(second_upper_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_second_lower_percentile_flux[index] = np.nanmean(second_lower_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                
                #print("test",np.nanmean(log_mean_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2]))
                rolling_log_mean_dissipation[index] = np.nanmean(log_mean_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_arith_mean_dissipation[index] = np.nanmean(arith_mean_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_median_dissipation[index] = np.nanmean(median_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_lower_percentile_dissip[index] =  np.nanmean(lower_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_upper_percentile_dissip[index] =  np.nanmean(upper_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_second_upper_percentile_dissip[index] = np.nanmean(second_upper_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_second_lower_percentile_dissip[index] = np.nanmean(second_lower_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
    
                #print(index,longitude_list[index],np.round(mean_flux[index-rolling_window_size//2:index+rolling_window_size//2],3))
    

            except (IndexError,ValueError):
                rolling_mean_flux[index] = np.nan
                rolling_median_flux[index] = np.nan  
                rolling_upper_percentile_flux[index] = np.nan
                rolling_lower_percentile_flux[index] = np.nan
                rolling_second_upper_percentile_flux[index] = np.nan
                rolling_second_lower_percentile_flux[index] = np.nan
                                
                rolling_log_mean_dissipation[index] = np.nan
                rolling_arith_mean_dissipation[index] = np.nan
                rolling_median_dissipation[index] = np.nan   
                rolling_lower_percentile_dissip[index] =  np.nan
                rolling_upper_percentile_dissip[index] =  np.nan
                rolling_second_upper_percentile_dissip[index] =  np.nan
                rolling_second_lower_percentile_dissip[index] =  np.nan 
                
    
    if density_interval == True:
        unit = "pot_density"
    else:
        unit = "pressure"
    np.savez("/home/ole/Thesis/Analysis/mss/data/"+cruisename+"_bathymetry_"+unit, longitude = bathymetry_longitude_list , bathymetry = bathymetry_list, interval = interval_list, unit = unit)
    
    print("total_number_of_valid_profiles",total_number_of_valid_profiles)             
    print("total_number_of_chosen_profiles",total_number_of_chosen_profiles)     
    print("number_of_fluxes_over_the_threshold\ttotal_number_of_fluxes\tratio")
    print("NaN",amount_of_missing_values,total_number_of_fluxes,100*amount_of_missing_values/total_number_of_fluxes,"%")
    print("0",number_of_zero_flux,total_number_of_fluxes,100*number_of_zero_flux/total_number_of_fluxes,"%")
    print(">",number_of_fluxes_over_the_threshold,total_number_of_fluxes,100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    print("Sum:",100*amount_of_missing_values/total_number_of_fluxes + 100*number_of_zero_flux/total_number_of_fluxes + 100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
        
        
       
    f_flux,flux_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
    bathymetrie_axes = flux_axarr.twinx()
    
    #define flux_axarr as the foreground
    flux_axarr.set_zorder(10)
    flux_axarr.patch.set_visible(False)
        
    
    flux_axarr.plot(longitude_list,mean_flux,"x", c = "tab:green", alpha = 0.4, label = "mean flux")                
    flux_axarr.plot(longitude_list,rolling_mean_flux, lw = 3, zorder = 3, c = "k", label = "rolling mean flux")#"tab:blue")
    #flux_axarr.plot(longitude_list,rolling_median_flux,"--",zorder = 3, c = "k")#"tab:blue")
    

    flux_axarr.fill_between(longitude_list,rolling_upper_percentile_flux,rolling_lower_percentile_flux, color = "tab:blue", zorder = 2, alpha = 0.7, label = str(flux_percentile)+"% percentile")
    flux_axarr.fill_between(longitude_list,rolling_upper_percentile_flux,rolling_second_upper_percentile_flux, color = "tab:blue", zorder = 2, alpha = 0.4, label = str(second_flux_percentile)+"% percentile")
    flux_axarr.fill_between(longitude_list,rolling_lower_percentile_flux,rolling_second_lower_percentile_flux, color = "tab:blue", zorder = 2, alpha = 0.4)
    

    bathymetrie_axes.set_ylim((np.nanmin(bathymetry_list)-5,np.nanmax(bathymetry_list)))
    bathymetrie_axes.invert_yaxis()
    bathymetrie_axes.set_ylabel("pressure [dbar]")
    
    #bathymetrie_axes.plot(bathymetry_longitude_list,bathymetry_list)
    bathymetrie_axes.fill_between(bathymetry_longitude_list,bathymetry_list, np.ones(len(bathymetry_list))*max(bathymetry_list),color = "lightgrey", zorder = 1, alpha = 0.8, label = "bathymetry")
    
    #bathymetrie_axes.plot(bathymetry_longitude_list,interval_list[:,0])
    #bathymetrie_axes.plot(bathymetry_longitude_list,interval_list[:,1])
    bathymetrie_axes.fill_between(bathymetry_longitude_list,interval_list[:,0],interval_list[:,1],color = "tab:red", alpha = 0.5)
    
    """
    mean_label = mlines.Line2D([], [], color='k', label='mean')
    median_label = mlines.Line2D([], [], ls = "--", color='k', label='median')
    down_label = mlines.Line2D([], [], color='tab:green', label='maximum downwards BB flux')
    up_label = mlines.Line2D([], [], color='tab:blue', label='maximum upwards BB flux')
    
    max_flux_label =  mpatches.Patch(color='tab:green', alpha = 0.6,label='downwards flux '+str(flux_percentile)+"% percentile")
    min_flux_label =  mpatches.Patch(color='tab:blue', alpha = 0.6, label='upwards flux '+str(flux_percentile)+"% percentile")
    second_max_flux_label =  mpatches.Patch(color='tab:green', alpha = 0.4,label='downwards flux '+str(second_flux_percentile)+"% percentile")
    second_min_flux_label =  mpatches.Patch(color='tab:blue', alpha = 0.4, label='upwards flux '+str(second_flux_percentile)+"% percentile")
    """
       
    bathymetrie_label =  mpatches.Patch(color='lightgrey', label='bathymetry')
    flux_axarr.legend(handles=[bathymetrie_label]) #loc=8

    #flux_axarr.set_ylim((-20,1))
    if cruisename == "emb177":
        flux_axarr.set_ylim((-50,1))    
    elif cruisename == "emb169":
        pass
        #flux_axarr.set_ylim((-50,1))  
                
    flux_axarr.set_xlabel(r"longitude [$\degree$]")    
    flux_axarr.set_ylabel(r"oxygen flux [mmol/(m$^2$*d]")
    flux_axarr.legend()
    
    f_flux.set_size_inches(18,10.5)
    f_flux.tight_layout() 
    f_flux.subplots_adjust(top=0.94)

    
    if flux_through_halocline == True:
        if density_interval == True:
            f_flux.suptitle(cruisename+": rolling mean Shih oxygen flux (over "+str(rolling_window_size)+" points) around the halocline ("+str(upper_bound_halocline_as_density)+r"< $\sigma$ < "+str(lower_bound_halocline_as_density)+r" kg/m$^3$)") 
            f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_running_mean_halocline_flux_density")           
        else:       
            f_flux.suptitle(cruisename+": rolling mean Shih oxygen flux (over "+str(rolling_window_size)+" points) around the halocline ("+str(upper_bound_halocline_in_db)+"-"+str(lower_bound_halocline_in_db)+"dbar)") 
            f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_running_mean_halocline_flux_db")
    else:
        f_flux.suptitle(cruisename+": rolling mean Shih oxygen flux in the lowermost "+str(height_above_ground)+" meters above ground")
        f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_running_mean_BBL_flux")

    
        
    ###############################################################################################################
    

    f_dissip,dissip_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
     #define dissip_axarr as the foreground
    dissip_axarr.set_zorder(10)
    dissip_axarr.patch.set_visible(False)
    
    bathymetrie_axes2 = dissip_axarr.twinx()
    bathymetrie_axes2.set_ylim((np.nanmin(bathymetry_list)-5,np.nanmax(bathymetry_list)))
    bathymetrie_axes2.invert_yaxis()
    bathymetrie_axes2.set_ylabel("pressure [dbar]")
    
    #bathymetrie_axes2.plot(bathymetry_longitude_list,bathymetry_list)    
    bathymetrie_axes2.fill_between(bathymetry_longitude_list,bathymetry_list, np.ones(len(bathymetry_list))*max(bathymetry_list),color = "lightgrey", zorder = 1, alpha = 0.8, label = "bathymetry")
    
    #bathymetrie_axes2.plot(bathymetry_longitude_list,interval_list[:,0])
    #bathymetrie_axes2.plot(bathymetry_longitude_list,interval_list[:,1])
    bathymetrie_axes2.fill_between(bathymetry_longitude_list,interval_list[:,0],interval_list[:,1],color = "tab:red", alpha = 0.5)
    
    dissip_axarr.plot(longitude_list,log_mean_dissipation,"x", c = "tab:green", alpha = 0.4, label = "mean dissipation")        
    dissip_axarr.plot(longitude_list,rolling_log_mean_dissipation, "k", label ="rolling logarithmic mean")
    dissip_axarr.plot(longitude_list,rolling_arith_mean_dissipation, "k--", label ="rolling arithmic mean")

    dissip_axarr.fill_between(longitude_list,rolling_upper_percentile_dissip,rolling_lower_percentile_dissip, color = "tab:blue", alpha = 0.7, label=str(dissip_percentile)+"% percentile") 
    dissip_axarr.fill_between(longitude_list,rolling_upper_percentile_dissip,rolling_second_upper_percentile_dissip, color = "tab:blue", alpha = 0.4, label=str(second_dissip_percentile)+"% percentile")
    dissip_axarr.fill_between(longitude_list,rolling_lower_percentile_dissip,rolling_second_lower_percentile_dissip, color = "tab:blue", alpha = 0.4) 

                
    dissip_axarr.set_ylabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
    dissip_axarr.set_xlabel(r"longitude [$\degree$E]")       
    """      
    bathymetrie_label =  mpatches.Patch(color='lightgrey', label='bathymetrie')
    dissip_mean_label = mlines.Line2D([], [], color= "k", label='mean dissipation $\epsilon$') #tab:blue
    dissip_median_label = mlines.Line2D([], [], color = "k", ls = "--", label='median dissipation $\epsilon$')
    dissip_percent_label =  mpatches.Patch(color='tab:blue', label=str(dissip_percentile)+"% percentile")
    dissip_second_percent_label =  mpatches.Patch(color='tab:blue', alpha = 0.4, label=str(second_dissip_percentile)+"% percentile")           
    dissip_axarr.legend(handles=[dissip_mean_label,dissip_median_label,dissip_percent_label,dissip_second_percent_label,bathymetrie_label])
    """
    #f_dissip.suptitle(cruisename+": mean dissipation in the lowermost "+str(height_above_ground)+" meters above ground")
    #f_dissip.suptitle(cruisename+": mean dissipation around the halocline (67-77dbar) ("+str(number_of_intervals)+" intervals)")

    dissip_axarr.legend()
    
    f_dissip.set_size_inches(18,10.5)
    f_dissip.tight_layout() 
    f_dissip.subplots_adjust(top=0.94)
    if flux_through_halocline == True:
        if density_interval == True:
            titlestring = cruisename+": rolling mean dissipation (over "+str(rolling_window_size)+" points) around the halocline ("+str(upper_bound_halocline_as_density)+r"< $\sigma$ < "+str(lower_bound_halocline_as_density)+r" kg/m$^3$)"          
        else:          
            titlestring = cruisename+": rolling mean dissipation (over "+str(rolling_window_size)+" points) around the halocline ("+str(upper_bound_halocline_in_db)+"-"+str(lower_bound_halocline_in_db)+"dbar)"
        f_dissip.suptitle(titlestring)    
        f_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_running_mean_halocline_dissipation")
    else:
        f_dissip.suptitle(cruisename+": rolling mean dissipation in the lowermost "+str(height_above_ground)+" meters above ground")
        f_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_running_mean_BBL_dissipation")
       
    plt.show()
    
    
    
    
    
