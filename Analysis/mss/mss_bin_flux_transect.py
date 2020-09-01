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

import geopy.distance as geo
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

    
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177"] #"/home/ole/windows/processed_mss/emb217"]#,"/home/ole/windows/processed_mss/emb169",
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217","/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb169"]

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
maximum_reasonable_flux = 500 #float('Inf') #200 #Fluxes with absolute values above this cut off value will be discarded
acceptable_slope = 2 #float('Inf') #acceptable bathymetry difference in dbar between two neighboring data points. 

flux_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_flux_percentile = 97.72
dissip_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_dissip_percentile = 97.72

number_of_dissipation_subplots = 1 #Decide if both the mean and the median subplots is shown or only the mean


f_flux,flux_axarr = plt.subplots(nrows = 2, ncols = 1, sharex = True) 

f_dissip,dissip_axarr = plt.subplots(nrows = 2, ncols = 1, sharex = True) 
textstr = ""

list_of_bad_profiles,reasons = np.loadtxt("./data/list_of_bad_profiles.txt", dtype=str, unpack=True)
         
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    number_of_fluxes_over_the_threshold = 0
    number_of_transects = 0
    total_number_of_fluxes = 0
    number_of_zero_flux = 0
    amount_of_missing_values = 0
    total_number_of_valid_profiles = 0    
    total_number_of_correct_profiles = 0
    
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
        upper_bound_halocline_as_density = 1006.9 #1006.9
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
    
        if "_".join((cruisename,transect_name)) in list_of_bad_profiles:
            print("_".join((cruisename,transect_name)),"skipped")
            continue
            
                
        print("\n",transect_name)
            
        
        data = np.load(datafile_path)
        
        try:
            number_of_profiles = data["number_of_profiles"] #
            lat = data["lat"] #Latitude of the profiles
            lon = data["lon"] #Longitude of the profiles
            distance = data["distance"] #distance from the starting profile (monotonically increasing)
            
            #interp_pressure = data["interp_pressure"]
            #oxygen_grid = data["oxygen_grid"]
            #salinity_grid = data["salinity_grid"]
            #consv_temperature_grid = data["consv_temperature_grid"]
            #density_grid = data["density_grid"]
            
            eps_pressure = data["bin_pressure"]
            eps_grid = data["bin_eps_grid"]
            corrected_eps_grid = data["corrected_bin_eps_grid"]
            eps_consv_temperature_grid = data["bin_consv_temperature_grid"]
            eps_oxygen_sat_grid = data["bin_oxygen_sat_grid"]
            eps_oxygen_grid = data["bin_oxygen_grid"] 
            
            eps_N_squared_grid = data["bin_N_squared_grid"]
            eps_density_grid = data["bin_density_grid"]
            eps_pot_density_grid = data["bin_pot_density_grid"]
            #eps_viscosity_grid = data["eps_viscosity_grid"]
            eps_Reynolds_bouyancy_grid = data["bin_Reynolds_bouyancy_grid"]
            corrected_eps_Reynolds_bouyancy_grid = data["corrected_bin_Reynolds_bouyancy_grid"]

            
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
            """
        
        except KeyError:
            print(transect_name,"Error during loading data")
            raise AssertionError
            
        number_of_transects+=1 
            
        print("Number of profiles:",number_of_profiles)
        
        #print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
        
        #calculate the idices of the bottom and some meters above that
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid,height_above_ground = height_above_ground)
        """
        bathymetry                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetry_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
        """
        bathymetry,list_of_bathymetry_indices = results[0]
        #BBL,list_of_BBL_indices = results[1] #not needed here
        BBL_range,list_of_BBL_range_indices = results[2]
        
        eps_N_grid = np.sqrt(eps_N_squared_grid)
        #ozmidov scale
        ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))
        
        #conversion from pressure coordinates to depth
        eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
        bathymetry_in_m = gsw.z_from_p(bathymetry,np.mean(lat))
        
        eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
        
        distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetry_in_m,(-1,1))
        distance_from_ground_grid[distance_from_ground_grid < 0] = np.nan
        boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)
        
        turbulent_diffusivity_Osborn_grid = thesis.get_turbulent_diffusivity_Osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
        turbulent_diffusivity_BB_grid = thesis.get_turbulent_diffusivity_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
        turbulent_diffusivity_Shih_grid = thesis.get_turbulent_diffusivity_Shih(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)

        """
        shear_velocity_grid = thesis.get_shear_velocity(eps_grid,distance_from_ground_grid)
        boolean_array_for_law_of_the_wall = distance_from_ground_grid < 1.5
        
        turbulent_diffusivity_Osborn_grid[boolean_array_for_law_of_the_wall] = law_of_the_wall_turbulent_diffusivity_Osborn_grid[boolean_array_for_law_of_the_wall]
        turbulent_diffusivity_BB_grid[boolean_array_for_law_of_the_wall] = law_of_the_wall_turbulent_diffusivity_BB_grid[boolean_array_for_law_of_the_wall]
        turbulent_diffusivity_Shih_grid[boolean_array_for_law_of_the_wall] = law_of_the_wall_turbulent_diffusivity_Shih_grid[boolean_array_for_law_of_the_wall]
        """
                
        oxygen_gradient_grid = thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        unit_conversion_grid = 86400*(1000/eps_density_grid) #to convert from m*micromol/(kg*s) to mmol/(m^2*d)
    
        oxygen_flux_Osborn_grid = - turbulent_diffusivity_Osborn_grid * oxygen_gradient_grid * unit_conversion_grid
        oxygen_flux_BB_grid =  - turbulent_diffusivity_BB_grid * oxygen_gradient_grid * unit_conversion_grid
        oxygen_flux_Shih_grid = - turbulent_diffusivity_Shih_grid * oxygen_gradient_grid * unit_conversion_grid
       
       
       
        
        for profile in range(number_of_profiles):
        
        
            if "_".join((cruisename,transect_name,str(profile))) in list_of_bad_profiles:
                print("_".join((cruisename,transect_name,str(profile))),"skipped")
                continue
                

            from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-upper_bound_halocline_as_density))     
            to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-lower_bound_halocline_as_density))
            
            """
            sea_floor_index =  -np.argmax(np.flip(~np.isnan(eps_grid[profile,:]))) #at the moment the index is negative
            sea_floor_index = eps_grid[profile,:].size + sea_floor_index #now defined as positive index
            
           
            if sea_floor_index == eps_grid[profile,:].size:
                if not np.isnan(eps_grid[profile,-1]): #if there are no NAN values towards the bottom
                    sea_floor_index = len(eps_pressure)-1 #set the last index as the index of the bottom
        
            #sea_floor_pressure = eps_pressure[sea_floor_index]
            
            to_index = max(from_index,min(sea_floor_index-3,to_index))
            """
                
            #find the correct position in the sorted list
            for index,value in enumerate(bathymetry_longitude_list):
                if value > lon[profile]:
                    list_position = index
                    break
                elif index == len(bathymetry_longitude_list)-1:
                    list_position = len(bathymetry_longitude_list)
                    break
                                
            if len(bathymetry_longitude_list) == 0:   
                bathymetry_list.append(bathymetry[profile])
                bathymetry_longitude_list.append(lon[profile])
                
                if from_index == to_index:
                    interval_list.append([np.nan,np.nan])                
                else:
                    interval_list.append([eps_pressure[from_index],eps_pressure[to_index]])
            else:
                bathymetry_list.insert(list_position,bathymetry[profile])
                bathymetry_longitude_list.insert(list_position,lon[profile])
                if from_index == to_index:
                    interval_list.insert(list_position,[np.nan,np.nan])
                else:
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
                #BB_flux_list.append(oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.append(oxygen_flux_Shih_grid[profile,from_index:to_index])
                Osborn_flux_list.append(oxygen_flux_Osborn_grid[profile,from_index:to_index])
                longitude_list.append(lon[profile])
            
            else:
                
                #Sort the current profile into the list            
                dissipation_list.insert(list_position,eps_grid[profile,from_index:to_index])
                #BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.insert(list_position,oxygen_flux_Shih_grid[profile,from_index:to_index])
                Osborn_flux_list.insert(list_position,oxygen_flux_Osborn_grid[profile,from_index:to_index])
                longitude_list.insert(list_position,lon[profile])
    
            assert(np.all(longitude_list == sorted(longitude_list)))

            total_number_of_valid_profiles+=1

        
    ###########################################################################################################################
    print(len(longitude_list),"used profiles")
    assert(len(longitude_list) != 0)
    assert(np.all(longitude_list == sorted(longitude_list)))
    assert(np.all(bathymetry_longitude_list == sorted(bathymetry_longitude_list)))        
    interval_list = np.asarray(interval_list)
    
    #compute mean and std over the saved intervals
    mean_Osborn_flux = [None] * total_number_of_valid_profiles
    mean_Shih_flux = [None] * total_number_of_valid_profiles
    median_flux = [None] * total_number_of_valid_profiles
    upper_percentile_flux = [None] * total_number_of_valid_profiles
    lower_percentile_flux = [None] * total_number_of_valid_profiles
    second_upper_percentile_flux = [None] * total_number_of_valid_profiles
    second_lower_percentile_flux = [None] * total_number_of_valid_profiles
  
    #bathymetry_mean = [None] * number_of_intervals
    #bathymetry_percentile = [None] * number_of_intervals
    
    log_mean_dissipation = [None] * total_number_of_valid_profiles
    arith_mean_dissipation = [None] * total_number_of_valid_profiles
    median_dissipation = [None] * total_number_of_valid_profiles
    lower_percentile_dissip = [None] * total_number_of_valid_profiles
    upper_percentile_dissip = [None] * total_number_of_valid_profiles

    second_lower_percentile_dissip = [None] * total_number_of_valid_profiles
    second_upper_percentile_dissip = [None] * total_number_of_valid_profiles

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
    for index in range(total_number_of_valid_profiles):
        temp_Shih_flux = Shih_flux_list[index]
        number_of_fluxes_over_the_threshold += np.sum(np.abs(temp_Shih_flux)>maximum_reasonable_flux)
        number_of_zero_flux += np.sum(np.abs(temp_Shih_flux)==0)
        amount_of_missing_values += np.sum(np.isnan(temp_Shih_flux))
        #count the number of flux data points        
        total_number_of_fluxes += temp_Shih_flux.size
        
        temp_Shih_flux[np.abs(temp_Shih_flux)>maximum_reasonable_flux] = np.nan
    
        mean_Shih_flux[index] = np.nanmean(temp_Shih_flux)
        median_flux[index] = np.nanmedian(temp_Shih_flux)
        
        temp_Osborn_flux = Osborn_flux_list[index]
        temp_Osborn_flux[np.abs(temp_Osborn_flux)>maximum_reasonable_flux] = np.nan  
        mean_Osborn_flux[index] = np.nanmean(temp_Osborn_flux)
        
        upper_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, flux_percentile)
        lower_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, 100-flux_percentile)
        second_upper_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, second_flux_percentile)
        second_lower_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, 100-second_flux_percentile)
                        
        """        
        mean_min_flux[index] = np.nanmean(oxygen_flux_statistic[index][:,0],axis=0)
        median_min_flux[index] = np.nanmedian(oxygen_flux_statistic[index][:,0],axis=0)

        upper_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], flux_percentile)
        lower_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], 100-flux_percentile)
        """        
        #bathymetry_mean[index] = np.nanmean(bathymetry_statistic[index])
        
        log_mean_dissipation[index] = np.nanmean(np.log10(dissipation_list[index]))
        arith_mean_dissipation[index] = np.log10(np.nanmean(dissipation_list[index]))
        median_dissipation[index] = np.log10(np.nanmedian(dissipation_list[index]))
        upper_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], dissip_percentile))
        lower_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], 100-dissip_percentile))
        second_upper_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], second_dissip_percentile))
        second_lower_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], 100-second_dissip_percentile))
                    
                       
                           
    rolling_mean_Shih_flux = [None] * total_number_of_valid_profiles
    rolling_mean_Osborn_flux = [None] * total_number_of_valid_profiles
    rolling_median_flux = [None] * total_number_of_valid_profiles
    rolling_upper_percentile_flux = [None] * total_number_of_valid_profiles
    rolling_lower_percentile_flux = [None] * total_number_of_valid_profiles
    rolling_second_upper_percentile_flux = [None] * total_number_of_valid_profiles
    rolling_second_lower_percentile_flux = [None] * total_number_of_valid_profiles
    
    rolling_log_mean_dissipation = [None] * total_number_of_valid_profiles
    rolling_arith_mean_dissipation = [None] * total_number_of_valid_profiles
    rolling_median_dissipation = [None] * total_number_of_valid_profiles
    rolling_lower_percentile_dissip = [None] * total_number_of_valid_profiles
    rolling_upper_percentile_dissip = [None] * total_number_of_valid_profiles
    rolling_second_upper_percentile_dissip = [None] * total_number_of_valid_profiles
    rolling_second_lower_percentile_dissip = [None] * total_number_of_valid_profiles
            
    #max_longitude_gap_index = np.argmax(np.diff(longitude_list))
    #max_longitude_gap = np.diff(longitude_list)[max_longitude_gap_index]
    #print("LONGITUDE GAP",max_longitude_gap)
    #compute rolling average
    for index in range(total_number_of_valid_profiles):

        #controls that the mean is not computed over too distant points
        number_of_nans_in_averaging_window = np.count_nonzero(np.isnan(mean_Shih_flux[index-(rolling_window_size//2):index+rolling_window_size//2])) 
        if number_of_nans_in_averaging_window > 3/4 * rolling_window_size:
 
            rolling_mean_Shih_flux[index] = np.nan
            rolling_mean_Osborn_flux[index] = np.nan
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
            rolling_mean_Shih_flux[index] = mean_Shih_flux[index]
            rolling_mean_Osborn_flux[index] = mean_Osborn_flux[index]
            rolling_median_flux[index] = mean_flux[index]



            #print(index,longitude_list[index],np.round(mean_flux[index],3))
            
        else:
            try:
                
                rolling_mean_Shih_flux[index] = np.nanmean(mean_Shih_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_mean_Osborn_flux[index] = np.nanmean(mean_Osborn_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
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
                rolling_mean_Shih_flux[index] = np.nan
                rolling_mean_Osborn_flux[index] = np.nan
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
                
    print(cruisename,"total number of transects =",number_of_transects)  
    print("total_number_of profiles",total_number_of_correct_profiles)        
    print("total_number_of_valid_profiles",total_number_of_valid_profiles)     
    print("number_of_fluxes_over_the_threshold\ttotal_number_of_fluxes\tratio")
    print("NaN",amount_of_missing_values,total_number_of_fluxes,100*amount_of_missing_values/total_number_of_fluxes,"%")
    print("0",number_of_zero_flux,total_number_of_fluxes,100*number_of_zero_flux/total_number_of_fluxes,"%")
    print(">",number_of_fluxes_over_the_threshold,total_number_of_fluxes,100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    print("Sum:",100*amount_of_missing_values/total_number_of_fluxes + 100*number_of_zero_flux/total_number_of_fluxes + 100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    
    distance_list = np.zeros(np.shape(longitude_list))
    first_point = (np.mean(lat),longitude_list[0]) 
    for i in range(len(longitude_list)):
        current_point = (np.mean(lat),longitude_list[i])
        distance_list[i] = geo.geodesic(first_point,current_point).km #Distance in km
        

    delta_X = thesis.central_differences(distance_list)
    
    print("\n\n\n",cruisename,"flux sum:")
    print("Osborn rolling mean",np.nansum(rolling_mean_Osborn_flux*delta_X),"Osborn profiles", np.nansum(mean_Osborn_flux*delta_X))
    print("Shih rolling mean",np.nansum(rolling_mean_Shih_flux*delta_X),"Shih profiles", np.nansum(mean_Shih_flux*delta_X))
    print("\n\n\n")
    
    
    np.savetxt("./data/"+cruisename+'_bin_flux_results.txt', np.transpose([longitude_list,distance_list,bathymetry_list,mean_Osborn_flux,rolling_mean_Osborn_flux,mean_Shih_flux,rolling_mean_Shih_flux]), header = "longitude\tdistance\tdepth[dbar]\traw Osborn\trolling mean Osborn\traw Shih\trolling mean Shih", fmt = "%3.8f")
    
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
        
    if cruisename == "emb169":
        color = "tab:green"
        label_name = "autumn cruise" 
                
    if cruisename == "emb177":
        color = "tab:blue"
        label_name = "winter cruise" 
        
        
    if cruisename == "emb217":
        color = "tab:red"      
        label_name = "summer cruise" 
    
    
    if textstr != "":
        textstr = textstr + "\n"
    
    textstr = textstr + label_name +" "+cruisename+ ":\t"+ str(upper_bound_halocline_as_density)+r" < σ < "+str(lower_bound_halocline_as_density)+r" kg/m$^3$"
    
    
                  
    flux_axarr[0].plot(longitude_list,rolling_mean_Shih_flux, lw = 2.5, zorder = 3, c = color)#, label = label_name)#"tab:blue")
    flux_axarr[0].plot(longitude_list,rolling_mean_Osborn_flux, ls = ":", lw = 2.5, zorder = 3, c = color)#, label = label_name)
    
    #flux_axarr[0].plot(longitude_list,mean_Shih_flux,"x", c = "tab:green", alpha = 0.4, label = "mean dissipation")        
    
    #flux_axarr.fill_between(longitude_list,rolling_upper_percentile_flux,rolling_lower_percentile_flux, color = "tab:blue", zorder = 2, alpha = 0.7, label = str(flux_percentile)+"% percentile")
    #flux_axarr.fill_between(longitude_list,rolling_upper_percentile_flux,rolling_second_upper_percentile_flux, color = "tab:blue", zorder = 2, alpha = 0.4, label = str(second_flux_percentile)+"% percentile")
    #flux_axarr.fill_between(longitude_list,rolling_lower_percentile_flux,rolling_second_lower_percentile_flux, color = "tab:blue", zorder = 2, alpha = 0.4)
    


    
    #bathymetry_axes.plot(bathymetry_longitude_list,bathymetry_list)
    if cruisename == "emb217":
        flux_axarr[1].set_ylim((np.nanmin(bathymetry_list)-5,np.nanmax(bathymetry_list)))
        flux_axarr[1].invert_yaxis()
        flux_axarr[1].set_ylabel("pressure [dbar]")
        flux_axarr[1].fill_between(bathymetry_longitude_list,bathymetry_list, np.ones(len(bathymetry_list))*max(bathymetry_list),color = "lightgrey", zorder = 1, alpha = 0.8, label = "bathymetry")
    
    #bathymetry_axes.plot(bathymetry_longitude_list,interval_list[:,0])
    #bathymetry_axes.plot(bathymetry_longitude_list,interval_list[:,1])
    flux_axarr[1].fill_between(bathymetry_longitude_list,interval_list[:,0],interval_list[:,1],color = color, alpha = 0.5, label = label_name+" "+cruisename)
    


    
        
    ###############################################################################################################
    


    
    #bathymetry_axes2.plot(bathymetry_longitude_list,bathymetry_list) 
    if cruisename == "emb217": 
        dissip_axarr[1].set_ylim((np.nanmin(bathymetry_list)-5,np.nanmax(bathymetry_list)))
        dissip_axarr[1].invert_yaxis()
        dissip_axarr[1].set_ylabel("pressure [dbar]")  
        dissip_axarr[1].fill_between(bathymetry_longitude_list,bathymetry_list, np.ones(len(bathymetry_list))*max(bathymetry_list),color = "lightgrey", zorder = 1, alpha = 0.8, label = "bathymetry")
    
    #bathymetry_axes2.plot(bathymetry_longitude_list,interval_list[:,0])
    #bathymetry_axes2.plot(bathymetry_longitude_list,interval_list[:,1])
    dissip_axarr[1].fill_between(bathymetry_longitude_list,interval_list[:,0],interval_list[:,1],color = color, alpha = 0.5, label = label_name+" "+cruisename)
    
    #dissip_axarr[0].plot(longitude_list,log_mean_dissipation,"x", c = "tab:green", alpha = 0.4, label = "mean dissipation")        
    #dissip_axarr.plot(longitude_list,rolling_log_mean_dissipation, "k", label ="rolling logarithmic mean")
    dissip_axarr[0].plot(longitude_list,rolling_arith_mean_dissipation, lw = 2.5, c = color, label = label_name+" "+cruisename)

    #dissip_axarr.fill_between(longitude_list,rolling_upper_percentile_dissip,rolling_lower_percentile_dissip, color = "tab:blue", alpha = 0.7, label=str(dissip_percentile)+"% percentile") 
    #dissip_axarr.fill_between(longitude_list,rolling_upper_percentile_dissip,rolling_second_upper_percentile_dissip, color = "tab:blue", alpha = 0.4, label=str(second_dissip_percentile)+"% percentile")
    #dissip_axarr.fill_between(longitude_list,rolling_lower_percentile_dissip,rolling_second_lower_percentile_dissip, color = "tab:blue", alpha = 0.4) 


###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

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

emb169_label = mpatches.Patch(color='tab:green', label='autumn cruise')
emb177_label = mpatches.Patch(color='tab:blue', label='winter cruise')
emb217_label = mpatches.Patch(color='tab:red', label='summer cruise')
Osborn_label = mlines.Line2D([], [], ls = ":", lw = 2.5, c = "k", label = "Osborn oxygen flux")
Shih_label = mlines.Line2D([], [], ls = "-", lw = 2.5, c = "k", label = "Shih oxygen flux")

flux_axarr[0].legend(handles=[emb217_label,emb177_label,emb169_label,Osborn_label,Shih_label],loc = "lower left")
 
bathymetry_label =  mpatches.Patch(color='lightgrey', label='bathymetry')
flux_axarr[1].legend(loc = "lower left")


flux_axarr[0].set_ylim((-85,1))    
        
flux_axarr[1].set_xlabel(r"longitude [$\degree$E]")    
flux_axarr[0].set_ylabel(r"oxygen flux [mmol/(m$^2$d)]")


f_flux.set_size_inches(16,10.5)
f_flux.tight_layout() 
f_flux.subplots_adjust(top=0.94)

props = dict(boxstyle='square', facecolor = "white")
flux_axarr[1].text(0.62, 0.05, textstr, transform=flux_axarr[1].transAxes, fontsize=14,verticalalignment='bottom', bbox=props, multialignment = "right")

f_flux.suptitle("Oxygen flux through the halocline along the transect: rolling isopycnal mean over "+str(rolling_window_size)+" points", weight = "bold") # (binned data)")
f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/flux_bin_presentation", dpi = 300)           

   
#------------------------------------------------------------------------------------------------------------#
                        
dissip_axarr[0].set_ylabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
dissip_axarr[1].set_xlabel(r"longitude [$\degree$E]")       
"""      
bathymetry_label =  mpatches.Patch(color='lightgrey', label='bathymetry')
dissip_mean_label = mlines.Line2D([], [], color= "k", label='mean dissipation $\epsilon$') #tab:blue
dissip_median_label = mlines.Line2D([], [], color = "k", ls = "--", label='median dissipation $\epsilon$')
dissip_percent_label =  mpatches.Patch(color='tab:blue', label=str(dissip_percentile)+"% percentile")
dissip_second_percent_label =  mpatches.Patch(color='tab:blue', alpha = 0.4, label=str(second_dissip_percentile)+"% percentile")           
dissip_axarr.legend(handles=[dissip_mean_label,dissip_median_label,dissip_percent_label,dissip_second_percent_label,bathymetry_label])
"""
#f_dissip.suptitle(cruisename+": mean dissipation in the lowermost "+str(height_above_ground)+" meters above ground")
#f_dissip.suptitle(cruisename+": mean dissipation around the halocline (67-77dbar) ("+str(number_of_intervals)+" intervals)")

dissip_axarr[0].legend(loc = "upper left")
dissip_axarr[1].legend(loc = "lower left")

f_dissip.set_size_inches(16,10.5)
f_dissip.tight_layout() 
f_dissip.subplots_adjust(top=0.94)
 
    
props = dict(boxstyle='square', facecolor = "white")
dissip_axarr[1].text(0.65, 0.05, textstr, transform=dissip_axarr[1].transAxes, fontsize=14,verticalalignment='bottom', bbox=props, multialignment = "right")
        
f_dissip.suptitle("rolling isopycnal mean dissipation (over "+str(rolling_window_size)+" points) around the halocline")# (binned data)")
f_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/dissip_bin_presentation", dpi = 300)

   
plt.show()
    
    
    
    
    
