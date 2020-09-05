#######################################################
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


rolling_window_size = 12


#plot halocline (and other isolines) in pressure and density coordinates. 
f_iso,iso_axis = plt.subplots(2, sharex = True)

#plot mean oxygen flux dependent on in pressure and density coordinates (once for shih and once for Osborn)
f_box, box_axis = plt.subplots(2, sharex = True, sharey = True)

height_above_ground = 5 #Size of the averaging interval above ground for the BBL, has no meaning if (flux_through_halocline == True)
maximum_reasonable_flux = float('Inf') #200 #Fluxes with absolute values above this cut off value will be discarded
acceptable_slope = 2 #float('Inf') #acceptable bathymetry difference in dbar between two neighboring data points. 

#flux_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
#second_flux_percentile = 97.72
#dissip_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
#second_dissip_percentile = 97.72

number_of_dissipation_subplots = 1 #Decide if both the mean and the median subplots is shown or only the mean

#datafile_paths = ["/home/ole/windows/processed_mss/emb177/TS1_8.npz","/home/ole/windows/processed_mss/emb177/TS1_10.npz","/home/ole/windows/processed_mss/emb169/TS112.npz"]
datafile_paths = ["/home/ole/windows/processed_mss/emb169/TS112.npz"]
#datafile_paths = ["/home/ole/windows/processed_mss/emb177/TS1_10.npz"]

list_of_bad_profiles,reasons = np.loadtxt("./data/list_of_bad_profiles.txt", dtype=str, unpack=True)

interior_longitude_interval = [20.53,20.57]
edge_longitude_interval = [20.57,20.64]
      
box_sizes = np.arange(0.1,4,0.2)
box_sizes_in_dbar = np.zeros(box_sizes.size)
example_transect_longitude = 20.57
        
for datafile_index,datafile_path in enumerate(datafile_paths):

    transect_halocline_depths = []
    transect_halocline_densities = []


    splitted_filename = datafile_path.split("/")
    cruisename = splitted_filename[4][0:6]
    transect_name = splitted_filename[-1][:-4]

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

        
        eps_pressure = data["bin_pressure"]
        eps_grid = data["bin_eps_grid"]
        corrected_eps_grid = data["corrected_bin_eps_grid"]
        eps_consv_temperature_grid = data["bin_consv_temperature_grid"]
        eps_oxygen_sat_grid = data["bin_oxygen_sat_grid"]
        eps_oxygen_grid = data["bin_oxygen_grid"] 
        eps_salinity_grid = data["bin_salinity_grid"] 
                
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
        
        
    print("Number of profiles:",number_of_profiles)
    
    #print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
    
    #calculate the idices of the bottom and some meters above that
    results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid,height_above_ground = height_above_ground)

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
   
   
    interior_Shih_flux_for_different_intervals = []
    interior_Osborn_flux_for_different_intervals = []

    edge_Shih_flux_for_different_intervals = []
    edge_Osborn_flux_for_different_intervals = []
    

    
    example_profile_index = np.argmin(np.abs(lon-example_transect_longitude))
    
    for width_index,density_box_width in enumerate(box_sizes):
    
    
        dissipation_list = []
        BB_flux_list = []
        Shih_flux_list = []
        Osborn_flux_list = []
        longitude_list = []
        bathymetry_list = []
        bathymetry_longitude_list = []
        interval_list = []    
    
        total_number_of_valid_profiles = 0
    
        for profile in range(number_of_profiles):
    
    
            if "_".join((cruisename,transect_name,str(profile))) in list_of_bad_profiles:
                if profile == example_profile_index: example_profile_index+=1
                print("_".join((cruisename,transect_name,str(profile))),"skipped")
                continue
            
            total_number_of_valid_profiles+=1
                
            halocline_depth,halocline_density = thesis.get_halocline_and_halocline_density(eps_pressure,eps_oxygen_sat_grid[profile],eps_salinity_grid[profile],eps_consv_temperature_grid[profile],eps_pot_density_grid[profile])
        
            transect_halocline_depths.append(halocline_depth)
            transect_halocline_densities.append(halocline_density)
            
            if np.isnan(halocline_depth):
                #find the correct position in the sorted list
                for index,value in enumerate(longitude_list):
                    if value > lon[profile]:
                        list_position = index
                        break
                    elif index == len(longitude_list)-1:
                        list_position = len(longitude_list)
                        break
                    
            
                print("no halocline",lon[profile])
                if len(longitude_list) == 0:    
                    dissipation_list.append(np.nan)
                    #BB_flux_list.append(oxygen_flux_BB_grid[profile,from_index:to_index])
                    Shih_flux_list.append(np.nan)
                    Osborn_flux_list.append(np.nan)
                    longitude_list.append(lon[profile])
                    bathymetry_list.append(bathymetry[profile])
                                                
                else:
                    
                    #Sort the current profile into the list            
                    dissipation_list.insert(list_position,np.nan)
                    #BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                    Shih_flux_list.insert(list_position,np.nan)
                    Osborn_flux_list.insert(list_position,np.nan)
                    longitude_list.insert(list_position,lon[profile])
                    bathymetry_list.insert(list_position,bathymetry[profile])
                continue
        
            #choose the vertical averaging interval dependent on the box size
            from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile]) - (halocline_density - density_box_width)))     
            to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])- (halocline_density + density_box_width)))
           
            if profile == example_profile_index:
                box_sizes_in_dbar[width_index] = np.abs(eps_pressure[to_index] - eps_pressure[from_index])
            
           
           
            #print(from_index,to_index)
            assert from_index <= to_index
            
            #from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-upper_bound_halocline_as_density))     
            #to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-lower_bound_halocline_as_density))
            
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
                bathymetry_list.append(bathymetry[profile])
                                            
            else:
                
                #Sort the current profile into the list            
                dissipation_list.insert(list_position,eps_grid[profile,from_index:to_index])
                #BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.insert(list_position,oxygen_flux_Shih_grid[profile,from_index:to_index])
                Osborn_flux_list.insert(list_position,oxygen_flux_Osborn_grid[profile,from_index:to_index])
                longitude_list.insert(list_position,lon[profile])
                bathymetry_list.insert(list_position,bathymetry[profile])
                
                #print(np.nanmean(oxygen_flux_Shih_grid[profile,from_index:to_index]))

            

        
        ###########################################################################################################################
        print(len(longitude_list),"used profiles")
        assert(len(longitude_list) != 0)
        assert(np.all(longitude_list == sorted(longitude_list)))      
        interval_list = np.asarray(interval_list)

        #compute mean and std over the saved intervals
        mean_Osborn_flux = [None] * total_number_of_valid_profiles
        mean_Shih_flux = [None] * total_number_of_valid_profiles
        arith_mean_dissipation = [None] * total_number_of_valid_profiles


        #average vertically in the denisty interval
        for index in range(len(longitude_list)):
            temp_Shih_flux = np.asarray(Shih_flux_list[index])    
            temp_Shih_flux[np.abs(temp_Shih_flux)>maximum_reasonable_flux] = np.nan
            mean_Shih_flux[index] = np.nanmean(temp_Shih_flux)
                        
            temp_Osborn_flux = np.asarray(Osborn_flux_list[index])
            temp_Osborn_flux[np.abs(temp_Osborn_flux)>maximum_reasonable_flux] = np.nan  
            mean_Osborn_flux[index] = np.nanmean(temp_Osborn_flux)
     
            arith_mean_dissipation[index] = np.log10(np.nanmean(dissipation_list[index]))

                        
                           
        longitude_list = np.asarray(longitude_list)
         
        left_index_box_edge = np.argmin(np.abs(longitude_list - edge_longitude_interval[0]))
        right_index_box_edge = np.argmin(np.abs(longitude_list - edge_longitude_interval[1]))
        
        edge_Shih_flux = np.nanmean(mean_Shih_flux[left_index_box_edge:right_index_box_edge])          
        edge_Osborn_flux = np.nanmean(mean_Osborn_flux[left_index_box_edge:right_index_box_edge])  

        left_index_box_interior = np.argmin(np.abs(longitude_list - interior_longitude_interval[0]))
        right_index_box_interior = np.argmin(np.abs(longitude_list - interior_longitude_interval[1]))
        interior_Shih_flux = np.nanmean(mean_Shih_flux[left_index_box_interior:right_index_box_interior])           
        interior_Osborn_flux = np.nanmean(mean_Osborn_flux[left_index_box_interior:right_index_box_interior])
 
        edge_Shih_flux_for_different_intervals.append(edge_Shih_flux)
        edge_Osborn_flux_for_different_intervals.append(edge_Osborn_flux)
           
        interior_Shih_flux_for_different_intervals.append(interior_Shih_flux)
        interior_Osborn_flux_for_different_intervals.append(interior_Osborn_flux)

        #print(density_box_width,interior_Shih_flux)

        
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
    
    box_axis[0].plot(np.arange(0.1,4,0.2),np.abs(edge_Osborn_flux_for_different_intervals),"k")
    box_axis[0].plot(np.arange(0.1,4,0.2),np.abs(interior_Osborn_flux_for_different_intervals),"r--")

    box_axis[1].plot(np.arange(0.1,4,0.2),np.abs(edge_Shih_flux_for_different_intervals),"k")
    box_axis[1].plot(np.arange(0.1,4,0.2),np.abs(interior_Shih_flux_for_different_intervals),"r--")
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

    assert np.all(transect_halocline_depths[0:total_number_of_valid_profiles] == transect_halocline_depths[total_number_of_valid_profiles:2*total_number_of_valid_profiles])

    iso_axis[0].plot(longitude_list,transect_halocline_depths[0:total_number_of_valid_profiles])
    iso_axis[1].plot(longitude_list,transect_halocline_densities[0:total_number_of_valid_profiles])
    iso_axis[0].plot(longitude_list,bathymetry_list,"k")



iso_axis[0].invert_yaxis()
box_m0 = box_axis[0].twiny()
#box_m1 = box_axis[1].twiny()
locs = box_axis[0].get_xticks()
box_m0.set_xticks(locs)
box_m0.set_xticklabels(box_sizes_in_dbar)
#box_m1.set_xticks(locs)
#box_m1.set_xticklabels(box_sizes_in_dbar)
plt.show()

    
    
    
    


