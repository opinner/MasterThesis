############################################################
#TODO
##############################################################
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
#matplotlib preferences:
MINI_SIZE = 9
SMALL_SIZE = 10.95
MEDIUM_SIZE = 12
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE, titleweight = "bold")     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MINI_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight = "bold")  # fontsize of the figure title

import geopy.distance as geo
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import scipy.stats as ss
import warnings
#warnings.filterwarnings('ignore')

    
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177"] #"/home/ole/windows/processed_mss/emb217"]#,"/home/ole/windows/processed_mss/emb169",
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb217"]
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177"]

rolling_window_size = 9 # for longitudinal averaging
density_box_width = 1.5 #in kg/m³ (for vertical averaging)

height_above_ground = 20 #Size of the averaging interval above ground for the BBL, has no meaning if (flux_through_halocline == True)
maximum_reasonable_flux = float('Inf') #200 #Fluxes with absolute values above this cut off value will be discarded
maximum_halocline_thickness = 20 #float('Inf') #30

#density_bin_edges = np.linspace(1004,1010.5,20)
density_step = 0.1
density_bin_edges = np.arange(1004,1011,density_step)
#density_bin_center = density_bin_edges[:-1] + density_step/2 
#assert len(density_bin_center) == len(density_bin_edges) - 1
number_of_density_bins = density_bin_edges.size #-1 
#density_bin_center = np.arange(1004,1010.5,0.2)


f_reb,reb_axarr = plt.subplots(nrows = 3, sharex = True, sharey = True)
f_eps,eps_axarr = plt.subplots(nrows = 3, sharex = True, sharey = True)
f_sum,sum_axarr = plt.subplots(nrows = 3, sharex = True)
f_cum,cum_axarr = plt.subplots(nrows = 3, ncols =2, sharex = True, sharey = True)
f_mean,mean_axarr = plt.subplots(nrows = 3, ncols= 2, sharex = True)
f_median,median_axarr = plt.subplots(nrows = 3, ncols= 2, sharex = True)
f_log_sum,log_sum_axarr = plt.subplots(nrows = 1, sharex = True)

list_of_bad_profiles,reasons = np.loadtxt("./data/list_of_bad_profiles.txt", dtype=str, unpack=True)

total_bathymetry_list = []
total_bathymetry_longitude_list = []
             
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    number_of_fluxes_over_the_threshold = 0
    number_of_transects = 0
    total_number_of_fluxes = 0
    number_of_zero_flux = 0
    amount_of_missing_values = 0
    total_number_of_valid_profiles = 0    
    
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruise_name = splitted_foldername[-1]
    
    print(cruise_name)
    
    """
    if cruise_name == "emb217":
        upper_bound_halocline_as_density = 1006.4 #1005.75
        lower_bound_halocline_as_density = 1008.5 #1006.25
    elif cruise_name == "emb177":
        upper_bound_halocline_as_density = 1006.9 #1006.9
        lower_bound_halocline_as_density = 1008.2 #1007.9   
    elif cruise_name == "emb169":
        upper_bound_halocline_as_density = 1006.5 
        lower_bound_halocline_as_density = 1008.6    
    """
    
    #define empty arrays to save the results in
    iso_pressure_list = []
    iso_dissipation_list = []
    #iso_BB_flux_list = []
    iso_Shih_flux_list = []
    iso_Osborn_flux_list = []
    longitude_list = []
    iso_interval_density_list = []
    iso_interval_pressure_list = []

    dissipation_list = []
    BB_flux_list = []
    Shih_flux_list = []
    Osborn_flux_list = []
    interval_pressure_list = []
     
    bathymetry_list = []
    halocline_position_list = []    
     
    hist_Reb_list = []
    hist_Osborn_flux_list = [] 
    hist_Shih_flux_list = []
    hist_dissipation_list = []
           
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".npz":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    #print(DATAFILENAMES)
    
    
    for transect_name in DATAFILENAMES:
    
        datafile_path = FOLDERNAME+"/"+transect_name
        
        transect_name = transect_name[:-4]
    
        if "_".join((cruise_name,transect_name)) in list_of_bad_profiles:
            print("_".join((cruise_name,transect_name)),"skipped")
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
            assert np.all(eps_grid[~np.isnan(eps_grid)] > 0)
            corrected_eps_grid = data["corrected_bin_eps_grid"]
            eps_consv_temperature_grid = data["bin_consv_temperature_grid"]
            eps_salinity_grid = data["bin_salinity_grid"]
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
        BBL,list_of_BBL_indices = results[1] #not needed here
        #BBL_range,list_of_BBL_range_indices = results[2]
        
        #print(np.asarray(list_of_bathymetry_indices) - np.asarray(list_of_BBL_indices))
        
        #eps_N_grid = np.sqrt(eps_N_squared_grid)
        #ozmidov scale
        #ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))
        
        #conversion from pressure coordinates to depth
        eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
        bathymetry_in_m = gsw.z_from_p(bathymetry,np.mean(lat))
        
        eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
        
        distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetry_in_m,(-1,1))
        distance_from_ground_grid[distance_from_ground_grid < 0] = np.nan
        #boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)
        
        turbulent_diffusivity_Osborn_grid = thesis.get_turbulent_diffusivity_Osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
        turbulent_diffusivity_Shih_grid = thesis.get_turbulent_diffusivity_Shih(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
        #turbulent_diffusivity_BB_grid = thesis.get_turbulent_diffusivity_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
          
                     
        shear_velocity_grid = thesis.get_shear_velocity(eps_grid,distance_from_ground_grid)
        
        #compute the mean shear velocity in the BBL as theoretically it should be constant
        shear_velocity = np.zeros(number_of_profiles)
        for profile in range(number_of_profiles):
            BBL_from = int(list_of_BBL_indices[profile])
            BBL_to = int(list_of_bathymetry_indices[profile])
        
            #print(list_of_BBL_indices[profile],list_of_bathymetry_indices[profile])
            BBL_shear_velocity = np.nanmean(shear_velocity_grid[BBL_from:BBL_to])
            
            law_of_the_wall_turbulent_diffusivity = thesis.von_Karman_constant * BBL_shear_velocity * distance_from_ground_grid[profile,BBL_from:BBL_to]
            
            #replace the corresponding bins with the turbulent diffusivity from the law of the wall
            turbulent_diffusivity_Osborn_grid[profile,BBL_from:BBL_to] = law_of_the_wall_turbulent_diffusivity
            turbulent_diffusivity_Shih_grid[profile,BBL_from:BBL_to] = law_of_the_wall_turbulent_diffusivity
        
        oxygen_gradient_grid = thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        unit_conversion_grid = 86400*(1000/eps_density_grid) #to convert from m*micromol/(kg*s) to mmol/(m^2*d)
    
        oxygen_flux_Osborn_grid = - turbulent_diffusivity_Osborn_grid * oxygen_gradient_grid * unit_conversion_grid
        oxygen_flux_Shih_grid = - turbulent_diffusivity_Shih_grid * oxygen_gradient_grid * unit_conversion_grid
        #oxygen_flux_BB_grid =  - turbulent_diffusivity_BB_grid * oxygen_gradient_grid * unit_conversion_grid

        
        # check if the absolute Osborn fluxes are higher than the absolute Shih fluxes (as they should be, because Gamma_Osborn >= Gamma_Shih)
        assert np.all(np.isnan(oxygen_flux_Osborn_grid) == np.isnan(oxygen_flux_Shih_grid))
        test_Osborn = oxygen_flux_Osborn_grid[~np.isnan(oxygen_flux_Osborn_grid)]
        test_Shih = oxygen_flux_Shih_grid[~np.isnan(oxygen_flux_Shih_grid)]
        assert np.all( np.abs(test_Osborn) >= np.abs(test_Shih))
        
        
        for profile in range(number_of_profiles):
        
            #skip profile is necessary
            if "_".join((cruise_name,transect_name,str(profile))) in list_of_bad_profiles:
                print("_".join((cruise_name,transect_name,str(profile))),"skipped")
                continue
            
            if cruise_name == "emb217" and transect_name == "TR1-8" and profile == 44:
                continue
                
            #################################################
            #all iso arrays are in density space and contain one float for every density bin of the profile 
            #the other arrays are in pressure space and contains only the  data points inside the desired density interval     
            #################################################
             
            #allocate arrays for the density bins   
            iso_density_bins = [ [] for _ in range(number_of_density_bins) ]
            iso_dissipation_density_bins = [ [] for _ in range(number_of_density_bins) ]
            iso_Shih_flux_density_bins = [ [] for _ in range(number_of_density_bins) ]
            iso_Osborn_flux_density_bins = [ [] for _ in range(number_of_density_bins) ]
            iso_pressure_density_bins = [ [] for _ in range(number_of_density_bins) ]
        

        
            #arrays in density space that will contain the mean variable for that density bin
            iso_Osborn_flux = []
            iso_Shih_flux = []
            iso_pressure = []
            iso_density = []
            iso_dissipation = []
            
            #calculate the halocline depth and density 
            halocline_depth,halocline_density = thesis.get_halocline_and_halocline_density(eps_pressure,eps_oxygen_sat_grid[profile],eps_salinity_grid[profile],eps_consv_temperature_grid[profile],eps_pot_density_grid[profile])
            
            if not np.isnan(halocline_density):
                #choose the vertical averaging interval dependent on the box size
                
                #used for the raw data, every profile is indepent of its longitudinally surrounding profiles
                from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile]) - (halocline_density - (density_box_width/2))))     
                to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])- (halocline_density + (density_box_width/2))))
                
                
                #check if the vertical interval is bigger than the maximum halocline thickness
                while True:
                    if abs(eps_pressure[from_index] - halocline_depth) < maximum_halocline_thickness/2:
                        break
                    else:
                        from_index += 1
                        
                while True:
                    if abs(eps_pressure[to_index] - halocline_depth) < maximum_halocline_thickness/2:
                        break
                    else:
                        to_index -= 1
                
                assert from_index < to_index                             
                
                #used for the isopcnal averaging of the whole profiles
                start_density_interval = halocline_density - density_box_width/2
                stop_density_interval = halocline_density + density_box_width/2
                
                #above the halocline
                #from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile]) - (halocline_density - density_box_width))) # - density_box_width/2)))     
                #to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])- (halocline_density)))
                
                #below the halocline
                #from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile]) - ())) # - density_box_width/2)))     
                #to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])- (halocline_density + density_box_width)))
                #start_density_interval = halocline_density
                #stop_density_interval = halocline_density + density_box_width          
                
                   
            else:
                from_index, to_index = (0,0)
                start_density_interval = np.nan
                stop_density_interval = np.nan

            
            
            #if from_index == to_index, only an empty list gets appended
            hist_Reb_list.extend(eps_Reynolds_bouyancy_grid[profile,from_index:to_index])
            #TODO only valid for a infinite cutoff value
            hist_Osborn_flux_list.extend(oxygen_flux_Osborn_grid[profile,from_index:to_index])
            hist_Shih_flux_list.extend(oxygen_flux_Shih_grid[profile,from_index:to_index])
            hist_dissipation_list.extend(eps_grid[profile,from_index:to_index])
            
            """
            #remove the last 3 indexes above the sea floor
            sea_floor_index =  -np.argmax(np.flip(~np.isnan(eps_grid[profile,:]))) #at the moment the index is negative
            sea_floor_index = eps_grid[profile,:].size + sea_floor_index #now defined as positive index
            
            if sea_floor_index == eps_grid[profile,:].size:
                if not np.isnan(eps_grid[profile,-1]): #if there are no NAN values towards the bottom
                    sea_floor_index = len(eps_pressure)-1 #set the last index as the index of the bottom
        
            assert sea_floor_index == bathymetry[profile]
            #sea_floor_pressure = eps_pressure[sea_floor_index]
            
            #to_index = max(from_index,min(sea_floor_index-3,to_index))
            """
            
            """
            #sort the profile into the density bins
            for value_index,value in enumerate(eps_pot_density_grid[profile]):
                    
                bin_number = np.nan
                
                #find the corresponding density bin
                for bin_index,density_bin_edge in enumerate(density_bin_edges):  
                    #print(value,density_bin_edge)      
                    if value <= density_bin_edge:
                        bin_number = bin_index
                        break            
            
                    #Test for edge cases?
                    if value > density_bin_edges[-1]:
                        raise AssertionError
            
                # add every point from that profile to the corresponding density bin
                if not np.isnan(bin_number): 
                        iso_Shih_flux_density_bins[bin_number].append(oxygen_flux_Shih_grid[profile,value_index])
                        iso_Osborn_flux_density_bins[bin_number].append(oxygen_flux_Osborn_grid[profile,value_index])
                        iso_pressure_density_bins[bin_number].append(eps_pressure[value_index])
                        iso_density_bins[bin_number].append(eps_pot_density_grid[profile,value_index])
                        iso_dissipation_density_bins[bin_number].append(eps_grid[profile,value_index])
                        
            #average the values inside the density bins
            for bin_index in range(number_of_density_bins):
                #print(np.size(density_bins[bin_index]))
                if len(iso_Shih_flux_density_bins[bin_index]) == 0: #test for empty list, ergo no points in that density bin
                    #continue
                    iso_dissipation.append(np.nan)
                    iso_Shih_flux.append(np.nan)
                    iso_Osborn_flux.append(np.nan)
                    iso_pressure.append(np.nan)
                    iso_density.append(np.nan)                
                    
                else:
                    
                    temp_iso_Osborn_flux = np.asarray(iso_Osborn_flux_density_bins[bin_index])
                    temp_iso_Osborn_flux[np.abs(temp_iso_Osborn_flux)>maximum_reasonable_flux] = np.nan
                    
                    temp_iso_Shih_flux = np.asarray(iso_Shih_flux_density_bins[bin_index])
                    temp_iso_Shih_flux[np.abs(temp_iso_Shih_flux)>maximum_reasonable_flux] = np.nan
    
                    iso_Osborn_flux.append(np.nanmean(temp_iso_Osborn_flux))
                    iso_Shih_flux.append(np.nanmean(temp_iso_Shih_flux))
                    iso_pressure.append(np.nanmean(iso_pressure_density_bins[bin_index]))
                    iso_density.append(np.mean(iso_density_bins[bin_index]))
                    iso_dissipation.append(np.nanmean(iso_dissipation_density_bins[bin_index]))
                    
                    
            

                        
            #bathymetry: find the correct position in the sorted longitude list
            for index,value in enumerate(total_bathymetry_longitude_list):
                if value > lon[profile]:
                    list_position = index
                    break
                elif index == len(total_bathymetry_longitude_list)-1:
                    list_position = len(total_bathymetry_longitude_list)
                    break
                                
            if len(total_bathymetry_longitude_list) == 0:   
                total_bathymetry_list.append(bathymetry[profile])
                total_bathymetry_longitude_list.append(lon[profile])
                

            else:
                total_bathymetry_list.insert(list_position,bathymetry[profile])
                total_bathymetry_longitude_list.insert(list_position,lon[profile])

            
                                            
            #fluxes: find the correct position in the sorted longitude list
            for index,value in enumerate(longitude_list):
                if value > lon[profile]:
                    list_position = index
                    break
                elif index == len(longitude_list)-1:
                    list_position = len(longitude_list)
                    break
            
                 
            if len(longitude_list) == 0: 
                iso_pressure_list.append(iso_pressure)
                iso_dissipation_list.append(iso_dissipation)
                iso_Shih_flux_list.append(iso_Shih_flux)
                iso_Osborn_flux_list.append(iso_Osborn_flux)
                iso_interval_density_list.append([start_density_interval,stop_density_interval])
                
                dissipation_list.append(eps_grid[profile,from_index:to_index])
                Shih_flux_list.append(oxygen_flux_Shih_grid[profile,from_index:to_index])
                Osborn_flux_list.append(oxygen_flux_Osborn_grid[profile,from_index:to_index])
                
                if from_index != to_index:
                    interval_pressure_list.append([eps_pressure[from_index],eps_pressure[to_index]])
                else:
                    interval_pressure_list.append([np.nan,np.nan])
                    
                longitude_list.append(lon[profile])
                halocline_position_list.append(halocline_depth) 
                bathymetry_list.append(bathymetry[profile])
                
            else:
                #Sort the current profile  after their longitude coordinates into the list    
                iso_dissipation_list.insert(list_position,iso_dissipation)
                iso_pressure_list.insert(list_position,iso_pressure)
                iso_Shih_flux_list.insert(list_position,iso_Shih_flux)
                iso_Osborn_flux_list.insert(list_position,iso_Osborn_flux)
                iso_interval_density_list.insert(list_position,[start_density_interval,stop_density_interval])
                          
                dissipation_list.insert(list_position,eps_grid[profile,from_index:to_index])
                #BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.insert(list_position,oxygen_flux_Shih_grid[profile,from_index:to_index])
                Osborn_flux_list.insert(list_position,oxygen_flux_Osborn_grid[profile,from_index:to_index])
                
                
                #if np.nanmean(oxygen_flux_Osborn_grid[profile,from_index:to_index]) < - 350:
                #    print("#"*50)
                #    print(cruise_name,transect_name,profile)
                #    print("#"*50)
                
                    
                if from_index != to_index:
                    interval_pressure_list.insert(list_position,[eps_pressure[from_index],eps_pressure[to_index]])
                else:
                    interval_pressure_list.insert(list_position,[np.nan,np.nan])    
                    
                longitude_list.insert(list_position,lon[profile])
                halocline_position_list.insert(list_position,halocline_depth) 
                bathymetry_list.insert(list_position,bathymetry[profile])
                    
            assert(np.all(longitude_list == sorted(longitude_list)))

            total_number_of_valid_profiles+=1

            #print(np.shape(iso_Shih_flux_list),number_of_density_bins,total_number_of_valid_profiles)
            """
                
    ###########################################################################################################################
    #print(len(longitude_list),"used profiles")
    #assert(len(longitude_list) != 0)
    #assert(np.all(longitude_list == sorted(longitude_list)))
    #assert(np.all(total_bathymetry_longitude_list == sorted(total_bathymetry_longitude_list)))        
    
    """
    #from here now on the lists should not be mutated any more
    iso_interval_density_list = np.asarray(iso_interval_density_list)
    iso_dissipation_list = np.asarray(iso_dissipation_list)
    iso_pressure_list = np.asarray(iso_pressure_list)
    iso_Shih_flux_list = np.asarray(iso_Shih_flux_list)
    iso_Osborn_flux_list = np.asarray(iso_Osborn_flux_list)
              
    dissipation_list = np.asarray(dissipation_list)
    Shih_flux_list = np.asarray(Shih_flux_list)
    Osborn_flux_list = np.asarray(Osborn_flux_list)
    interval_pressure_list = np.asarray(interval_pressure_list)
    
    #compute mean and std over the saved intervals (that will not be isopycnical averaged)
    mean_Osborn_flux = [None] * total_number_of_valid_profiles
    mean_Shih_flux = [None] * total_number_of_valid_profiles
    arith_mean_dissipation = [None] * total_number_of_valid_profiles

    
    #remove fluxes over the cut off value
    for index in range(total_number_of_valid_profiles):
        temp_Shih_flux = Shih_flux_list[index]
        number_of_fluxes_over_the_threshold += np.sum(np.abs(temp_Shih_flux)>maximum_reasonable_flux)
        number_of_zero_flux += np.sum(np.abs(temp_Shih_flux)==0)
        amount_of_missing_values += np.sum(np.isnan(temp_Shih_flux))
        #count the number of flux data points        
        total_number_of_fluxes += temp_Shih_flux.size
        
        temp_Shih_flux[np.abs(temp_Shih_flux)>maximum_reasonable_flux] = np.nan
        mean_Shih_flux[index] = np.nanmean(temp_Shih_flux)
        
        temp_Osborn_flux = Osborn_flux_list[index]
        temp_Osborn_flux[np.abs(temp_Osborn_flux)>maximum_reasonable_flux] = np.nan  
        mean_Osborn_flux[index] = np.nanmean(temp_Osborn_flux)
        
        arith_mean_dissipation[index] = np.log10(np.nanmean(dissipation_list[index]))
        
    #remove fluxes over the cut off value in the arrays for isopycnal averaging

                       
                           
    iso_rolling_mean_Shih_flux = [None] * total_number_of_valid_profiles
    iso_rolling_mean_Osborn_flux = [None] * total_number_of_valid_profiles
    iso_rolling_mean_dissipation = [None] * total_number_of_valid_profiles
    iso_mean_pressure = [None] * total_number_of_valid_profiles
    
    
    #compute rolling longitudinal average
    for index in range(total_number_of_valid_profiles):

        left = index-(rolling_window_size//2)
        right = index+rolling_window_size//2

        #controls that the mean is not computed over too distant points
        number_of_nans_in_averaging_window = np.count_nonzero(np.isnan(mean_Shih_flux[left:right])) 
        
        if number_of_nans_in_averaging_window > 0.5 * rolling_window_size or right >= total_number_of_valid_profiles or left <= 0:
 
            iso_mean_pressure[index] = np.nan*np.ones(number_of_density_bins)
            iso_rolling_mean_Shih_flux[index] = np.nan*np.ones(number_of_density_bins)
            iso_rolling_mean_Osborn_flux[index] = np.nan*np.ones(number_of_density_bins)
            iso_rolling_mean_dissipation[index] = np.nan*np.ones(number_of_density_bins)
            continue

        if rolling_window_size ==1:
        
            iso_mean_pressure[index] = iso_pressure_list[index,:]
            iso_rolling_mean_Shih_flux[index] = iso_Shih_flux_list[index,:]
            iso_rolling_mean_Osborn_flux[index] = iso_Osborn_flux_list[index,:]
            iso_rolling_mean_dissipation[index] = iso_dissipation_list[index,:]
            
        else:
            try:
                iso_mean_pressure[index] = np.nanmean(iso_pressure_list[left:right,:],axis = 0)
                iso_rolling_mean_Shih_flux[index] = np.nanmean(iso_Shih_flux_list[left:right,:],axis = 0)
                iso_rolling_mean_Osborn_flux[index] = np.nanmean(iso_Osborn_flux_list[left:right,:],axis = 0)
                iso_rolling_mean_dissipation[index] = np.nanmean(iso_dissipation_list[left:right,:],axis = 0)


            except (IndexError,ValueError):
                raise AssertionError
                iso_mean_pressure[index] = np.nan*np.ones(number_of_density_bins)
                iso_rolling_mean_Shih_flux[index] = np.nan*np.ones(number_of_density_bins)
                iso_rolling_mean_Osborn_flux[index] = np.nan*np.ones(number_of_density_bins)
                iso_rolling_arith_mean_dissipation[index] = np.nan*np.ones(number_of_density_bins)
     
    #Assure that the isopycnal averaging did not change the overall shape    
    #for element in iso_mean_pressure:
    #    print(np.shape(element),element) 
    #print(np.shape(iso_mean_pressure),np.shape(iso_pressure_list))
    assert np.shape(iso_mean_pressure) == np.shape(iso_pressure_list)
    assert np.shape(iso_rolling_mean_Shih_flux) == np.shape(iso_Shih_flux_list)
    
    #convert the arrays from lists to numpy arrays
    iso_rolling_mean_Shih_flux = np.asarray(iso_rolling_mean_Shih_flux)
    iso_rolling_mean_Osborn_flux = np.asarray(iso_rolling_mean_Osborn_flux)
    iso_rolling_mean_dissipation = np.asarray(iso_rolling_mean_dissipation)
    iso_mean_pressure = np.asarray(iso_mean_pressure)
    
    #allocate new arrays for the halocline results
    iso_vertical_mean_Shih_flux = np.asarray([None] * total_number_of_valid_profiles)
    iso_vertical_mean_Osborn_flux = np.asarray([None] * total_number_of_valid_profiles)
    iso_vertical_mean_dissipation = np.asarray([None] * total_number_of_valid_profiles)
     
    #extract the desired vertical interval  
    for profile in range(total_number_of_valid_profiles):      
       
        #look up the vertical density interval for this profile
        interval_start, interval_stop = iso_interval_density_list[profile]

        #check if the interval are indeed numbers
        if not np.isnan(interval_start) and not np.isnan(interval_stop) and not np.isnan(halocline_position_list[profile]):
        
            #compute which bins are in that profile
            start_interval_index = np.nanargmin(np.abs(density_bin_edges - interval_start))
            stop_interval_index = np.nanargmin(np.abs(density_bin_edges - interval_stop))
            
            try:
                while True:
                    #if the mean position of that density bin is too far away from the halocline, start with the next density bin and check again
                    if abs(iso_mean_pressure[profile,start_interval_index] - halocline_position_list[profile]) < maximum_halocline_thickness/2:
                        break
                    else:
                        start_interval_index += 1
                           
                while True:
                    #if the mean position of that density bin is too far away from the halocline, stop at the previous density bin and check again
                    if abs(iso_mean_pressure[profile,stop_interval_index] - halocline_position_list[profile]) < maximum_halocline_thickness/2:
                        break
                    else:
                        stop_interval_index -= 1
            
                #the interval should still start higher than it stops    
                #print(start_interval_index,stop_interval_index,iso_mean_pressure[profile,start_interval_index],iso_mean_pressure[profile,stop_interval_index])
                assert start_interval_index < stop_interval_index
                assert iso_mean_pressure[profile,start_interval_index] < iso_mean_pressure[profile,stop_interval_index]
                #average over the bins inside that interval 
                iso_vertical_mean_Shih_flux[profile] = np.nanmean(iso_rolling_mean_Shih_flux[profile,start_interval_index:stop_interval_index+1]) #TODO Why the plus 1???
                iso_vertical_mean_Osborn_flux[profile] = np.nanmean(iso_rolling_mean_Osborn_flux[profile,start_interval_index:stop_interval_index+1])
                iso_vertical_mean_dissipation[profile] = np.log10(np.nanmean(iso_rolling_mean_dissipation[profile,start_interval_index:stop_interval_index+1]))

                iso_interval_pressure_list.append([iso_mean_pressure[profile,start_interval_index],iso_mean_pressure[profile,stop_interval_index+1]])

            
            except (IndexError,AssertionError):
                iso_vertical_mean_Shih_flux[profile] = np.nan
                iso_vertical_mean_Osborn_flux[profile] = np.nan
                iso_vertical_mean_dissipation[profile] = np.nan
                iso_interval_pressure_list.append([np.nan,np.nan])                 


        else:
            iso_vertical_mean_Shih_flux[profile] = np.nan
            iso_vertical_mean_Osborn_flux[profile] = np.nan
            iso_vertical_mean_dissipation[profile] = np.nan
            iso_interval_pressure_list.append([np.nan,np.nan])  
    
    
    iso_interval_pressure_list = np.asarray(iso_interval_pressure_list)
    assert np.shape(iso_interval_pressure_list) == np.shape(interval_pressure_list)
    
    print(cruise_name,"total number of transects =",number_of_transects)        
    print("total_number_of_valid_profiles",total_number_of_valid_profiles)     
    #print(np.shape(iso_vertical_mean_Shih_flux),np.shape(mean_Shih_flux))
    """
    
    """
    print("number_of_fluxes_over_the_threshold\ttotal_number_of_fluxes\tratio")
    print("NaN",amount_of_missing_values,total_number_of_fluxes,100*amount_of_missing_values/total_number_of_fluxes,"%")
    print("0",number_of_zero_flux,total_number_of_fluxes,100*number_of_zero_flux/total_number_of_fluxes,"%")
    print(">",number_of_fluxes_over_the_threshold,total_number_of_fluxes,100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    print("Sum:",100*amount_of_missing_values/total_number_of_fluxes + 100*number_of_zero_flux/total_number_of_fluxes + 100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    """

    """
    distance_list = np.zeros(np.shape(longitude_list))
    first_point = (np.mean(lat),longitude_list[0]) 
    for i in range(len(longitude_list)):
        current_point = (np.mean(lat),longitude_list[i])
        distance_list[i] = geo.geodesic(first_point,current_point).km #Distance in km
    """
    
        
    """
    delta_X = thesis.central_differences(distance_list)
    
    print("\n\n\n",cruise_name,"flux sum:")
    print("Osborn rolling mean",np.nansum(rolling_mean_Osborn_flux*delta_X),"Osborn profiles", np.nansum(mean_Osborn_flux*delta_X))
    print("Shih rolling mean",np.nansum(rolling_mean_Shih_flux*delta_X),"Shih profiles", np.nansum(mean_Shih_flux*delta_X))
    print("\n\n\n")
    """

    
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
    
    
    #the colorscheme ['#d95f02','#7570b3','#1b9e77'] stems from colorbrewer (colorbrewer2.org) to be friendly to color blindness and colored printing
    for figure_index,color,label_name,cruise in zip([2,1,0],['#d95f02','#7570b3','#1b9e77'],["summer cruise emb217","winter cruise emb177","autumn cruise emb169"],["emb217","emb177","emb169"]):
        if cruise_name == cruise:
            break
 
 
    hist_Reb_list = np.asarray(hist_Reb_list)
    hist_Osborn_flux_list = np.asarray(hist_Osborn_flux_list)
    hist_Shih_flux_list = np.asarray(hist_Shih_flux_list)
    
    print(np.shape(hist_Reb_list))
   
    #remove Nans from data
    data_mask = np.isnan(hist_Reb_list)
    hist_Reb_list = hist_Reb_list[~data_mask]
    hist_Osborn_flux_list = hist_Osborn_flux_list[~data_mask]
    hist_Shih_flux_list = hist_Shih_flux_list[~data_mask]

    data_mask = np.isnan(hist_Osborn_flux_list)
    hist_Reb_list = hist_Reb_list[~data_mask]
    hist_Osborn_flux_list = hist_Osborn_flux_list[~data_mask]
    hist_Shih_flux_list = hist_Shih_flux_list[~data_mask]
        
    assert np.all(~np.isnan(hist_Shih_flux_list))
   
    #create logarithmic equidistant bins
    bins = np.logspace(np.log10(min(hist_Reb_list)),np.log10(max(hist_Reb_list)),50)
    #bins = np.linspace(min(hist_Reb_list),max(hist_Reb_list),100)
    
    reb_hist_values, bin_edges, _binnumber = ss.binned_statistic(hist_Reb_list,hist_Osborn_flux_list, statistic = "count", bins = bins)
    Osborn_hist_values, bin_edges, _binnumber = ss.binned_statistic(hist_Reb_list,hist_Osborn_flux_list, statistic = "sum", bins = bins)
    Shih_hist_values, bin_edges, _binnumber = ss.binned_statistic(hist_Reb_list,hist_Shih_flux_list, statistic = "sum", bins = bins)
    

    
    #take the mean and median inside the bans and set all nans to 0
    mean_Osborn_hist_values, bin_edges, _binnumber = ss.binned_statistic(hist_Reb_list,hist_Osborn_flux_list, statistic = "mean", bins = bins)
    mean_Osborn_hist_values[np.isnan(mean_Osborn_hist_values)] = 0
    mean_Shih_hist_values, bin_edges, _binnumber = ss.binned_statistic(hist_Reb_list,hist_Shih_flux_list, statistic = "mean", bins = bins)
    mean_Shih_hist_values[np.isnan(mean_Shih_hist_values)] = 0
    median_Osborn_hist_values, bin_edges, _binnumberX = ss.binned_statistic(hist_Reb_list,hist_Osborn_flux_list, statistic = "median", bins = bins)
    median_Osborn_hist_values[np.isnan(median_Osborn_hist_values)] = 0
    median_Shih_hist_values, bin_edges, _binnumber = ss.binned_statistic(hist_Reb_list,hist_Shih_flux_list, statistic = "median", bins = bins)
    median_Shih_hist_values[np.isnan(median_Shih_hist_values)] = 0
    
    """
    sorted_bin_number = [x for _,x in sorted(zip(median_Osborn_hist_values,_binnumberX))]
    for index,value in enumerate(sorted(median_Osborn_hist_values)):
        if sorted_bin_number[index] == len(median_Osborn_hist_values) -1:
            plt.plot(index,sorted(median_Osborn_hist_values)[index])
    plt.show()
    """
              
    reb_axarr[figure_index].set_title(label_name)
    reb_axarr[figure_index].hist(bin_edges[:-1], bin_edges, weights= reb_hist_values, log = True, histtype='bar', ec='black', color = color)
    cum_axarr[figure_index,0].hist(bin_edges[:-1], bin_edges, weights= -Osborn_hist_values/sum(-Osborn_hist_values), histtype='bar', color = color, alpha = 0.6, fill = True, hatch='xx', cumulative=True)
    cum_axarr[figure_index,1].hist(bin_edges[:-1], bin_edges, weights= -Shih_hist_values/sum(-Shih_hist_values), histtype='bar', color = color, fill = True, edgecolor='k', cumulative=True)
    sum_axarr[figure_index].hist(bin_edges[:-1], bin_edges, weights= -Osborn_hist_values, histtype='bar', color = color, alpha = 0.6, fill = True, hatch='xx')
    sum_axarr[figure_index].hist(bin_edges[:-1], bin_edges, weights= -Shih_hist_values, histtype='bar', color = color, fill = True, edgecolor='k')
    
    #mean_axarr[figure_index,0].plot(bin_edges[:-1],-mean_Osborn_hist_values, color = color)
    #mean_axarr[figure_index,1].plot(bin_edges[:-1],-mean_Shih_hist_values, color = color)
    #median_axarr[figure_index,0].plot(bin_edges[:-1],-median_Osborn_hist_values, color = color)
    #median_axarr[figure_index,1].plot(bin_edges[:-1],-median_Shih_hist_values, color = color)

    #mean_axarr[figure_index,0].bar(bin_edges[:-1],-mean_Osborn_hist_values, color = color, width = np.diff(bin_edges), edgecolor='k')
    #mean_axarr[figure_index,1].bar(bin_edges[:-1],-mean_Shih_hist_values, color = color, width = np.diff(bin_edges), edgecolor='k')
    #median_axarr[figure_index,0].bar(bin_edges[:-1],-median_Osborn_hist_values, color = color, width = np.diff(bin_edges), edgecolor='k')
    #median_axarr[figure_index,1].bar(bin_edges[:-1],-median_Shih_hist_values, color = color, width = np.diff(bin_edges), edgecolor='k')
        
    mean_axarr[figure_index,0].hist(bin_edges[:-1], bin_edges, weights= -mean_Osborn_hist_values, histtype='bar', ec='black', color = color)
    mean_axarr[figure_index,1].hist(bin_edges[:-1], bin_edges, weights= -mean_Shih_hist_values, histtype='bar', ec='black', color = color)
    median_axarr[figure_index,0].hist(bin_edges[:-1], bin_edges, weights= -median_Osborn_hist_values, histtype='bar', ec='black', color = color)
    median_axarr[figure_index,1].hist(bin_edges[:-1], bin_edges, weights= -median_Shih_hist_values, histtype='bar', ec='black', color = color)
        
    if cruise_name == "emb177":
        log_sum_axarr.hist(bin_edges[:-1], bin_edges, weights= -Osborn_hist_values, histtype='bar', color = color, alpha = 0.6, fill = True, hatch='xx', label = "Osborn fluxes")
        log_sum_axarr.hist(bin_edges[:-1], bin_edges, weights= -Shih_hist_values, histtype='bar', color = color, fill = True, edgecolor='k', label = "Shih fluxes")
        log_sum_axarr.set_xscale('log')
        log_sum_axarr.legend(loc = "upper left")   
        log_sum_axarr.set_ylabel("Summed downward oxygen fluxes [mmol/(m$^2$d)]")
        log_sum_axarr.set_xlabel("Reynolds Buoyancy Number")
        
    mean_axarr[figure_index,0].set_title(label_name+" Osborn")
    mean_axarr[figure_index,1].set_title(label_name+" Shih")

    median_axarr[figure_index,0].set_title(label_name+" Osborn")
    median_axarr[figure_index,1].set_title(label_name+" Shih")
        
    reb_axarr[figure_index].set_xscale('log')
    mean_axarr[figure_index,0].set_xscale('log')
    median_axarr[figure_index,0].set_xscale('log')
    sum_axarr[figure_index].set_xscale('log')
    cum_axarr[figure_index,0].set_xscale('log')
    eps_axarr[figure_index].set_xscale('log')
        
    eps_axarr[figure_index].set_title(label_name)
    eps_bins = np.logspace(np.log10(min(hist_dissipation_list)),np.log10(max(hist_dissipation_list)),30)
    eps_axarr[figure_index].hist(hist_dissipation_list, bins = eps_bins ,histtype='bar', ec='black', color = color)

    
    


    
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
width = 6.2012
height = width * 1.618

#beamer figure sizes
width = 1.7*4.252 #6.2012
height = 1.5*3.7341 #* 4/3 #1.618


f_log_sum.set_size_inches(width, height)
f_reb.set_size_inches(width, height)
f_eps.set_size_inches(width, height)
f_sum.set_size_inches(width, height)

emb169_label = mpatches.Patch(color='#1b9e77', label='autumn cruise emb169')
emb177_label = mpatches.Patch(color='#7570b3', label='winter cruise emb177')
emb217_label = mpatches.Patch(color='#d95f02', label='summer cruise emb217')
Osborn_label = mlines.Line2D([], [], ls = ":", lw = 2.5, c = "k", label = "Osborn oxygen flux")
Shih_label = mlines.Line2D([], [], ls = "-", lw = 2.5, c = "k", label = "Shih oxygen flux")

f_reb.set_size_inches(width,height)
f_reb.subplots_adjust(top=0.92,bottom=0.09,left=0.12,right=0.975,hspace=0.058,wspace=0.185)
f_eps.set_size_inches(width,height)
f_eps.subplots_adjust(top=0.92,bottom=0.09,left=0.12,right=0.975,hspace=0.058,wspace=0.185)
f_sum.set_size_inches(width,height)
f_sum.subplots_adjust(top=0.92,bottom=0.09,left=0.12,right=0.975,hspace=0.058,wspace=0.185)
f_log_sum.set_size_inches(width,height)
f_log_sum.tight_layout() #subplots_adjust(top=0.95,bottom=0.09,left=0.12,right=0.975,hspace=0.058,wspace=0.185)   
f_mean.set_size_inches(width,height)
f_mean.subplots_adjust(top=0.92,bottom=0.09,left=0.12,right=0.975,hspace=0.22,wspace=0.185)
f_median.set_size_inches(width,height)
f_median.subplots_adjust(top=0.92,bottom=0.09,left=0.12,right=0.975,hspace=0.22,wspace=0.185)
  
f_reb.suptitle(r"Occurences of $Re_{\mathrm{b}}$")
f_eps.suptitle(r"Occurences of dissipation rate $\epsilon$")
f_sum.suptitle(r"Summed downward oxygen fluxes")
f_median.suptitle(r"Median downward oxygen fluxes")
f_mean.suptitle(r"Mean downward oxygen fluxes")
#f_log_sum.suptitle(r"Summed downward oxygen fluxes of emb177")

#f_reb.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/beamer_Reb_hist.png", dpi = 600)
#f_eps.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/beamer_eps_hist.png", dpi = 600)
#f_sum.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/beamer_log_flux_Reb.png", dpi = 600)
f_log_sum.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/beamer_flux_Reb.png", dpi = 600)
   
plt.show()
    
    
    
    
    
