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

import scipy.stats as ss 
import geopy.distance as geo
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker

import warnings
warnings.filterwarnings('ignore')

    
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177"] #"/home/ole/windows/processed_mss/emb217"]#,"/home/ole/windows/processed_mss/emb169",
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb217"]
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177"]

rolling_window_size = 9 # for longitudinal averaging
density_box_width = 1 #in kg/m³ (for vertical averaging)

height_above_ground = 20 #Size of the averaging interval above ground for the BBL, has no meaning if (flux_through_halocline == True)
maximum_reasonable_flux = float('Inf') #200 #Fluxes with absolute values above this cut off value will be discarded
maximum_halocline_thickness = 20 #float('Inf') #30

#density_bin_edges = np.linspace(1004,1010.5,20)
density_step = 0.05
density_bin_edges = np.arange(1004,1011,density_step)
density_bin_center = density_bin_edges[:-1] + density_step/2 
#assert len(density_bin_center) == len(density_bin_edges) - 1
number_of_density_bins = density_bin_edges.size -1 
#density_bin_center = np.arange(1004,1010.5,0.2)


list_of_bad_profiles,reasons = np.loadtxt("./data/list_of_bad_profiles.txt", dtype=str, unpack=True)

total_bathymetry_list = []
total_bathymetry_longitude_list = []

f_osborn,osborn_axarr = plt.subplots(nrows = 3, sharex = True)
f_shih,shih_axarr = plt.subplots(nrows = 3, sharex = True)
             
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
    halocline_density_list = []
    halocline_bin_index_list = []     
       
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

            #sort density profiles
            #eps_density_grid = thesis.sort_2D_array_with_nans(eps_density_grid)
            #eps_pot_density_grid = thesis.sort_2D_array_with_nans(eps_pot_density_grid)

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
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_pot_density_grid,eps_oxygen_grid,height_above_ground = height_above_ground)
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
        
        oxygen_gradient_grid = thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth) #in units of micromol/(m*l)
        unit_conversion_grid = 86400 #to convert from m*micromol/(l*s) to mmol/(m^2*d)  ; l to m^3 and micro to milli cancel each other out (factor of 1000)
    
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
            
            #outsider with really high fluxes
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
            halocline_depth,halocline_density,halocline_index = thesis.get_halocline_and_halocline_density(eps_pressure,eps_oxygen_sat_grid[profile],eps_salinity_grid[profile],eps_consv_temperature_grid[profile],eps_pot_density_grid[profile])
            

            #find the corresponding density bin for the halocline_density
            if not np.isnan(halocline_density):
                assert halocline_density < max(density_bin_center) and halocline_density > min(density_bin_center)
                halocline_bin_index =  int(np.nanargmin(np.abs(density_bin_center - halocline_density))) 
            else:
                halocline_bin_index = np.nan 

            
            if not np.isnan(halocline_density):
                #choose the vertical averaging interval dependent on the box size
                
                #used for the raw data, every profile is indepent of its longitudinally surrounding profiles
                from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile,:]) - (halocline_density - (density_box_width/2))))     
                to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile,:]) - (halocline_density + (density_box_width/2))))
                
                
                #check if the vertical interval is bigger than the maximum halocline thickness
                #if yes incrementally adjust the interval to fewer data points
                while True:
                    if halocline_depth - eps_pressure[from_index] <= maximum_halocline_thickness/2:
                        break
                    elif from_index >= halocline_index:
                        break
                    else:
                        from_index += 1
                        
                        
                while True:
                    if eps_pressure[to_index] - halocline_depth <= maximum_halocline_thickness/2:
                        break
                    elif to_index <= halocline_index:
                        break
                    else:
                        to_index -= 1
                
                #print(from_index,to_index, np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile,:]) - halocline_density)))

                assert from_index < to_index   
                
                """
                try:
                    if profile == 35: 
                        raise AssertionError   
                                               
                except AssertionError:
                    print(profile,lon[profile]) #np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile,:]) - (halocline_density - (density_box_width/2)))),
                    plt.plot(eps_pot_density_grid[profile,:],-eps_pressure)
                    plt.axhline(-eps_pressure[np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile,:]) - halocline_density))], c= "k")
                    plt.axhline(-eps_pressure[np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile,:]) - (halocline_density - (density_box_width/2))))], c = "red")
                    plt.axhline(-eps_pressure[np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile,:]) - (halocline_density + (density_box_width/2))))], c= "blue")
                    
                    plt.plot(eps_pot_density_grid[profile+1,:],-eps_pressure,":")
                    plt.plot(eps_pot_density_grid[profile-15,:],-eps_pressure,"--")
                    #plt.axhline(-eps_pressure[np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile+1,:]) - halocline_density))], c= "k")
                    #plt.axhline(-eps_pressure[np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile+1,:]) - (halocline_density - (density_box_width/2))))], c = "red")
                    #plt.axhline(-eps_pressure[np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile+1,:]) - (halocline_density + (density_box_width/2))))], c= "blue")
                    
                    plt.show()
                    raise
                """
                    
                #used for the isopcynal averaging of the whole profiles
                start_density_interval = halocline_bin_index - int((density_box_width/density_step)//2)   #int(np.argmin(np.abs(density_bin_center -(halocline_density - density_box_width/2))))
                stop_density_interval = halocline_bin_index  + int((density_box_width/density_step)//2)+1 #int(np.argmin(np.abs(density_bin_center -(halocline_density + density_box_width/2))))
                
                   
            else:
                from_index, to_index = (0,0)
                start_density_interval = np.nan
                stop_density_interval = np.nan

            
            
            
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
            
            #TODO replace with binned_statistic            
            
            ###############################################################################################################################
            
            #remove all nans in
            #thesis.remove_nans_from_2_arrays(a,b):
            
            #transform data from pressure coordinates to density coordinates
            try:
                density_input, Osborn_flux_input = thesis.remove_nans_from_2_arrays(eps_pot_density_grid[profile],oxygen_flux_Osborn_grid[profile])
                iso_Osborn_flux, _bin_edges, _bin_number = ss.binned_statistic(density_input,Osborn_flux_input, statistic = "mean", bins = density_bin_edges)

                density_input, Shih_flux_input = thesis.remove_nans_from_2_arrays(eps_pot_density_grid[profile],oxygen_flux_Shih_grid[profile])
                iso_Shih_flux, _bin_edges, _bin_number = ss.binned_statistic(density_input,Shih_flux_input, statistic = "mean", bins = density_bin_edges)
                
                density_input, pressure_input = thesis.remove_nans_from_2_arrays(eps_pot_density_grid[profile],eps_pressure)
                iso_pressure, _bin_edges, _bin_number = ss.binned_statistic(density_input,pressure_input, statistic = "mean", bins = density_bin_edges)
                
                density_input, eps_input = thesis.remove_nans_from_2_arrays(eps_pot_density_grid[profile],eps_grid[profile])
                iso_dissipation, _bin_edges, _bin_number = ss.binned_statistic(density_input,eps_input, statistic = "mean", bins = density_bin_edges)
                
                density_input = eps_pot_density_grid[profile][~np.isnan(eps_pot_density_grid[profile])]
                iso_density, _bin_edges, _bin_number = ss.binned_statistic(density_input,density_input, statistic = "mean", bins = density_bin_edges)
                
            #Error handling in the case of all NaN profiles
            except ValueError:
                #print(cruise_name,transect_name,profile,np.all(np.isnan(eps_pot_density_grid[profile])),np.all(np.isnan(oxygen_flux_Osborn_grid[profile])))
                if np.all(np.isnan(oxygen_flux_Osborn_grid[profile])):
                    continue
                
                else:
                    iso_Osborn_flux = np.nan*np.ones(number_of_density_bins)
                    iso_Shih_flux = np.nan*np.ones(number_of_density_bins)
                    iso_pressure = np.nan*np.ones(number_of_density_bins)
                    iso_dissipation = np.nan*np.ones(number_of_density_bins)
                    iso_density = np.nan*np.ones(number_of_density_bins)
            
            #print(len(density_bin_edges),np.shape(iso_Osborn_flux))
                 
                         
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
            
                #density coordinates
                iso_pressure_list.append(iso_pressure)
                iso_dissipation_list.append(iso_dissipation)
                iso_Shih_flux_list.append(iso_Shih_flux)
                iso_Osborn_flux_list.append(iso_Osborn_flux)
                iso_interval_density_list.append([start_density_interval,stop_density_interval])
                
                #pressure coordinates
                dissipation_list.append(eps_grid[profile,from_index:to_index])
                Shih_flux_list.append(oxygen_flux_Shih_grid[profile,from_index:to_index])
                Osborn_flux_list.append(oxygen_flux_Osborn_grid[profile,from_index:to_index])
                
                if from_index != to_index:
                    interval_pressure_list.append([eps_pressure[from_index],eps_pressure[to_index]])
                else:
                    interval_pressure_list.append([np.nan,np.nan])
                    
                longitude_list.append(lon[profile])
                halocline_position_list.append(halocline_depth)
                halocline_density_list.append(halocline_density)
                bathymetry_list.append(bathymetry[profile])
                halocline_bin_index_list.append(halocline_bin_index)
                
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
                
                """
                if np.nanmean(oxygen_flux_Osborn_grid[profile,from_index:to_index]) < - 350:
                    print("#"*50)
                    print(cruise_name,transect_name,profile)
                    print("#"*50)
                """
                    
                if from_index != to_index:
                    interval_pressure_list.insert(list_position,[eps_pressure[from_index],eps_pressure[to_index]])
                else:
                    interval_pressure_list.insert(list_position,[np.nan,np.nan])    
                    
                longitude_list.insert(list_position,lon[profile])
                halocline_position_list.insert(list_position,halocline_depth) 
                halocline_density_list.insert(list_position,halocline_density)
                halocline_bin_index_list.insert(list_position,halocline_bin_index)
                bathymetry_list.insert(list_position,bathymetry[profile])
                    
            assert(np.all(longitude_list == sorted(longitude_list)))

            total_number_of_valid_profiles+=1

            #print(np.shape(iso_Shih_flux_list),number_of_density_bins,total_number_of_valid_profiles)
                
    ###########################################################################################################################
    print(len(longitude_list),"used profiles")
    assert(len(longitude_list) != 0)
    assert(np.all(longitude_list == sorted(longitude_list)))
    assert(np.all(total_bathymetry_longitude_list == sorted(total_bathymetry_longitude_list)))        
    
    
    #from here now on the lists should not be mutated any more
    iso_interval_density_list = np.asarray(iso_interval_density_list)
    iso_dissipation_list = np.asarray(iso_dissipation_list)
    iso_pressure_list = np.asarray(iso_pressure_list)
    iso_Shih_flux_list = np.asarray(iso_Shih_flux_list)
    iso_Osborn_flux_list = np.asarray(iso_Osborn_flux_list)
    
    box_sizes = np.arange(0,3,2*density_step) #in density units [kg/m³]
    #box_sizes = np.arange(0,3.1,0.5) 

    Osborn_flux_density_interval_list = [] 
    Shih_flux_density_interval_list = []
    pressure_density_interval_list = []
    
    up = 0
    down = 1
    for density_box_width in box_sizes: 
    
        transect_Osborn_fluxes_for_current_interval = []
        transect_Shih_fluxes_for_current_interval = []
        pressure_width_for_current_interval = []
        #print(density_box_width,density_box_width//density_step)
        #print(up,down,"\n")
              
        for profile_index in range(np.shape(iso_Shih_flux_list)[0]):
            current_center_bin = halocline_bin_index_list[profile_index]
            #current_center_density = halocline_density_list[profile_index] 
            
            #up = int(np.argmin(np.abs(density_bin_center -(current_center_density - density_box_width/2))))
            #down = int(np.argmin(np.abs(density_bin_center -(current_center_density + density_box_width/2))))
   
   
            if not np.isnan(current_center_bin):
                #vertical mean dependent on the interval
                transect_Osborn_fluxes_for_current_interval.append(np.nanmean(iso_Osborn_flux_list[profile_index,int(current_center_bin - up) : int(current_center_bin + down)]))
                transect_Shih_fluxes_for_current_interval.append(np.nanmean(iso_Shih_flux_list[profile_index,int(current_center_bin - up) : int(current_center_bin + down)]))
                pressure_width_for_current_interval.append(iso_pressure_list[profile_index,int(current_center_bin + down)] - iso_pressure_list[profile_index,int(current_center_bin - up)])

        #mean over all profiles per cruise
        Osborn_flux_density_interval_list.append(np.nanmean(transect_Osborn_fluxes_for_current_interval))
        Shih_flux_density_interval_list.append(np.nanmean(transect_Shih_fluxes_for_current_interval))
        pressure_density_interval_list.append(np.nanmean(pressure_width_for_current_interval))
        
        up +=1
        down +=1
        
    print(pressure_density_interval_list)    
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

    osborn_axarr[figure_index].plot(box_sizes,np.abs(Osborn_flux_density_interval_list), c = color, label = label_name)
    osborn_axarr[figure_index].axvline(0.5,c = "k", ls = ":", alpha = 0.8)
    osborn_paxis = osborn_axarr[figure_index].twinx()
    osborn_paxis.plot(box_sizes,sorted(pressure_density_interval_list), "--", c = color, alpha = 0.7)
    osborn_paxis.set_ylabel(r"$\Delta p$ [dbar]")
    
    shih_axarr[figure_index].plot(box_sizes,np.abs(Shih_flux_density_interval_list), c = color, label = label_name)
    shih_axarr[figure_index].axvline(0.5,c = "k", ls = ":", alpha = 0.8)
    shih_paxis = shih_axarr[figure_index].twinx()
    shih_paxis.set_ylabel(r"$\Delta p$ [dbar]")
    shih_paxis.plot(box_sizes,sorted(pressure_density_interval_list), "--", c = color, alpha = 0.7)
       
    osborn_axarr[figure_index].legend()
    shih_axarr[figure_index].legend()

osborn_axarr[1].set_ylabel(r"|$\langle$ Osborn Oxygen flux $\rangle$| [mmol/m²/d]")
osborn_axarr[2].set_xlabel(r"$\Delta \sigma$ [kg/m³]")

shih_axarr[1].set_ylabel(r"|$\langle$ Shih Oxygen flux $\rangle$| [mmol/m²/d]")
shih_axarr[2].set_xlabel(r"$\Delta \sigma$ [kg/m³]")

f_osborn.suptitle("Impact of the vertical averaging interval")
f_shih.suptitle("Impact of the vertical averaging interval")

width = 6.2012
height = width * 0.7

f_osborn.set_size_inches(width,height)
f_osborn.subplots_adjust(top=0.93,bottom=0.117,left=0.130,right=0.899,hspace=0.2,wspace=0.2) #top=0.827,bottom=0.136,left=0.114,right=0.946,hspace=0.7,wspace=0.2)
f_osborn.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/vertical_interval_Osborn.png", dpi = 400)

f_shih.set_size_inches(width,height)
f_shih.subplots_adjust(top=0.93,bottom=0.117,left=0.104,right=0.899,hspace=0.2,wspace=0.2) #top=0.827,bottom=0.136,left=0.114,right=0.946,hspace=0.7,wspace=0.2)
f_shih.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/vertical_interval_Shih.png", dpi = 400)
    
plt.show()
