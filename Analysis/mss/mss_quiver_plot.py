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

    
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb217"]

rolling_window_size = 9 # for longitudinal averaging
density_box_width = 0.5 #in kg/m³ (for vertical averaging)

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

f_169,axarr_169 = plt.subplots(1)
f_177,axarr_177 = plt.subplots(1)
f_217,axarr_217 = plt.subplots(1)
f_all, all_axarr = plt.subplots(3, sharex = True, sharey = True) 
fi, iaxarr = plt.subplots(1)

#textstr = ""

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
    
    
    #define empty arrays to save the results in
    iso_pressure_list = []
    iso_dissipation_list = []
    #iso_BB_flux_list = []
    iso_Shih_flux_list = []
    iso_Osborn_flux_list = []
    longitude_list = []
    iso_interval_density_list = []
    iso_interval_pressure_list = []

    density_list = []
    dissipation_list = []
    BB_flux_list = []
    Shih_flux_list = []
    Osborn_flux_list = []
    interval_pressure_list = []
     
    bathymetry_list = []
    halocline_position_list = []    
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
            
        if cruise_name == "emb217" and transect_name == "TR1-8" and profile == 44:
            continue
            
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
            
        
        bool1 = cruise_name == "emb169" and transect_name == "TS112"
        bool2 = cruise_name == "emb177" and transect_name == "TS1_3"
        bool3 = cruise_name == "emb217" and transect_name == "TR1-4"
        if np.any([bool1,bool2,bool3]): 
            isopycnals_density = eps_pot_density_grid
            isopycnals_lon = lon
            print("done")
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

                assert from_index <= to_index   
                
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
                start_density_interval = int(np.argmin(np.abs(density_bin_center - (halocline_density - density_box_width/2))))
                stop_density_interval = int(np.argmin(np.abs(density_bin_center - (halocline_density + density_box_width/2))))
                
                   
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
                #density_list.append(eps_pot_density_grid[profile])
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
                halocline_bin_index_list.append(halocline_bin_index)
                
            else:
                #Sort the current profile  after their longitude coordinates into the list    
                iso_dissipation_list.insert(list_position,iso_dissipation)
                iso_pressure_list.insert(list_position,iso_pressure)
                iso_Shih_flux_list.insert(list_position,iso_Shih_flux)
                iso_Osborn_flux_list.insert(list_position,iso_Osborn_flux)
                iso_interval_density_list.insert(list_position,[start_density_interval,stop_density_interval])
                     
                #density_list.insert(list_position,eps_pot_density_grid[profile])          
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
    
    
    #print("-"*40)
    #print(np.shape(density_list), np.shape(eps_pressure), np.shape(longitude_list))
    #print("-"*40)
    #plt.plot(np.nanmean(density_list, axis = 0), -eps_pressure)
    
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
                #print(np.shape(iso_pressure_list),np.shape(iso_Shih_flux_list),len(density_bin_edges))
                
                iso_mean_pressure[index] = np.nanmean(iso_pressure_list[left:right,:],axis = 0)
                iso_rolling_mean_Shih_flux[index] = np.nanmean(iso_Shih_flux_list[left:right,:],axis = 0)
                iso_rolling_mean_Osborn_flux[index] = np.nanmean(iso_Osborn_flux_list[left:right,:],axis = 0)
                iso_rolling_mean_dissipation[index] = np.nanmean(iso_dissipation_list[left:right,:],axis = 0)


            except (IndexError,ValueError):
                raise AssertionError('Accessing the array did not work')
                #iso_mean_pressure[index] = np.nan*np.ones(number_of_density_bins)
                #iso_rolling_mean_Shih_flux[index] = np.nan*np.ones(number_of_density_bins)
                #iso_rolling_mean_Osborn_flux[index] = np.nan*np.ones(number_of_density_bins)
                #iso_rolling_arith_mean_dissipation[index] = np.nan*np.ones(number_of_density_bins)
     
    """   
    for pressure_profile in iso_mean_pressure:
        try:
            assert np.all(np.diff(pressure_profile) > 0)
        except AssertionError:    
            print(np.diff(pressure_profile)) 
            raise
            
    #print(np.shape(iso_mean_pressure),np.shape(iso_pressure_list))
    """
    
    #Assure that the isopycnal averaging did not change the overall shape 
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
    iso_vertical_mean_dissipation = np.asarray([None] * total_number_of_valid_profiles)
    iso_vertical_mean_pressure = np.asarray([None] * total_number_of_valid_profiles)
    
    off_set_list = range(0,6)
    quiver_array = np.zeros((len(off_set_list),total_number_of_valid_profiles)) #3 rows of length "total_number_of_valid_profiles"
    quiver_pressure = np.zeros((len(off_set_list),total_number_of_valid_profiles))
    
    #isopycnal_density_list =  np.zeros(len(off_set_list))
    
    """
    for offset in off_set_list:
        
        iso_vertical_mean_Shih_flux = np.asarray([None] * total_number_of_valid_profiles)
        iso_vertical_mean_Osborn_flux = np.asarray([None] * total_number_of_valid_profiles)
        iso_vertical_mean_dissipation = np.asarray([None] * total_number_of_valid_profiles)
        iso_vertical_mean_dissipation = np.asarray([None] * total_number_of_valid_profiles)
        iso_vertical_mean_pressure = np.asarray([None] * total_number_of_valid_profiles)
    
        #extract the desired vertical interval  
        for profile in range(total_number_of_valid_profiles):      
           
           
            #look up the indices of the vertical density interval for this profile
            interval_start, interval_stop = iso_interval_density_list[profile]
            halocline_bin_index = halocline_bin_index_list[profile]
            
            #check if the interval are indeed numbers
            if not np.isnan(interval_start) and not np.isnan(interval_stop) and not np.isnan(halocline_position_list[profile]):
              
                try:
                    
                    start_interval_index,stop_interval_index = int(interval_start),int(interval_stop)
                    #for the original interval
                    if offset == 0:
                    
                        #check if the vertical interval is bigger than the maximum halocline thickness
                        #if yes incrementally adjust the interval to fewer data points
                        while True:
                            if iso_mean_pressure[profile,halocline_bin_index] - iso_mean_pressure[profile,start_interval_index] <= maximum_halocline_thickness/2:
                                break
                            elif start_interval_index >= halocline_bin_index:
                                break
                            else:
                                start_interval_index += 1
                                
                                
                        while True:
                            if iso_mean_pressure[profile,stop_interval_index] - iso_mean_pressure[profile,halocline_bin_index] <= maximum_halocline_thickness/2:
                                break
                            elif stop_interval_index <= halocline_bin_index:
                                break
                            else:
                                stop_interval_index -= 1
                
                    #for all further intervals
                    else:
                        start_interval_index,stop_interval_index = stop_interval_index, stop_interval_index + offset*int(1/density_step)
                    
                    
                    #the interval should still start higher than it stops    
                    #try:
                    #print(type(iso_interval_pressure_list))
                    #print(start_interval_index,stop_interval_index,"\n")
                    assert start_interval_index < stop_interval_index
                    #assert iso_mean_pressure[profile,start_interval_index] < iso_mean_pressure[profile,stop_interval_index]
                      
                         
                        
                    #average over the bins inside that interval
                    iso_vertical_mean_Shih_flux[profile] = np.nanmean(iso_rolling_mean_Shih_flux[profile,start_interval_index:stop_interval_index])
                    iso_vertical_mean_Osborn_flux[profile] = np.nanmean(iso_rolling_mean_Osborn_flux[profile,start_interval_index:stop_interval_index])
                    iso_vertical_mean_dissipation[profile] = np.log10(np.nanmean(iso_rolling_mean_dissipation[profile,start_interval_index:stop_interval_index]))
                    iso_vertical_mean_pressure[profile] = np.nanmean(iso_mean_pressure[profile,start_interval_index:stop_interval_index])

                    #print( np.nanmean(iso_mean_pressure[profile,start_interval_index:stop_interval_index]))

                    #appendix = [iso_mean_pressure[profile,start_interval_index],iso_mean_pressure[profile,stop_interval_index]]
                    #assert np.all(~np.isnan(appendix))
                    #iso_interval_pressure_list.append(appendix)


                except IndexError:
                    iso_vertical_mean_Shih_flux[profile] = np.nan
                    iso_vertical_mean_Osborn_flux[profile] = np.nan
                    iso_vertical_mean_dissipation[profile] = np.nan
                    iso_vertical_mean_pressure[profile] = np.nan
                    #iso_interval_pressure_list.append([np.nan,np.nan])                 


     
                except AssertionError:
                
                    iso_vertical_mean_Shih_flux[profile] = np.nan
                    iso_vertical_mean_Osborn_flux[profile] = np.nan
                    iso_vertical_mean_dissipation[profile] = np.nan
                    #iso_interval_pressure_list.append([np.nan,np.nan])  
                    iso_vertical_mean_pressure[profile] = np.nan

            else:
                iso_vertical_mean_Shih_flux[profile] = np.nan
                iso_vertical_mean_Osborn_flux[profile] = np.nan
                iso_vertical_mean_dissipation[profile] = np.nan
                #iso_interval_pressure_list.append([np.nan,np.nan])  
                iso_vertical_mean_pressure[profile] = np.nan

        
        #iso_interval_pressure_list = np.asarray(iso_interval_pressure_list)
        #assert np.shape(iso_interval_pressure_list) == np.shape(interval_pressure_list)
 
        #quiver_array.append(iso_vertical_mean_Shih_flux)
        #quiver_pressure.append(iso_vertical_mean_pressure)
    
        quiver_array[offset,:] = iso_vertical_mean_Shih_flux
        quiver_pressure[offset,:] = iso_vertical_mean_pressure
    """
    
    #the colorscheme ['#d95f02','#7570b3','#1b9e77'] stems from colorbrewer (colorbrewer2.org) to be friendly to color blindness and colored printing
    for figure_index,color,label_name,cruise in zip([2,1,0],['#d95f02','#7570b3','#1b9e77'],["summer cruise emb217","winter cruise emb177","autumn cruise emb169"],["emb217","emb177","emb169"]):
        if cruise_name == cruise:
            break
            
   
        
    """            
    def running_mean(x, N):
        length = len(x)
        print("length",length)
        for n in range(length):
            if n > length-N:
                x[n] = np.nanmean(x[n-(N//2):n])
            elif n < N:
                x[n] = np.nanmean(x[n:n+(N//2)])
            else:
                x[n] = np.nanmean(x[n-(N//2):n+(N//2)]) 
        return x
    

    density_list = np.asarray(density_list)
    running_mean_density_list = np.zeros(np.shape(density_list))
    for i in range(len(longitude_list)):
        running_mean_density_list[:,i] = running_mean(density_list[:,i]
    
    def get_isopcycnal(density_array,density_value, eps_pressure):
        isopcnal = []
        for i in range(len(longitude_list)):
            isopycnal.append(eps_pressure[np.nanargmin(np.abs(density_array[:,i])]
    
        return 
    
    """
    
    
    #plot isolines
    #old_x = longitude_list
    #old_y = np.zeros(len(longitude_list))
    ibin = int(np.nanmean(halocline_bin_index_list))
    print(cruise_name,"halocline_density", density_bin_center[ibin])
    halocline_density = density_bin_center[ibin]

    if figure_index == 0:
        axarr_169.contour(isopycnals_lon, eps_pressure, isopycnals_density.T, levels = [halocline_density], colors = "k", linewidths = 3)
    if figure_index == 1:
        axarr_177.contour(isopycnals_lon, eps_pressure, isopycnals_density.T, levels = [halocline_density], colors = "k", linewidths = 3)
    if figure_index == 2:
        axarr_217.contour(isopycnals_lon, eps_pressure, isopycnals_density.T, levels = [halocline_density], colors = "k", linewidths = 3)
        
    #print(cruise_name,"total number of transects =",number_of_transects)        
    #print("total_number_of_valid_profiles",total_number_of_valid_profiles)     


    #if cruise_name == "emb169":
    #    for i in range(len(off_set_list)):
    #        iaxarr.plot(longitude_list,quiver_pressure[i,:], alpha = 0.6)
            
            
    iaxarr.invert_yaxis()
   
    
    #assert np.all(quiver_pressure[0,:] != quiver_pressure[1,:])
    
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
    
    """
    quiver_array = np.asarray(quiver_array)
    print(type(quiver_array[0]),type(longitude_list))
    print(quiver_array[0].dtype)
    print(type(quiver_array[0]))
    plt.plot(quiver_array[0])
    plt.show()
    test1 = np.isnan(quiver_array[0])
    print(test1)
    print(type(test1))
    print(type(test1[0]))
    
    test2 = np.isnan(longitude_list)
    """
        
    longitude_bin_edges = np.linspace(20.5,20.63,30)
    V = np.zeros((len(off_set_list),len(longitude_bin_edges)-1))    
    X = np.zeros((len(off_set_list),len(longitude_bin_edges)-1)) 
    Y = np.zeros((len(off_set_list),len(longitude_bin_edges)-1)) 
    
    levels = []
    levels.append(density_bin_center[ibin-20])
    levels.append(density_bin_center[ibin-10])
    
    for i in range(len(off_set_list)):
        
        levels.append(density_bin_center[ibin+i*10])
        
        quiver_array[i,:] = np.nanmean(iso_rolling_mean_Shih_flux[:,ibin+i*10-5:ibin+i*10+5], axis = 1)
        #quiver_pressure[i,:] = running_mean(np.nanmean(iso_mean_pressure[:,ibin+i*10-5:ibin+i*10+5], axis = 1),20)        
        
        
        flux_row, flux_longitude_list = thesis.remove_nans_from_2_arrays(quiver_array[i],longitude_list)
        V[i,:] = ss.binned_statistic(flux_longitude_list,flux_row, statistic = "mean", bins = longitude_bin_edges)[0]

        if i <= 3:
            lon_offset = (i-2)*0.0005
        else:
            lon_offset = 0

        x_positions = ss.binned_statistic(flux_longitude_list,flux_longitude_list, statistic = "mean", bins = longitude_bin_edges)[0]
        X[i,:] = lon_offset + x_positions


        #for every flux bin calculate the pressure position in the profile closest to the longitude value. 
               
        y_positions = []
        for position in x_positions:
            if np.isnan(position):
                y_positions.append([np.nan])
            else:
                chosen_profile_index = np.nanargmin(np.abs(isopycnals_lon - position))
                chosen_pressure_index = np.nanargmin(np.abs(isopycnals_density[chosen_profile_index,:] - density_bin_center[ibin+i*10]))
                y_position = eps_pressure[chosen_pressure_index]
                sea_floor_depth = total_bathymetry_list[np.nanargmin(np.abs(total_bathymetry_longitude_list - isopycnals_lon[chosen_profile_index]))]
                if y_position > sea_floor_depth:
                    y_position = sea_floor_depth -10
                
                y_positions.append([y_position])
        
        Y[i,:] = np.asarray(y_positions).flatten()
        
        #pressure_row, pressure_longitude_list = thesis.remove_nans_from_2_arrays(quiver_pressure[i],longitude_list)
        #ss.binned_statistic(flux_longitude_list,flux_longitude_list, statistic = "mean", bins = longitude_bin_edges)[0],x,y)
        
        #Y[i,:] = np.nanmean(quiver_pressure[i]) *np.ones(len(longitude_bin_edges)-1)

        
    
    U = np.zeros(np.shape(V))   
                       
    if figure_index == 0:
        axarr_169.contour(isopycnals_lon, eps_pressure, isopycnals_density.T, levels = sorted(levels), colors = "k", linewidth = 0.5, alpha = 0.8)#, linestyles =  "dotted")
        #axarr_169.axvline(20.6083)
        #axarr_169.axvline(20.571)
        axarr_169.axvspan(20.2,20.571, hatch="\\",edgecolor="tab:blue", alpha = 0.2, zorder = 0)
        axarr_169.axvspan(20.571,20.6091, hatch="//",edgecolor="tab:blue", alpha = 0.5, fill=False, zorder = 0)
        axarr_169.axvspan(20.6091,20.7, hatch="\\",edgecolor="tab:blue", alpha = 0.2, zorder = 0)
        
    if figure_index == 1:
        axarr_177.contour(isopycnals_lon, eps_pressure, isopycnals_density.T, levels = sorted(levels), colors = "k", linewidth = 0.5, alpha = 0.8)
        axarr_177.axvspan(20.2,20.60, hatch="\\",edgecolor="tab:blue", alpha = 0.2, zorder = 0)
        axarr_177.axvspan(20.60,20.62, hatch="//",edgecolor="tab:blue", alpha = 0.5, fill=False, zorder = 0)
        axarr_177.axvspan(20.62,20.7, hatch="\\",edgecolor="tab:blue", alpha = 0.2, zorder = 0)
    if figure_index == 2:
        axarr_217.contour(isopycnals_lon, eps_pressure, isopycnals_density.T, levels = sorted(levels), colors = "k", linewidth = 0.5, alpha = 0.8)
        axarr_217.axvspan(20.2,20.571, hatch="\\",edgecolor="tab:blue", alpha = 0.2, zorder = 0)
        axarr_217.axvspan(20.571,20.6091, hatch="//",edgecolor="tab:blue", alpha = 0.5, fill=False, zorder = 0)
        axarr_217.axvspan(20.6091,20.7, hatch="\\",edgecolor="tab:blue", alpha = 0.2, zorder = 0)
            
            
    for edge in longitude_bin_edges:
        all_axarr[figure_index].axvline(edge)
     
    scale = 8   
    pivot = "mid" 
    if figure_index == 0:
        V[3,-6] = np.nan
        Q169 = axarr_169.quiver(X,Y,U,V, units = "xy", zorder = 10, scale = scale, pivot = pivot, color = color)
        axarr_169.set_title("Shih oxygen flux during the "+label_name)
        axarr_169.quiverkey(Q169, 0.95, 0.20, -20, label = "Shih flux \n-20 mmol/(m²d)", coordinates = "axes", labelpos = "W", labelsep= 0.15, angle = 90, color = color)
    elif figure_index == 1:
        V[3,-5] = np.nan
        #V[3:,-5] = -200
        Q177 = axarr_177.quiver(X,Y,U,V, units = "xy", zorder = 10, scale = 6, pivot = pivot, color = color)
        axarr_177.set_title("Shih oxygen flux during the "+label_name)
        axarr_177.quiverkey(Q177, 0.95, 0.25, -20, label = "Shih flux \n-20 mmol/(m²d)", coordinates = "axes", labelpos = "W", labelsep= 0.15, angle = 90, color = color)
    elif figure_index == 2:
        V[2,-4] = np.nanmean(V[2:,-4])
        V[3:,-4] = np.nan
        
        Q217 = axarr_217.quiver(X,Y,U,V, units = "xy", zorder = 10, scale = 4, pivot = pivot, color = color)
        axarr_217.set_title("Shih oxygen flux during the "+label_name)
        axarr_217.quiverkey(Q217, 0.95, 0.25, -10, label = "Shih flux \n-10 mmol/(m²d)", coordinates = "axes", labelpos = "W", labelsep= 0.15, angle = 90, color = color)
        
    all_axarr[figure_index].quiver(X,Y,U,V, units = "xy", zorder = 10, scale = scale, pivot = pivot, headwidth = 2, headlength = 3.5, color = color)



    #plt.show()

    #flux_axarr[1].fill_between(longitude_list,interval_pressure_list[:,0],interval_pressure_list[:,1],color = color, alpha = 0.5, label = label_name)
    #flux_axarr[1].fill_between(longitude_list,iso_interval_pressure_list[:,0],iso_interval_pressure_list[:,1],color = color, alpha = 0.5, label = label_name)
  
    #flux_axarr[1].fill_between(total_bathymetry_longitude_list,total_bathymetry_list, np.ones(len(total_bathymetry_list))*max(total_bathymetry_list),color = "lightgrey", zorder = 1, alpha = 0.8, label = "bathymetry")
    
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
###############################################################################################################
width = 6.2012
height = width / 1.618 #* 0.8 #* 4/3 #

"""
axarr_169.set_xlim(20.5,20.63)
axarr_169.set_ylim(60,125)
axarr_177.set_xlim(20.53,20.665)
axarr_177.set_ylim(43,125)
axarr_217.set_xlim(20.47,20.63)
axarr_217.set_ylim(60,135)
"""

axarr_169.set_xlim(20.554,20.63)
axarr_169.set_ylim(41,115)

axarr_177.set_xlim(20.589,20.6344)
axarr_177.set_ylim(46.78,79.44)

axarr_217.set_xlim(20.55,20.6260)
axarr_217.set_ylim(63,112)

#all_axarr[0].set_xlim(20.47,20.6360)
#all_axarr[0].set_ylim(53.310,94.3507)

axarr_169.set_ylabel("pressure [dbar]")
axarr_169.set_xlabel(r"longitude [$\degree$E]")
axarr_177.set_ylabel("pressure [dbar]")
axarr_177.set_xlabel(r"longitude [$\degree$E]")
axarr_217.set_ylabel("pressure [dbar]")
axarr_217.set_xlabel(r"longitude [$\degree$E]")


axarr_169.invert_yaxis()
axarr_177.invert_yaxis()
axarr_217.invert_yaxis()
all_axarr[0].invert_yaxis()

axarr_169.fill_between(total_bathymetry_longitude_list,total_bathymetry_list, np.ones(len(total_bathymetry_list))*max(total_bathymetry_list), color = "lightgrey", zorder = 5, alpha = 1, label = "bathymetry")
axarr_177.fill_between(total_bathymetry_longitude_list,total_bathymetry_list, np.ones(len(total_bathymetry_list))*max(total_bathymetry_list), color = "lightgrey", zorder = 5, alpha = 1, label = "bathymetry")
axarr_217.fill_between(total_bathymetry_longitude_list,total_bathymetry_list, np.ones(len(total_bathymetry_list))*max(total_bathymetry_list), color = "lightgrey", zorder = 5, alpha = 1, label = "bathymetry")

all_axarr[0].fill_between(total_bathymetry_longitude_list,total_bathymetry_list, np.ones(len(total_bathymetry_list))*max(total_bathymetry_list),color = "lightgrey", zorder = 5, alpha = 1, label = "bathymetry")
all_axarr[1].fill_between(total_bathymetry_longitude_list,total_bathymetry_list, np.ones(len(total_bathymetry_list))*max(total_bathymetry_list),color = "lightgrey", zorder = 5, alpha = 1, label = "bathymetry")
all_axarr[2].fill_between(total_bathymetry_longitude_list,total_bathymetry_list, np.ones(len(total_bathymetry_list))*max(total_bathymetry_list),color = "lightgrey", zorder = 5, alpha = 1, label = "bathymetry")

f_169.set_size_inches(width,height)
f_177.set_size_inches(width,height)
f_217.set_size_inches(width,height)

f_169.tight_layout()
f_169.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/emb169_quiverplot_regions", dpi = 400)
f_177.tight_layout()
f_177.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/emb177_quiverplot_regions", dpi = 400)
f_217.tight_layout()
f_217.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/emb217_quiverplot_regions", dpi = 400)
   
plt.show()
    
    
    
    
    
