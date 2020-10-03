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
import warnings
warnings.filterwarnings('ignore')

#red sequential color scheme
red_iterator = iter(['#fed98e','#fe9929','#cc4c02'])
blue_iterator = iter(['#cbc9e2','#9e9ac8','#6a51a3'])
green_iterator = iter(['#b2e2e2','#66c2a4','#238b45'])

rolling_window_size = 12


#plot halocline (and other isolines) in pressure and density coordinates. 
f_iso,iso_axis = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = "col")

#plot mean oxygen flux dependent on in pressure and density coordinates (once for shih and once for Osborn)
f_box, box_axis = plt.subplots(nrows = 3, sharex = True) #, sharey = True)

maximum_reasonable_flux = float('Inf') #200 #Fluxes with absolute values above this cut off value will be discarded



number_of_dissipation_subplots = 1 #Decide if both the mean and the median subplots is shown or only the mean

datafile_paths = ["/home/ole/windows/processed_mss/emb177/TS1_10.npz","/home/ole/windows/processed_mss/emb217/TR1-10.npz","/home/ole/windows/processed_mss/emb169/TS112.npz"]
#datafile_paths = ["/home/ole/windows/processed_mss/emb169/TS112.npz"]
#datafile_paths = ["/home/ole/windows/processed_mss/emb177/TS1_10.npz"]

emb169_datafile_paths = ["/home/ole/windows/processed_mss/emb169/TS11.npz","/home/ole/windows/processed_mss/emb169/TS12.npz","/home/ole/windows/processed_mss/emb169/TS18.npz"]
emb177_datafile_paths = ["/home/ole/windows/processed_mss/emb177/TS1_5.npz","/home/ole/windows/processed_mss/emb177/TS1_10.npz","/home/ole/windows/processed_mss/emb177/TS1_12.npz"]
emb217_datafile_paths = ["/home/ole/windows/processed_mss/emb217/TR1-7.npz","/home/ole/windows/processed_mss/emb217/TR1-9.npz","/home/ole/windows/processed_mss/emb217/TR1-10.npz"]



list_of_bad_profiles,reasons = np.loadtxt("./data/list_of_bad_profiles.txt", dtype=str, unpack=True)

interior_longitude_interval = [20.53,20.57]
edge_longitude_interval = [20.57,20.62]
      
box_sizes = np.arange(0.1,3,0.05) #in density units [kg/m³]
#box_sizes = np.arange(0.1,3,0.3) 
box_sizes_in_dbar = np.zeros(box_sizes.size)
 
label_list = []

datafile_paths =  emb169_datafile_paths +  emb177_datafile_paths +  emb217_datafile_paths
for file_index,datafile_path in enumerate(datafile_paths):

    transect_halocline_depths = []
    transect_halocline_densities = []

    splitted_filename = datafile_path.split("/")
    cruise_name = splitted_filename[5][0:6]
    transect_name = splitted_filename[-1][:-4]

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
    results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid)
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
    
    #get the turbulent diffusivity im Osborn Model with three different parametrizations
    turbulent_diffusivity_Osborn_grid = thesis.get_turbulent_diffusivity_Osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
    turbulent_diffusivity_BB_grid = thesis.get_turbulent_diffusivity_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
    turbulent_diffusivity_Shih_grid = thesis.get_turbulent_diffusivity_Shih(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)           
   
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
    oxygen_flux_BB_grid =  - turbulent_diffusivity_BB_grid * oxygen_gradient_grid * unit_conversion_grid
    oxygen_flux_Shih_grid = - turbulent_diffusivity_Shih_grid * oxygen_gradient_grid * unit_conversion_grid
   
   
    interior_Shih_flux_for_different_intervals = []
    interior_Osborn_flux_for_different_intervals = []

    edge_Shih_flux_for_different_intervals = []
    edge_Osborn_flux_for_different_intervals = []
    
    #profile for the example conversion from density to pressure space 
    example_profile_index = np.argmin(abs(lon-20.572))
    
    
    for width_index,density_box_width in enumerate(box_sizes):
    
    
        dissipation_list = []
        BB_flux_list = []
        Shih_flux_list = []
        Osborn_flux_list = []
        longitude_list = []
        bathymetry_list = []
    
        total_number_of_valid_profiles = 0
    
        for profile in range(number_of_profiles):
    
            
            if "_".join((cruise_name,transect_name,str(profile))) in list_of_bad_profiles:
                if profile == example_profile_index: example_profile_index+=1
                #print("_".join((cruise_name,transect_name,str(profile))),"skipped")
                continue
            
            #print("\n")
             
            total_number_of_valid_profiles+=1
                
            halocline_depth,halocline_density,halocline_index = thesis.get_halocline_and_halocline_density(eps_pressure,eps_oxygen_sat_grid[profile],eps_salinity_grid[profile],eps_consv_temperature_grid[profile],eps_pot_density_grid[profile])
        
            
            if np.isnan(halocline_depth):
                #find the correct position in the sorted list
                for index,value in enumerate(longitude_list):
                    if value > lon[profile]:
                        list_position = index
                        break
                    elif index == len(longitude_list)-1:
                        list_position = len(longitude_list)
                        break
                    
            
                #print("no halocline",lon[profile])
                
                #Append Nans
                if len(longitude_list) == 0:    
                    dissipation_list.append(np.nan)
                    #BB_flux_list.append(oxygen_flux_BB_grid[profile,from_index:to_index])
                    Shih_flux_list.append(np.nan)
                    Osborn_flux_list.append(np.nan)
                    longitude_list.append(lon[profile])
                    bathymetry_list.append(bathymetry[profile])
                    if width_index == 0:    
                        transect_halocline_depths.append(halocline_depth)
                        transect_halocline_densities.append(halocline_density)
                                                
                else:
                    
                    #Sort the current profile into the list            
                    dissipation_list.insert(list_position,np.nan)
                    #BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                    Shih_flux_list.insert(list_position,np.nan)
                    Osborn_flux_list.insert(list_position,np.nan)
                    longitude_list.insert(list_position,lon[profile])
                    bathymetry_list.insert(list_position,bathymetry[profile])
                    if width_index == 0: 
                        transect_halocline_depths.insert(list_position,halocline_depth)
                        transect_halocline_densities.insert(list_position,halocline_density)
                continue
        
            #choose the vertical averaging interval dependent on the box size
            #around the halocline
            from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile]) - (halocline_density - density_box_width/2)))     
            to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])- (halocline_density + density_box_width/2)))
            
            #from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile]) - (halocline_density - 0)))     
            #to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])- (halocline_density + density_box_width)))
                        
            maximum_halocline_thickness = 20
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

                    
            if profile == example_profile_index:
                box_sizes_in_dbar[width_index] = np.abs(eps_pressure[to_index] - eps_pressure[from_index])
            
                        
            if cruise_name == "emb177" and file_index%3 == 0:
                pass
                #print(from_index,to_index)
           

                          
            """
            if profile == example_profile_index and file_index%3 == 0:
                box_sizes_in_dbar[width_index] = np.abs(eps_pressure[to_index] - eps_pressure[from_index])
                test,test_ax = plt.subplots(1)
                test_ax.plot(eps_pot_density_grid[profile],-eps_pressure)
                test_ax.axhline(-eps_pressure[from_index], ls = "--", c = "k")
                test_ax.axhline(-eps_pressure[halocline_index], c = "r")
                test_ax.axhline(-eps_pressure[to_index], c = "k")
                test_ax.set_title(str(np.abs(eps_pressure[to_index] - eps_pressure[from_index])))
            """
           
            assert from_index <= to_index

                
                                            
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
                if width_index == 0: 
                    transect_halocline_depths.append(halocline_depth)
                    transect_halocline_densities.append(halocline_density)                           
            else:
                
                #Sort the current profile into the list            
                dissipation_list.insert(list_position,eps_grid[profile,from_index:to_index])
                #BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.insert(list_position,oxygen_flux_Shih_grid[profile,from_index:to_index])
                Osborn_flux_list.insert(list_position,oxygen_flux_Osborn_grid[profile,from_index:to_index])
                longitude_list.insert(list_position,lon[profile])
                bathymetry_list.insert(list_position,bathymetry[profile])
                if width_index == 0: 
                    transect_halocline_depths.insert(list_position,halocline_depth)
                    transect_halocline_densities.insert(list_position,halocline_density)
                #print(np.nanmean(oxygen_flux_Shih_grid[profile,from_index:to_index]))

            

        
        ###########################################################################################################################
        #print(len(longitude_list),"used profiles")
        assert(len(longitude_list) != 0)
        assert(np.all(longitude_list == sorted(longitude_list)))      

        #compute mean and std over the saved intervals
        mean_Osborn_flux = [None] * total_number_of_valid_profiles
        mean_Shih_flux = [None] * total_number_of_valid_profiles
        arith_mean_dissipation = [None] * total_number_of_valid_profiles


        #average vertically in the density interval
        for index in range(len(longitude_list)):
            temp_Shih_flux = np.asarray(Shih_flux_list[index])    
            #temp_Shih_flux[np.abs(temp_Shih_flux)>maximum_reasonable_flux] = np.nan
            mean_Shih_flux[index] = np.nanmean(temp_Shih_flux)
                        
            temp_Osborn_flux = np.asarray(Osborn_flux_list[index])
            #temp_Osborn_flux[np.abs(temp_Osborn_flux)>maximum_reasonable_flux] = np.nan  
            mean_Osborn_flux[index] = np.nanmean(temp_Osborn_flux)
     
            arith_mean_dissipation[index] = np.log10(np.nanmean(dissipation_list[index]))

                        
                           
        longitude_list = np.asarray(longitude_list)
        
        #get indices for the longitude interval 
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
    
    print(cruise_name)
    #'#d95f02','#7570b3','#1b9e77'
    
    for figure_index,color_iterator,labelname,cruise in zip([2,1,0],[red_iterator,blue_iterator,green_iterator],["summer cruise emb217","winter cruise emb177","autumn cruise emb169"],["emb217","emb177","emb169"]):
        if cruise_name == cruise:
            break

    
    color = next(color_iterator)
            
    #box_axis[figure_index].plot(box_sizes,np.abs(edge_Osborn_flux_for_different_intervals), c = color, label = " ".join((cruise_name,transect_name)))
    #box_axis[figure_index].plot(box_sizes,np.abs(interior_Osborn_flux_for_different_intervals), "--", c = color)

    box_axis[figure_index].plot(box_sizes,np.abs(edge_Shih_flux_for_different_intervals), c = color, label = " ".join((cruise_name,transect_name)))
    box_axis[figure_index].plot(box_sizes,np.abs(interior_Shih_flux_for_different_intervals), "--", c = color)

    iso_axis[figure_index,0].plot(longitude_list,transect_halocline_depths, c = color, label = " ".join((cruise_name,transect_name)))
    iso_axis[figure_index,1].plot(longitude_list,transect_halocline_densities, c = color, label = " ".join((cruise_name,transect_name)))
    iso_axis[figure_index,0].plot(longitude_list,bathymetry_list,"k")



    label_list.append(mpatches.Patch(color= color, label= " ".join((transect_name))))


    #only for the first time per cruise
    if file_index%3 == 0:
        box_m0 = box_axis[figure_index].twiny()
        #box_m1 = box_axis[1].twiny()
        locs = box_axis[figure_index].get_xticks()
        #box_m0.set_xticks(locs)
        #box_m0.set_xticklabels(box_sizes_in_dbar, alpha = 0.7)
        
        #make ticks smaller, little transparent and closer to the axis
        plt.setp(box_m0.get_xticklabels(), alpha=0.8)
        box_m0.tick_params(axis='x', which='major', pad= -2 )
        box_m0.tick_params(axis='both', which='major', labelsize= MINI_SIZE)
        box_m0.tick_params(axis='both', which='minor', labelsize= MINI_SIZE)
        
        #dummy plot
        box_m0.plot(box_sizes_in_dbar,np.ones(len(box_sizes_in_dbar)), alpha = 0)
        
        
        print(box_sizes_in_dbar)
        
        box_sizes_in_dbar = -1 * np.zeros(box_sizes.size)
        
        #set label only at the top
        if figure_index == 0:
            box_m0.set_xlabel("Averaging interval in meter")
    
    
    if file_index == 0:
        #iso_axis[figure_index,0].legend(handles = label_list, loc = "lower center", ncol=3)
        #bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
        #iso_axis[0,1].legend(handles = label_list, loc = "lower center", ncol=3)
        #bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
        
        label_list = []
        label_list.append(mlines.Line2D([0], [0], color= "k", ls = "-", label='basin edge'))
        label_list.append(mlines.Line2D([], [], color = "k", ls = "--", label='basin interior'))
        box_axis[figure_index].legend(handles = label_list, loc = "upper right") #, bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
        box_axis[figure_index].legend(handles = label_list, loc = "upper right") #, bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
        #bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
        label_list = []
    



iso_axis[0,0].invert_yaxis()
iso_axis[0,1].invert_yaxis()





#box_axis[1].set_ylabel(r"|$\langle$ Osborn OF $\rangle$|")
box_axis[1].set_ylabel(r"|$\langle$ Shih OF $\rangle$|")
box_axis[2].set_xlabel("Averaging interval in kg/m³")

iso_axis[0,0].set_title("Halocline depths")
iso_axis[0,1].set_title("Halocline densities")
iso_axis[1,0].set_ylabel("pressure [dbar]")
iso_axis[1,1].set_ylabel("potential density [kg/m³]")
iso_axis[2,0].set_xlabel(r"longitude [$\degree$E]")
iso_axis[2,1].set_xlabel(r"longitude [$\degree$E]")

#box_axis[0].legend(handles = label_list, loc = "upper right")
#box_axis[1].legend(handles = label_list, loc = "upper right")

f_iso.set_size_inches(6.2012,6.2012/1.618)
f_iso.subplots_adjust(top=0.924,bottom=0.154,left=0.111,right=0.984,hspace=0.467,wspace=0.36)
f_iso.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/halocline_depths_and_densities.png", dpi = 600)

f_box.suptitle("Impact of the bounded vertical averaging interval")
f_box.set_size_inches(6.2012,6.2012/1.618)
f_box.subplots_adjust(top=0.827,bottom=0.136,left=0.114,right=0.946,hspace=0.7,wspace=0.2)
f_box.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/vertical_averaging_comparison_bounded.png", dpi = 600)

plt.show()
    
    
    
    


