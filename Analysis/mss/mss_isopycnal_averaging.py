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

cmap_RdBu = plt.get_cmap('RdBu_r')
cmap_RdBu.set_bad(color = 'lightgrey')
cmap_hot = plt.get_cmap('hot_r')
cmap_hot.set_bad(color = 'lightgrey')


import matplotlib.tri as tri
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
    
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217"],"/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217","/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177"]

rolling_window_size = 10

maximum_reasonable_flux = 500 #float('Inf') #200 #Fluxes above this value will be discarded
acceptable_slope = 2 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 

density_axis = np.linspace(1004,1010.5,20) #maybe change to a non equidistant array?
#density_axis = np.arange(1004,1010.5,0.02)

 
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
     

    pressure_list = []
    dissipation_list = []
    BB_flux_list = []
    Shih_flux_list = []
    Osborn_flux_list = []
    longitude_list = []
    bathymetry_list = []
    
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
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid)
        """
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
        """
        bathymetry,list_of_bathymetrie_indices = results[0]
        #BBL,list_of_BBL_indices = results[1] #not needed here
        BBL_range,list_of_BBL_range_indices = results[2]
        
        eps_N_grid = np.sqrt(eps_N_squared_grid)
        #ozmidov scale
        ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))
        
        #conversion from pressure coordinates to depth
        eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
        bathymetrie_in_m = gsw.z_from_p(bathymetry,np.mean(lat))
        
        eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
        
        distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetrie_in_m,(-1,1))
        boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)

        
        oxygen_flux_Osborn_grid = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_Shih_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0
       

        #compare the pressure value of the lowest valid data point per profile th the neighbouring profiles to determine outliers
        list_of_short_profiles = thesis.get_list_of_short_profiles(number_of_profiles,bathymetry,acceptable_slope)
        
        #hard coded to account for one short profile, that my condtions doesn't recognize
        #if the boolean array contains at least one true
        if np.any(lon == 20.466176666666666):
            list_of_short_profiles.append(np.argmax(lon == 20.466176666666666))
        
        for profile in range(number_of_profiles):
        
        
            number_of_bins = density_axis.size
            density_bins = [ [] for _ in range(number_of_bins) ]
            dissipation_density_bins = [ [] for _ in range(number_of_bins) ]
            Shih_flux_density_bins = [ [] for _ in range(number_of_bins) ]
            Osborn_flux_density_bins = [ [] for _ in range(number_of_bins) ]
            pressure_density_bins = [ [] for _ in range(number_of_bins) ]
        
            iso_Osborn_flux = []
            iso_Shih_flux = []
            iso_pressure = []
            iso_density = []
            iso_dissipation = []
        
            #if the current profile is too short, skip it
            if profile in list_of_short_profiles:
                print(str(lon[profile])+": short profile")
                continue
                
            
            if np.any(eps_oxygen_sat_grid[profile,:] < -0.05):
                temp = eps_oxygen_sat_grid[profile,:]
                print(cruisename,transect_name," is skipped due to ",len(temp[temp<0])," negative oxygen values")
                if len(temp[temp<0]) < 10:
                    print(temp[temp<0])
                continue
                
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
                                    
   
           
            #sort the profile into the density bins
            for value_index,value in enumerate(eps_pot_density_grid[profile]):
            
                bin_number = np.nan
                
                for bin_index,density_bin in enumerate(density_axis):  
                    #print(value,density_bin)      
                    if value <= density_bin:
                        bin_number = bin_index
                        break            
            
                    #Test for edge cases?
                    if value > density_axis[-1]:
                        raise AssertionError
            
                if not np.isnan(bin_number): 
                        Shih_flux_density_bins[bin_number].append(oxygen_flux_Shih_grid[profile,value_index])
                        Osborn_flux_density_bins[bin_number].append(oxygen_flux_Osborn_grid[profile,value_index])
                        pressure_density_bins[bin_number].append(eps_pressure[value_index])
                        density_bins[bin_number].append(eps_pot_density_grid[profile,value_index])
                        dissipation_density_bins[bin_number].append(eps_grid[profile,value_index])
                        
            #average the values inside the density bins
            for bin_index in range(len(density_bins)):
                #print(np.size(density_bins[bin_index]))
                if len(Shih_flux_density_bins[bin_index]) == 0: #test for empty list, ergo no points in that density bin
                    #continue
                    iso_dissipation.append(np.nan)
                    iso_Shih_flux.append(np.nan)
                    iso_Osborn_flux.append(np.nan)
                    iso_pressure.append(np.nan)
                    iso_density.append(np.nan)                
                    
                else:
                    iso_Osborn_flux.append(np.nanmean(Osborn_flux_density_bins[bin_index]))
                    iso_Shih_flux.append(np.nanmean(Shih_flux_density_bins[bin_index]))
                    iso_pressure.append(np.nanmean(pressure_density_bins[bin_index]))
                    iso_density.append(np.mean(density_bins[bin_index]))
                    iso_dissipation.append(np.nanmean(dissipation_density_bins[bin_index]))
                   
                #print(np.shape(density_bins[bin_index]),np.nanmean(pressure_density_bins[bin_index]),density_axis[bin_index])
                
            """
            if profile == np.argmin(np.abs(lon-20.56)):
                test,test_axis = plt.subplots(1)

                test_axis.plot(iso_density,iso_pressure)
                test_axis.plot(eps_pot_density_grid[profile],eps_pressure,"r.")
                test_axis.plot(iso_density,iso_pressure,"ko")
                test_axis.invert_yaxis()
                for line in density_axis:
                    test_axis.axvline(line)
                plt.show()
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
                bathymetry_list.append(bathymetry[profile])  
                pressure_list.append(iso_pressure)
                dissipation_list.append(iso_dissipation)
                #BB_flux_list.append(oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.append(iso_Shih_flux)
                Osborn_flux_list.append(iso_Osborn_flux)
                longitude_list.append(lon[profile])
            
            
            else:
                
                #Sort the current profile into the list   
                bathymetry_list.insert(list_position,bathymetry[profile])    
                dissipation_list.insert(list_position,iso_dissipation)
                pressure_list.insert(list_position,iso_pressure)
                Shih_flux_list.insert(list_position,iso_Shih_flux)
                Osborn_flux_list.insert(list_position,iso_Osborn_flux)
                longitude_list.insert(list_position,lon[profile])


            #print(longitude_list)        
            assert(np.all(longitude_list == sorted(longitude_list)))

            total_number_of_chosen_profiles+=1

        
        #print("removed",outlier_count,"profiles as outliers")
        #print("removed",count_of_short_profiles,"profiles as they did not reach the sea floor")

        
    ###########################################################################################################################
    assert(len(longitude_list) != 0)
    assert(np.all(longitude_list == sorted(longitude_list)))     

    
    ###########################################################################################################################
    
    
    
    #compute rolling isopycnal average
                       
    rolling_mean_Shih_flux = [None] * total_number_of_chosen_profiles                      
    rolling_mean_Osborn_flux = [None] * total_number_of_chosen_profiles
    rolling_arith_mean_dissipation = [None] * total_number_of_chosen_profiles
    mean_pressure = [None] * total_number_of_chosen_profiles

    print(np.shape(pressure_list))
    print(np.shape(dissipation_list))
    print(np.shape(Osborn_flux_list))

    pressure_list = np.asarray(pressure_list)
    dissipation_list = np.asarray(dissipation_list)
    Shih_flux_list = np.asarray(Shih_flux_list)
    Osborn_flux_list = np.asarray(Osborn_flux_list)
                
    #loop over all profiles
    for index in range(0,total_number_of_chosen_profiles): 

        #print(index-(rolling_window_size//2),index+rolling_window_size//2)
        if index >= rolling_window_size//2 and index <= (total_number_of_chosen_profiles+1-rolling_window_size//2):
        #+1 to get the last element through an array slice (Example: a=[0,1], a[:2] = [0,1], a[2] = IndexError)
            
            #isopycnal average over all density bins at teh same time
            mean_pressure[index] = np.nanmean(pressure_list[index-(rolling_window_size//2):index+rolling_window_size//2,:],axis = 0)
            rolling_arith_mean_dissipation[index] = np.nanmean(dissipation_list[index-(rolling_window_size//2):index+rolling_window_size//2,:],axis = 0)
            rolling_mean_Osborn_flux[index] = np.nanmean(Osborn_flux_list[index-(rolling_window_size//2):index+rolling_window_size//2,:],axis = 0)
            rolling_mean_Shih_flux[index] = np.nanmean(Shih_flux_list[index-(rolling_window_size//2):index+rolling_window_size//2,:],axis = 0)
        else:
            mean_pressure[index] = np.nan*np.ones(density_axis.size)
            rolling_arith_mean_dissipation[index] = np.nan*np.ones(density_axis.size)
            rolling_mean_Osborn_flux[index] = np.nan*np.ones(density_axis.size)
            rolling_mean_Shih_flux[index] = np.nan*np.ones(density_axis.size)
     
     
    print("TEST:",np.nanmax(Shih_flux_list),np.nanmin(Shih_flux_list),np.nanmax(rolling_mean_Shih_flux),np.nanmin(rolling_mean_Shih_flux))
            
    #Assure that the isopycnal averaging did not change the overall shape        
    assert(np.shape(mean_pressure) == np.shape(pressure_list))
    rolling_mean_Shih_flux = np.asarray(rolling_mean_Shih_flux)                     
    rolling_mean_Osborn_flux = np.asarray(rolling_mean_Osborn_flux)
    rolling_arith_mean_dissipation = np.asarray(rolling_arith_mean_dissipation)
    mean_pressure = np.asarray(mean_pressure)
                       
    print("total_number_of_chosen_profiles",total_number_of_chosen_profiles)     
    print("number_of_fluxes_over_the_threshold\ttotal_number_of_fluxes\tratio")
    #print("NaN",amount_of_missing_values,total_number_of_fluxes,100*amount_of_missing_values/total_number_of_fluxes,"%")
    #print("0",number_of_zero_flux,total_number_of_fluxes,100*number_of_zero_flux/total_number_of_fluxes,"%")
    #print(">",number_of_fluxes_over_the_threshold,total_number_of_fluxes,100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    #print("Sum:",100*amount_of_missing_values/total_number_of_fluxes + 100*number_of_zero_flux/total_number_of_fluxes + 100*number_of_fluxes_over_the_threshold/total_number_of_fluxes,"%")
    
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
    x_array = longitude_list
    y_array = pressure_list


    x = (np.ones(np.shape(y_array))*np.reshape(x_array,(-1,1))).flatten() #1D to a 2D array back to an 1D array to use the 
    y1 = np.asarray(pressure_list).flatten()
    y2 = np.asarray(mean_pressure).flatten()
    z1 = np.asarray(dissipation_list).flatten()   
    z2 = np.asarray(rolling_arith_mean_dissipation).flatten()

    #print(np.shape(x),np.shape(y1),np.shape(z1),np.shape(z2))
    missing_data_mask1 = ~np.isnan(y1)
    missing_data_mask2 = ~np.isnan(y2)
    
    #remove NaN values as the triangulation is not possible with them
    x1= x[missing_data_mask1]
    x2= x[missing_data_mask2]
    y1= y1[missing_data_mask1]
    y2= y2[missing_data_mask2]
    z1= z1[missing_data_mask1]
    z2= z2[missing_data_mask2]
    #print(np.shape(x),np.shape(y),np.shape(z1),np.shape(z2))
    
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    isbad1 = np.isnan(np.log10(z1))
    triang1 = tri.Triangulation(x1, y1)
    mask1 = np.any(np.where(isbad1[triang1.triangles], True, False), axis=1)
    triang1.set_mask(mask1)
    
    isbad2 = np.isnan(np.log10(z2))
    triang2 = tri.Triangulation(x2, y2)
    mask2 = np.any(np.where(isbad2[triang2.triangles], True, False), axis=1)
    triang2.set_mask(mask2)

    cmap_hot = plt.get_cmap('hot_r')
    cmap_hot.set_bad(color = 'lightgrey')

    f1, ax = plt.subplots(2,2, sharex=True, sharey=True)
    shading= "flat" #'gouraud'
    img00 = ax[0,0].tripcolor(triang1,np.log10(z1), shading=shading,vmin = -9, vmax = -5, cmap = cmap_hot)
    img01 = ax[0,1].tricontourf(triang1,np.log10(z1), levels = 100, vmin = -9, vmax = -5, cmap = cmap_hot, extend = "both") # choose 20 contour levels, just to show how good its interpolation is
    
    print(np.shape(longitude_list),np.shape(bathymetry_list))
    ax[0,0].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax[0,1].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    
    img10 = ax[1,0].tripcolor(triang2,np.log10(z2), shading=shading,vmin = -9, vmax = -5, cmap = cmap_hot)
    img11 = ax[1,1].tricontourf(triang2,np.log10(z2), levels = 100, vmin = -9, vmax = -5, cmap = cmap_hot, extend = "both") # choose 20 contour levels, just to show how good its interpolation is
    
    thesis.colorbar(img00)
    thesis.colorbar(img01, ax[0,1])
    thesis.colorbar(img10)
    thesis.colorbar(img11, ax[1,1])
    
    print(np.shape(longitude_list),np.shape(bathymetry_list))
    ax[1,0].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax[1,1].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    
    ax[0,0].invert_yaxis()
    #ax[1].plot(x,y, 'k,')
    #ax[0].plot(x,y, 'k,')
    ax[0,0].set_title("tripcolor, non averaged")
    ax[0,1].set_title("tripcontourf, non averaged")
    ax[1,0].set_title("tripcolor, isopycnally averaged")
    ax[1,1].set_title("tripcontourf, isopycnally averaged")
    
    f1.suptitle(cruisename+": Comparison based on dissipation measurements")
    f1.set_size_inches(18,10.5)
    f1.tight_layout()
    f1.subplots_adjust(top=0.94)
    f1.savefig("./iso_comparison_dissip_"+cruisename,dpi=300)
    
    ###############################################################################################################
    x = (np.ones(np.shape(y_array))*np.reshape(x_array,(-1,1))).flatten() #1D to a 2D array back to an 1D array to use the 
    y1 = np.asarray(pressure_list).flatten()
    y2 = np.asarray(mean_pressure).flatten()
    z1 = np.asarray(Shih_flux_list).flatten()   
    z2 = np.asarray(rolling_mean_Shih_flux).flatten()

    #print(np.shape(x),np.shape(y1),np.shape(z1),np.shape(z2))
    missing_data_mask1 = ~np.isnan(y1)
    missing_data_mask2 = ~np.isnan(y2)
    
    #remove NaN values as the triangulation is not possible with them
    x1= x[missing_data_mask1]
    x2= x[missing_data_mask2]
    y1= y1[missing_data_mask1]
    y2= y2[missing_data_mask2]
    z1= z1[missing_data_mask1]
    z2= z2[missing_data_mask2]
    #print(np.shape(x),np.shape(y),np.shape(z1),np.shape(z2))
    
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    isbad1 = np.isnan(z1)
    triang1 = tri.Triangulation(x1, y1)
    mask1 = np.any(np.where(isbad1[triang1.triangles], True, False), axis=1)
    triang1.set_mask(mask1)
    
    isbad2 = np.isnan(z2)
    triang2 = tri.Triangulation(x2, y2)
    mask2 = np.any(np.where(isbad2[triang2.triangles], True, False), axis=1)
    triang2.set_mask(mask2)

    f2, ax2 = plt.subplots(2,2, sharex=True, sharey=True)
    shading= "flat" #'gouraud'
    levels = np.linspace(-30,30,20)
        
    img00 = ax2[0,0].tripcolor(triang1,z1, shading=shading,vmin = -30, vmax = +30, cmap = cmap_RdBu)
    img01 = ax2[0,1].tricontourf(triang1,z1, levels = levels, vmin = -30, vmax = +30, cmap = cmap_RdBu, extend = "both") # choose 20 contour levels, just to show how good its interpolation is
    
    print(np.shape(longitude_list),np.shape(bathymetry_list))
    ax2[0,0].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax2[0,1].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    
    img10 = ax2[1,0].tripcolor(triang2,z2, shading=shading,vmin = -30, vmax = +30, cmap = cmap_RdBu)

    img11 = ax2[1,1].tricontourf(triang2,z2, levels = levels, vmin = -30, vmax = +30, cmap = cmap_RdBu, extend = "both") # choose 20 contour levels, just to show how good its interpolation is
    
    thesis.colorbar(img00)
    thesis.colorbar(img01, ax2[0,1])
    thesis.colorbar(img10)
    thesis.colorbar(img11, ax2[1,1])
    
    print(np.shape(longitude_list),np.shape(bathymetry_list))
    ax2[1,0].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax2[1,1].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    
    ax2[0,0].invert_yaxis()
    #ax2[1].plot(x,y, 'k,')
    #ax2[0].plot(x,y, 'k,')
    ax2[0,0].set_title("tripcolor, non averaged")
    ax2[0,1].set_title("tripcontourf, non averaged")
    ax2[1,0].set_title("tripcolor, isopycnally averaged")
    ax2[1,1].set_title("tripcontourf, isopycnally averaged")
        
    f2.suptitle(cruisename+": Comparison based on oxygen fluxes")
    f2.set_size_inches(18,10.5)
    f2.tight_layout()
    f2.subplots_adjust(top=0.94)
    f2.savefig("./iso_comparison_flux_"+cruisename,dpi=300)
    ################################################################################################################
    
    x = (np.ones(np.shape(y_array))*np.reshape(x_array,(-1,1))).flatten() #1D to a 2D array back to an 1D array to use the 
    y1 = np.asarray(pressure_list).flatten()
    y2 = np.asarray(mean_pressure).flatten()
    z1 = np.asarray(dissipation_list).flatten()   
    z2 = np.asarray(rolling_arith_mean_dissipation).flatten()

    #print(np.shape(x),np.shape(y1),np.shape(z1),np.shape(z2))
    missing_data_mask1 = ~np.isnan(y1)
    missing_data_mask2 = ~np.isnan(y2)
    
    #remove NaN values as the triangulation is not possible with them
    x1= x[missing_data_mask1]
    x2= x[missing_data_mask2]
    y1= y1[missing_data_mask1]
    y2= y2[missing_data_mask2]
    z1= z1[missing_data_mask1]
    z2= z2[missing_data_mask2]
    #print(np.shape(x),np.shape(y),np.shape(z1),np.shape(z2))
    
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    isbad1 = np.isnan(np.log10(z1))
    triang1 = tri.Triangulation(x1, y1)
    mask1 = np.any(np.where(isbad1[triang1.triangles], True, False), axis=1)
    triang1.set_mask(mask1)
    
    isbad2 = np.isnan(np.log10(z2))
    triang2 = tri.Triangulation(x2, y2)
    mask2 = np.any(np.where(isbad2[triang2.triangles], True, False), axis=1)
    triang2.set_mask(mask2)
    
    f3, ax3 = plt.subplots(1,2, sharex=True, sharey=True)
    print(type(pressure_list),np.shape(pressure_list))
    ax3[0].tricontourf(triang1,np.log10(z1), levels = 20, vmin = -9, vmax = -5, cmap = cmap_hot)
    ax3[0].plot(longitude_list,pressure_list[:,density_axis.size//2])
    ax3[1].tricontourf(triang2,np.log10(z2), levels = 20, vmin = -9, vmax = -5, cmap = cmap_hot) # choose 20 contour levels, just to show how good its interpolation is
    ax3[1].plot(longitude_list,mean_pressure[:,density_axis.size//2])
    
    #ax3[0].colorbar()
    print(np.shape(longitude_list),np.shape(bathymetry_list))
    ax3[0].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax3[1].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax3[0].invert_yaxis()
    ax3[0].plot(x1,y1, 'k,')
    ax3[1].plot(x2,y2, 'k,')
    f3.suptitle(cruisename+": normal,isopycnal averaged")
    f3.set_size_inches(18,10.5)
    f3.tight_layout()
    f3.subplots_adjust(top=0.94)
    #f3.savefig("./normal_iso_comparison_flux_"+cruisename)        
    #TODO test for monotonous pressure axis
    #TODO plot fluxes
    #TODO compare normal and averaged pressure axis
    
    ################################################################################################################    
      

        
    x_array = longitude_list
    y_array = pressure_list

    x = (np.ones(np.shape(y_array))*np.reshape(x_array,(-1,1))).flatten() #1D to a 2D array back to an 1D array to use the 
    y = np.asarray(mean_pressure).flatten()
    z1 = np.asarray(rolling_arith_mean_dissipation).flatten()   
    z2 = np.asarray(rolling_mean_Shih_flux).flatten()

    #print(np.shape(x),np.shape(y),np.shape(z1),np.shape(z2))
    missing_data_mask = ~np.isnan(y)
    
    #remove NaN values as the triangulation is not possible with them
    x= x[missing_data_mask]
    y= y[missing_data_mask]
    z1= z1[missing_data_mask]
    z2= z2[missing_data_mask]
    #print(np.shape(x),np.shape(y),np.shape(z1),np.shape(z2))
    
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    isbad1 = np.isnan(np.log10(z1))
    triang1 = tri.Triangulation(x, y)
    mask1 = np.any(np.where(isbad1[triang1.triangles], True, False), axis=1)
    triang1.set_mask(mask1)
    
    isbad2 = np.isnan(z2)
    triang2 = tri.Triangulation(x, y)
    mask2 = np.any(np.where(isbad2[triang2.triangles], True, False), axis=1)
    triang2.set_mask(mask2)

    f_flux,flux_axarr = plt.subplots(nrows = 2, ncols = 1, sharey = True, sharex = True) 
    flux_levels = np.linspace(-15,15,20)
    dissip_levels = np.linspace(-9,-5,20)
        
    img_dissip = flux_axarr[0].tricontourf(triang1,np.log10(z1), levels = dissip_levels, cmap = cmap_hot, extend = "both")
    img_flux = flux_axarr[1].tricontourf(triang2,z2, levels = flux_levels, cmap = cmap_RdBu, extend = "both")
    
    print(np.shape(longitude_list),np.shape(bathymetry_list))
    flux_axarr[0].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    flux_axarr[1].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    flux_axarr[0].invert_yaxis()
    
    
    thesis.colorbar(img_dissip,flux_axarr[0]).set_label(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
    thesis.colorbar(img_flux, flux_axarr[1]).set_label(r"oxygen flux [mmol/(m$^2$*d]")
    
    flux_axarr[0].set_ylabel(r"pressure [dbar]")      
    flux_axarr[1].set_xlabel(r"longitude [$\degree$]")    
    flux_axarr[1].set_ylabel(r"pressure [dbar]")

    f_flux.suptitle(cruisename)
    f_flux.set_size_inches(18,10.5)
    f_flux.tight_layout() 
    f_flux.subplots_adjust(top=0.94)
    f_flux.savefig("./isopycnal_average_"+cruisename,dpi=300)

    """
        
    ###############################################################################################################
    

    f_dissip,dissip_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 


                
    dissip_axarr.set_ylabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
    dissip_axarr.set_xlabel(r"longitude [$\degree$E]")       
    """
       
    plt.show()
    

    
    
    
