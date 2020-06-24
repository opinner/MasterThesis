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

rolling_window_size = 12

maximum_reasonable_flux = 500 #float('Inf') #200 #Fluxes above this value will be discarded
acceptable_slope = 2 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 

density_axis = np.linspace(1004,1010.5,10) #maybe change to a non equidistant array?


 
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
    
    """
    if cruisename == "emb217":
        upper_bound_halocline_as_density = 1006.4 #1005.75
        lower_bound_halocline_as_density = 1008.5 #1006.25
    elif cruisename == "emb177":
        upper_bound_halocline_as_density = 1006.9 #1007.4
        lower_bound_halocline_as_density = 1008.2 #1007.9   
    elif cruisename == "emb169":
        pass
        #upper_bound_halocline_as_density = 1006.9 #1007.4
        #lower_bound_halocline_as_density = 1008.2 #1007.9  
    """       

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
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid,height_above_ground = height_above_ground)
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
                
            
            if np.nanmean(eps_oxygen_sat_grid[profile]) < 0:
                print(cruisename,transect_name," is skipped due to negative oxygen values")
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
                if len(Shih_flux_density_bins[bin_index]) == 0: #test for empty list
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
            if lon
            test,test_axis = plt.subplots(1)

            test_axis.plot(iso_density,iso_pressure)
            test_axis.plot(iso_density,iso_pressure,"o")
            test_axis.plot(eps_pot_density_grid[profile],eps_pressure,"r.")
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
          

    x_array = longitude_list
    y_array = pressure_list


    x = (np.ones(np.shape(y_array))*np.reshape(x_array,(-1,1))).flatten()
    y = np.asarray(pressure_list).flatten()
    z = np.asarray(dissipation_list).flatten()

    print(np.shape(x),np.shape(y),np.shape(z))
    missing_data_mask = ~np.isnan(y)
    
    #remove NaN values as the triangulation is not possible with them
    x= x[missing_data_mask]
    y= y[missing_data_mask]
    z= z[missing_data_mask]
    print(np.shape(x),np.shape(y),np.shape(z))
    
    # Create the Triangulation; no triangles so Delaunay triangulation created.
    isbad = np.isnan(np.log10(z))
    triang = tri.Triangulation(x, y)
    mask = np.any(np.where(isbad[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)

        
    cmap_hot = plt.get_cmap('hot_r')
    cmap_hot.set_bad(color = 'lightgrey')

    f, ax = plt.subplots(1,2, sharex=True, sharey=True)
    ax[0].tripcolor(triang,np.log10(z), shading='gouraud',vmin = -9.5, vmax = -5, cmap = cmap_hot)
    ax[1].tricontourf(triang,np.log10(z), levels = 20, vmin = -9.5, vmax = -5, cmap = cmap_hot) # choose 20 contour levels, just to show how good its interpolation is
    
    print(np.shape(longitude_list),np.shape(bathymetry_list))
    ax[0].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax[1].fill_between(longitude_list,bathymetry_list,np.ones(np.shape(bathymetry_list))*np.max(bathymetry_list), color = "lightgrey")
    ax[0].invert_yaxis()
    #ax[1].plot(x,y, 'ko ')
    #ax[0].plot(x,y, 'ko ')
       
    plt.show()

    
    ###########################################################################################################################
    
    
    
    #compute rolling isopycnal average
                       
                         
    rolling_mean_flux = [None] * total_number_of_chosen_profiles
    rolling_arith_mean_dissipation = [None] * total_number_of_chosen_profiles

    #loop over all profiles
    for index in range(total_number_of_chosen_profiles):
   
   
        #loop over all density bins 
        if rolling_window_size ==1:
            rolling_mean_flux[index] = mean_flux[index]
            rolling_median_flux[index] = mean_flux[index]


            
        else:
            try:
                
                rolling_mean_flux[index] = np.nanmean(mean_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_arith_mean_dissipation[index] = np.nanmean(arith_mean_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2])

    
                #print(index,longitude_list[index],np.round(mean_flux[index-rolling_window_size//2:index+rolling_window_size//2],3))
    

            except (IndexError,ValueError):
                rolling_mean_flux[index] = np.nan

                rolling_arith_mean_dissipation[index] = np.nan

    
    
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
    
    flux_axarr.legend(handles=[bathymetrie_label]) #loc=8

    flux_axarr.set_ylim((-20,1))
    if cruisename == "emb177":
        flux_axarr.set_ylim((-50,1))    
            
    flux_axarr.set_xlabel(r"longitude [$\degree$]")    
    flux_axarr.set_ylabel(r"oxygen flux [mmol/(m$^2$*d]")
    flux_axarr.legend()
    
    f_flux.set_size_inches(18,10.5)
    f_flux.tight_layout() 
    f_flux.subplots_adjust(top=0.94)


    
        
    ###############################################################################################################
    

    f_dissip,dissip_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 


                
    dissip_axarr.set_ylabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
    dissip_axarr.set_xlabel(r"longitude [$\degree$E]")       

       
    plt.show()
    

    
    
    
