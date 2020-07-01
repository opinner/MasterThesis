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
plt.rc('savefig', dpi=300)

cmap_RdBu = plt.get_cmap('RdBu_r')
cmap_RdBu.set_bad(color = 'lightgrey')
cmap_hot = plt.get_cmap('hot_r')
cmap_hot.set_bad(color = 'lightgrey')

import matplotlib.ticker as ticker
tick_spacing = 20
import matplotlib.tri as tri
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
    
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217","/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177"]

maximum_reasonable_flux = 500 #float('Inf') #200 #Fluxes above this value will be discarded
acceptable_slope = 2 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 

density_axis = np.linspace(1004,1010.5,50) #maybe change to a non equidistant array?

averaging_intervals_borders = [20.55,20.62]
#averaging_intervals_borders = np.linspace(20.48,20.7,4)
 
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
    temperature_list = []
    salinity_list = []
    oxygen_sat_list = []
    density_list = []
    
    
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
            eps_salinity_grid = data["eps_salinity_grid"] 
            
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
            
            #collects all values in a certain density interval
            density_bins = [ [] for _ in range(number_of_bins) ]
            dissipation_density_bins = [ [] for _ in range(number_of_bins) ]
            Shih_flux_density_bins = [ [] for _ in range(number_of_bins) ]
            Osborn_flux_density_bins = [ [] for _ in range(number_of_bins) ]
            pressure_density_bins = [ [] for _ in range(number_of_bins) ]
            salinity_density_bins = [ [] for _ in range(number_of_bins) ]
            temperature_density_bins = [ [] for _ in range(number_of_bins) ]
            oxygen_sat_density_bins = [ [] for _ in range(number_of_bins) ]
                                            
            #collects only the mean in a certain density interval
            iso_Osborn_flux = []
            iso_Shih_flux = []
            iso_pressure = []
            iso_density = []
            iso_dissipation = []
            iso_temperature = []
            iso_salinity = []
            iso_oxygen_sat = []
            
            #if the current profile is too short, skip it
            if profile in list_of_short_profiles:
                print(str(lon[profile])+": short profile")
                continue
                
            
            #if np.nanmean(eps_oxygen_sat_grid[profile]) < 0:
            #    print(cruisename,transect_name," is skipped due to negative oxygen values")
            #    continue  
               
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
                        salinity_density_bins[bin_number].append(eps_salinity_grid[profile,value_index])
                        temperature_density_bins[bin_number].append(eps_consv_temperature_grid[profile,value_index])
                        oxygen_sat_density_bins[bin_number].append(eps_oxygen_sat_grid[profile,value_index])
            
                        
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
                    iso_temperature.append(np.nan)     
                    iso_salinity.append(np.nan)     
                    iso_oxygen_sat.append(np.nan)     
            
                else:
                    iso_Osborn_flux.append(np.nanmean(Osborn_flux_density_bins[bin_index]))
                    iso_Shih_flux.append(np.nanmean(Shih_flux_density_bins[bin_index]))
                    iso_pressure.append(np.nanmean(pressure_density_bins[bin_index]))
                    iso_density.append(np.mean(density_bins[bin_index]))
                    iso_dissipation.append(np.nanmean(dissipation_density_bins[bin_index]))
                    iso_temperature.append(np.nanmean(temperature_density_bins[bin_index])) 
                    iso_salinity.append(np.nanmean(salinity_density_bins[bin_index]))    
                    iso_oxygen_sat.append(np.nanmean(oxygen_sat_density_bins[bin_index]))    
                                       
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
                temperature_list.append(iso_temperature)
                salinity_list.append(iso_salinity)
                oxygen_sat_list.append(iso_oxygen_sat)
                density_list.append(iso_density)
            
            else:
                
                #Sort the current profile into the list   
                bathymetry_list.insert(list_position,bathymetry[profile])    
                dissipation_list.insert(list_position,iso_dissipation)
                pressure_list.insert(list_position,iso_pressure)
                Shih_flux_list.insert(list_position,iso_Shih_flux)
                Osborn_flux_list.insert(list_position,iso_Osborn_flux)
                longitude_list.insert(list_position,lon[profile])
                temperature_list.insert(list_position,iso_temperature)
                salinity_list.insert(list_position,iso_salinity)
                oxygen_sat_list.insert(list_position,iso_oxygen_sat)
                density_list.insert(list_position,iso_density)
                
            #print(longitude_list)        
            assert(np.all(longitude_list == sorted(longitude_list)))

            total_number_of_chosen_profiles+=1

        
        #print("removed",outlier_count,"profiles as outliers")
        #print("removed",count_of_short_profiles,"profiles as they did not reach the sea floor")

        
    ###########################################################################################################################
    assert(len(longitude_list) != 0)
    assert(np.all(longitude_list == sorted(longitude_list)))     

    if cruisename == "emb217":
        np.savetxt("./bathymetry_emb217.txt",np.transpose([longitude_list,bathymetry_list]),header = "longitude [degree East] \t depth in pressure coordinates [dbar]", fmt = ["%2f","%2.2f"])
    ###########################################################################################################################
    
    
    
    #compute isopycnal average in the longitude intervals                    
    binned_mean_Shih_flux = [None] * (len(averaging_intervals_borders)+1)                  
    binned_mean_Osborn_flux = [None] * (len(averaging_intervals_borders)+1) 
    binned_arith_mean_dissipation = [None] * (len(averaging_intervals_borders)+1) 
    binned_mean_salinity =  [None] * (len(averaging_intervals_borders)+1) 
    binned_mean_temperature =  [None] * (len(averaging_intervals_borders)+1) 
    binned_mean_oxygen_sat =   [None] * (len(averaging_intervals_borders)+1) 
    binned_mean_pressure = [None] * (len(averaging_intervals_borders)+1) 
    binned_mean_density = [None] * (len(averaging_intervals_borders)+1) 
    
    print(np.shape(pressure_list))
    print(np.shape(dissipation_list))
    print(np.shape(Osborn_flux_list))

    longitude_list = np.asarray(longitude_list) #1D
    pressure_list = np.asarray(pressure_list) #2D
    dissipation_list = np.asarray(dissipation_list) #2D
    Shih_flux_list = np.asarray(Shih_flux_list) #2D
    Osborn_flux_list = np.asarray(Osborn_flux_list) #2D
    temperature_list = np.asarray(temperature_list)
    salinity_list = np.asarray(salinity_list)
    oxygen_sat_list = np.asarray(oxygen_sat_list)           
    density_list = np.asarray(density_list)

    averaging_interval_indices = [] 
    for interval_border in averaging_intervals_borders:
        averaging_interval_indices.append(np.argmin(np.abs(longitude_list-interval_border)))
        
    #append the length of the longitude list to close the last interval    
    averaging_interval_indices.append(len(longitude_list))
    
    #print(longitude_list)
    #print(averaging_interval_indices)
    
    #average inside the longitude intervals
    start = 0    
    for interval_index, border_index in enumerate(averaging_interval_indices):
        stop = border_index
      
        #print(start,stop)
        #isopycnal average over all density bins at the same time
        binned_mean_pressure[interval_index] = np.nanmean(pressure_list[start:stop,:],axis = 0)
        binned_arith_mean_dissipation[interval_index] = np.nanmean(dissipation_list[start:stop,:],axis = 0)
        binned_mean_Osborn_flux[interval_index] = np.nanmean(Osborn_flux_list[start:stop,:],axis = 0)
        binned_mean_Shih_flux[interval_index] = np.nanmean(Shih_flux_list[start:stop,:],axis = 0)
        binned_mean_salinity[interval_index] = np.nanmean(salinity_list[start:stop,:],axis = 0)
        binned_mean_temperature[interval_index] = np.nanmean(temperature_list[start:stop,:],axis = 0)
        binned_mean_oxygen_sat[interval_index] = np.nanmean(oxygen_sat_list[start:stop,:],axis = 0)
        binned_mean_density[interval_index] = np.nanmean(density_list[start:stop,:],axis = 0)
                
        start = stop

     
    #print("TEST:",np.nanmax(Shih_flux_list),np.nanmin(Shih_flux_list),np.nanmax(binned_mean_Shih_flux),np.nanmin(binned_mean_Shih_flux))
      
    binned_mean_Shih_flux = np.asarray(binned_mean_Shih_flux)                     
    binned_mean_Osborn_flux = np.asarray(binned_mean_Osborn_flux)
    binned_arith_mean_dissipation = np.asarray(binned_arith_mean_dissipation)
    binned_mean_pressure = np.asarray(binned_mean_pressure)
    binned_mean_salinity = np.asarray(binned_mean_salinity)
    binned_mean_temperature = np.asarray(binned_mean_temperature)
    binned_mean_oxygen_sat = np.asarray(binned_mean_oxygen_sat)
    binned_mean_density = np.asarray(binned_mean_density)
        
    print(np.shape(binned_mean_Shih_flux))
                       
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

    f,axarr = plt.subplots(nrows = 3, ncols = len(averaging_interval_indices), sharey = True, sharex = "row") 

    for index in range(len(averaging_interval_indices)):
        
        
        nanmask = ~np.isnan(binned_mean_pressure[index,:])
        cleaned_mean_pressure = binned_mean_pressure[index,:][nanmask]
        cleaned_salinity = binned_mean_salinity[index,:][nanmask]
        cleaned_temperature = binned_mean_temperature[index,:][nanmask]
        cleaned_density = binned_mean_density[index,:][nanmask]
        cleaned_dissipation = binned_arith_mean_dissipation[index,:][nanmask]
        cleaned_oxygen = binned_mean_oxygen_sat[index,:][nanmask]
        cleaned_Osborn_flux = binned_mean_Osborn_flux[index,:][nanmask]
        cleaned_Shih_flux = binned_mean_Shih_flux[index,:][nanmask]
        
        salinity_plot_axis = axarr[0,index].twiny()
        density_plot_axis =  axarr[0,index].twiny()
        density_plot_axis.axis("off")
        
        dissipation_plot_axis = axarr[1,index].twiny()
        axarr[1,index].patch.set_visible(False)
        dissipation_plot_axis.set_zorder(0)
        dissipation_plot_axis.set_zorder(-10)

    
        print(cleaned_mean_pressure)

        
        
        axarr[0,index].plot(cleaned_temperature,cleaned_mean_pressure, c = "tab:red")
        salinity_plot_axis.plot(cleaned_salinity,cleaned_mean_pressure, c = "tab:green")
        density_plot_axis.plot(cleaned_density,cleaned_mean_pressure,c = "k")
        density_plot_axis.plot(cleaned_density,cleaned_mean_pressure,".",c = "k")
                    
        axarr[1,index].plot(cleaned_oxygen,cleaned_mean_pressure,c = "tab:blue")
                
        axarr[2,index].plot(cleaned_Osborn_flux,cleaned_mean_pressure,".", c = "tab:green", alpha = 0.4)
        axarr[2,index].plot(cleaned_Osborn_flux,cleaned_mean_pressure,ls = "-", c = "tab:green", alpha = 0.4)            
        axarr[2,index].plot(cleaned_Shih_flux,cleaned_mean_pressure,".", c = "tab:blue")
        axarr[2,index].plot(cleaned_Shih_flux,cleaned_mean_pressure,c = "tab:blue")

        dissipation_plot_axis.plot(np.log10(cleaned_dissipation),cleaned_mean_pressure,c = "k", alpha = 0.4)

        axarr[2,index].set_xlim(-20,20)
        
        
    axarr[0,0].invert_yaxis()
    axarr[0,0].yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))      
    #axarr.set_ylabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
    #axarr.set_xlabel(r"longitude [$\degree$E]")     
    f.suptitle(cruisename)  
    f.set_size_inches(18,10.5)
    f.tight_layout() 
    f.subplots_adjust(top=0.915)

    f.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+"isopycnal_averaged_overview_"+cruisename)
    
    
    b,bxarr = plt.subplots(nrows = 1, ncols = 1) 
    bxarr.invert_yaxis()
    bxarr.plot(longitude_list,bathymetry_list)
    for border in averaging_intervals_borders:
        bxarr.axvline(border, c = "k")
    bxarr.xaxis.set_major_locator(ticker.MultipleLocator(0.01))   
    b.set_size_inches(18,10.5)
    b.tight_layout() 
    b.subplots_adjust(top=0.915)
    b.suptitle(cruisename)
    b.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+"isopycnal_averaged_overview_"+cruisename+"_intervals")   
    
    #plt.show()
    

