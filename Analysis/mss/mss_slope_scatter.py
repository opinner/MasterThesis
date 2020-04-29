############################################################
#this program loads all profiles from a cruise, removes outliers
#and makes a scatter plot of dissipation in the lowermost 10m 
#and the slope of the basin. 

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

import geopy.distance as geo
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
#warnings.filterwarnings('ignore')
    
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169","/home/ole/share-windows/processed_mss/emb177"]

height_above_ground = 10
maximum_reasonable_flux = 150
acceptable_slope = 5 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 
flux_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_flux_percentile = 97.72
dissip_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_dissip_percentile = 97.72

  
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print(cruisename)
             

    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".npz":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    #print(DATAFILENAMES)
    
    number_of_transects = len(DATAFILENAMES)+1
    
    #2 borders = 3 intervals
    oxygen_flux_BB_statistic = [None] * number_of_transects
    oxygen_flux_Shih_statistic = [None] * number_of_transects
    oxygen_flux_Osborn_statistic = [None] * number_of_transects
    oxygen_flux_statistic = [None] * number_of_transects

    slope_statistic = [None] * number_of_transects
    bathymetrie_statistic = [None] * number_of_transects
    dissipation_statistic = [None] * number_of_transects
    lon_statistic = [None] * number_of_transects
    list_of_transect_names = [None] * number_of_transects
    
    #all_dissipation_statistic = np.asarray([])
    
    profile_count = 0
    
    for transect_index,DATAFILENAME in enumerate(DATAFILENAMES):
    
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        transect_name = DATAFILENAME[:-4]
    
        #skip the short "S206" transects
        if transect_name[0:4] == "S106":
            print(transect_name,"skipped")
            continue
            
        print("\n",transect_name)
            
        list_of_transect_names[transect_index] = transect_name
        
        data = np.load(datafile_path)
        
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
        eps_oxygen_grid = data["eps_oxygen_grid"] 
        
        eps_N_squared_grid = data["eps_N_squared_grid"]
        eps_density_grid = data["eps_density_grid"]
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
            
            
        #boundary_check_grid = np.zeros(np.shape(ozmidov_scale_grid))
        #check if ozimidov scale is bigger than the distance from the ground
        #for i in range(number_of_profiles):
        #density_grid[ozmidov_scale_grid,:])))
        
        """
        Gamma_Osborn_eps_grid = thesis.Osborn(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity 
        turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid<0] = np.nan
        oxygen_flux_osborn_grid = - turbulent_diffusivity_Osborn_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_osborn_grid = oxygen_flux_osborn_grid*86400*(1000/eps_density_grid)        
        
        Gamma_BB_eps_grid = thesis.BB(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_BB_grid = Gamma_BB_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity    
        turbulent_diffusivity_BB_grid[turbulent_diffusivity_BB_grid<0] = np.nan
        oxygen_flux_BB_grid = - turbulent_diffusivity_BB_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_BB_grid = oxygen_flux_BB_grid*86400*(1000/eps_density_grid)
        
        Gamma_Skif_eps_grid = thesis.Skif(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity
        turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid<0] = np.nan
        oxygen_flux_Skif_grid = - turbulent_diffusivity_Skif_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400*(1000/eps_density_grid)
        """
        
        oxygen_flux_osborn_grid = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_Skif_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)

        
        was_the_last_profile_removed = False

        list_of_short_profiles = []
        count_of_short_profiles = 0
                
        for profile in range(number_of_profiles):
            #check if the profile was stopped too early by comparing it to the predecessor and succesor. If yes, skipt it

            try:
                slope = (bathymetrie[profile]-bathymetrie[profile+1])    
                next_slope = (bathymetrie[profile]-bathymetrie[profile-1])
            except IndexError:
                try:
                    slope = (bathymetrie[profile]-bathymetrie[profile+1])    
                    next_slope = (bathymetrie[profile]-bathymetrie[profile+2])
                except IndexError:
                    slope = (bathymetrie[profile]-bathymetrie[profile-1])
                    next_slope = (bathymetrie[profile]-bathymetrie[profile-2])
                           
            #only remove a profile if the slope to the next and to the overnext point is too high and the last profile was not removed
            if abs(slope)>acceptable_slope and abs(next_slope)>acceptable_slope: 
                
                #as the short profiles were stopped earlier, their bathymetrie value has to be smaller
                if slope <= 0 and next_slope <= 0:
                    count_of_short_profiles +=1
                    list_of_short_profiles.append(profile)
                
            else:
                was_the_last_profile_removed = False
                
        
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0
        lon_without_outliers = []
        bathymetrie_without_outliers = []
        
        transect_oxygen_flux_statistic = []


        
        
        for profile in range(number_of_profiles):
        
            profile_count+=1
            from_index = int(list_of_BBL_range_indices[profile])     #np.argmin(abs(eps_pressure-67)) 
            to_index = int(list_of_bathymetrie_indices[profile])     #np.argmin(abs(eps_pressure-77))
            
            
            #if the current profile is too short, skip it
            if profile in list_of_short_profiles:
                print(str(lon[profile])+": short profile")
                continue
            
            index_profile_before = profile-1
            index_profile_after = profile+1
            
            #assumes no two bad profiles next to each other
            if index_profile_before in list_of_short_profiles:
                index_profile_before -= 1
                
            if index_profile_after in list_of_short_profiles:
                index_profile_after += 1       
            
            #TODO
            #here only botched edge case handling
            if (index_profile_before < 0) or (index_profile_after >= number_of_profiles):
                continue
                
                
            profile_distance = geo.geodesic((lat[index_profile_before],lon[index_profile_before]),(lat[index_profile_after],lon[index_profile_after])).km*1000
               
            #slope in percent, TODO assumes a 1:1 ratio between dbar and meter    
            slope =  (bathymetrie[index_profile_after] - bathymetrie[index_profile_before]) /profile_distance
            
            #If the transect is from east to west, flip the sign of the slope
            if (lon[index_profile_before] < lon[index_profile_after]):
                #pass
                slope = -1*slope   
            
            
            #slope = np.abs(slope)
                
            #print(np.round(lon[profile],3),transect_index),np.round(np.nanmedian(np.log10(eps_grid[profile,30:-30])),3),np.shape(np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))))
                       

            #if the profile contains only nan values, profile is skipped
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,:])): #from_index:to_index
                print("NaN profile")
                continue

                               
            #check for an outlier profile by looking at unusual high dissipation troughout the whole watercolumn
            if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+2*spread_of_profile_medians):      
                #print("\toutlier")
                outlier_count += 1
                continue
           

                
            """
            #if the water colum portion contains only nan values, save only the bathymetrie then skip it
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,from_index:to_index])):
                #if the list is empty
                if np.any(bathymetrie_statistic[transect_index]) == None:
                    bathymetrie_statistic[transect_index] = [bathymetrie[profile]]            
                else:
                    #concatenate all further profiles to the ones already in the array
                    bathymetrie_statistic[transect_index].append(bathymetrie[profile])
                continue
            """
            
            
            
            print(eps_pressure[from_index],eps_pressure[to_index],from_index,to_index)

            
            max_flux = np.nanmax(oxygen_flux_BB_grid[profile,from_index:to_index])
            
            #if the flux is not reasonable, that means too high, replace it with the highest flux below the threshold
            if max_flux>maximum_reasonable_flux:
                temp_array = oxygen_flux_BB_grid[profile,from_index:to_index]
                temp_array[temp_array>maximum_reasonable_flux] = np.nan
                temp_array = temp_array[~np.isnan(temp_array)]
                max_flux = np.max(temp_array)

            if max_flux<0:
               max_flux = 0
                
            #same but for the minimum flux
            min_flux = np.nanmin(oxygen_flux_BB_grid[profile,from_index:to_index])              
            if min_flux< (-maximum_reasonable_flux):
                temp_array = oxygen_flux_BB_grid[profile,from_index:to_index]
                temp_array[temp_array<(-maximum_reasonable_flux)] = np.nan
                temp_array = temp_array[~np.isnan(temp_array)]
                min_flux = np.min(temp_array)         
            
            if min_flux>0:
               min_flux = 0
               
            min_max_array = np.asarray([[min_flux,max_flux]])
            #print(min_max_array)
            
            if len(transect_oxygen_flux_statistic) ==0:
                transect_oxygen_flux_statistic = min_max_array   
            else:
                transect_oxygen_flux_statistic = np.concatenate((transect_oxygen_flux_statistic,min_max_array),axis=0)
               
            #lon_without_outliers.append(lon[profile])
            #bathymetrie_without_outliers.append(bathymetrie[profile])
            
            #all_dissipation_statistic = np.append(all_dissipation_statistic,eps_grid[profile,from_index:to_index])
                            
            #Sort the profiles into the intervals
            #if the list is empty
            if np.any(dissipation_statistic[transect_index]) == None:
  
                #fill it with a reshaped profile
                oxygen_flux_statistic[transect_index] = min_max_array
                dissipation_statistic[transect_index] = [np.nanmean(eps_grid[profile,from_index:to_index])]
                #dissipation_statistic[transect_index] = [eps_grid[profile,from_index:to_index]]
                slope_statistic[transect_index] = [slope] 
                lon_statistic[transect_index] = [lon[profile]] 
                bathymetrie_statistic[transect_index] = [bathymetrie[profile]]
                                                 
            else:
                #concatenate all further profiles to the ones already in the array
                oxygen_flux_statistic[transect_index] = np.concatenate((oxygen_flux_statistic[transect_index],min_max_array),axis=0)
                dissipation_statistic[transect_index].append(np.nanmean(eps_grid[profile,from_index:to_index]))
                #dissipation_statistic[transect_index] = np.concatenate((dissipation_statistic[transect_index],eps_grid[profile,from_index:to_index]),axis=0)
                slope_statistic[transect_index].append(slope)
                lon_statistic[transect_index].append(lon[profile])               
                bathymetrie_statistic[transect_index].append(bathymetrie[profile])
                    


        print("removed",outlier_count,"profiles as outliers")
        print("removed",count_of_short_profiles,"profiles as they did not reach the sea floor")


    
    
    print("\n\n\n Results: \n\n")
    for index in range(number_of_transects):
        print(index,np.shape(slope_statistic[index]),np.shape(dissipation_statistic[index]))
    
    
    f_dissip,dissip_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
    
    for index in range(number_of_transects):
    
        #check if container from that transect contains no data
        if slope_statistic[index] == None:
            continue
        
        #point_positions = np.mean()
        #errorbars_up = 
        #error
        #dissip_axarr.errorbar()
        #dissip_axarr.semilogy(slope_statistic[index],dissipation_statistic[index],".", label = list_of_transect_names[index])
        #dissip_axarr.semilogy(lon_statistic[index],dissipation_statistic[index], label = list_of_transect_names[index])
        dissip_axarr.plot(lon_statistic[index],np.asarray(bathymetrie_statistic[index]),"k.", label = list_of_transect_names[index])
        #break
    
  
    """           
    dissip_axarr.set_ylabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
    dissip_axarr.set_xlabel("slope [%]")             
    dissip_axarr.legend(loc = "upper left")
    """
    
    #f_dissip.suptitle(cruisename+": mean dissipation in the lowermost "+str(height_above_ground)+" meters above ground")
    #f_dissip.suptitle(cruisename+": mean dissipation around the halocline (67-77dbar) ("+str(number_of_transects)+" intervals)")
    
    f_dissip.suptitle("Depth of all valid MSS profiles measured during "+cruisename)
    dissip_axarr.invert_yaxis()
    dissip_axarr.set_ylabel(r"pressure [dbar]")   
    dissip_axarr.set_xlabel(r'Longitude [$^\circ$E]')
    
    f_dissip.set_size_inches(7.2,7.2/1.61)
            
    #f_dissip.set_size_inches(18,10.5)
    f_dissip.tight_layout()
    f_dissip.subplots_adjust(top=0.925)
    #f_dissip.subplots_adjust(top=0.95)
    #f_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_transects)+"_intervals_mean_dissipation")
    #f_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+"slope_dissipation")
        
    
    f_check,check_axarr = plt.subplots(nrows = 2, ncols = 1, sharex = True) 
    for index in range(number_of_transects):
    
        #check if container from that transect contains no data
        if slope_statistic[index] == None:
            continue
        
        #point_positions = np.mean()
        #errorbars_up = 
        #error
        #dissip_axarr.errorbar()
        check_axarr[0].plot(lon_statistic[index],slope_statistic[index],label = list_of_transect_names[index])
        
        cumsum = np.cumsum(slope_statistic[index])
        if (lon_statistic[index][-1] > lon_statistic[index][0]):
            cumsum = - np.cumsum(slope_statistic[index])    
        check_axarr[1].plot(lon_statistic[index],cumsum-min(cumsum),label = list_of_transect_names[index])
        #break
    

    check_axarr[0].legend()    
    check_axarr[1].invert_yaxis()
    plt.show()
    
   
