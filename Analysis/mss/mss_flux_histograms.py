############################################################
#this program loads all profiles from a cruise, removes outliers,
#retireves the fluxes (all 3 parametrizations)in a given density or pressure intervall in choosable longitude intervals
#and plots them as histograms

#Enables to compare fluxes in the interior and on the slope and parametrizations

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

import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
    
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177"]

averaging_intervals_borders = [20.55,20.625]

flux_through_halocline = True #set to True if the flux trough the Halocline instead of the BBL should be computed 
height_above_ground = 10 #Size of the averaging interval above ground for the BBL, has no meaning if (flux_through_halocline == True)
maximum_reasonable_flux = 200 #Fluxes above this value will be discarded
acceptable_slope = 20 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 

flux_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_flux_percentile = 97.72
dissip_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_dissip_percentile = 97.72

number_of_dissipation_subplots = 2 #Decide if both the mean and the median subplots is shown or only the mean

number_of_intervals = len(averaging_intervals_borders)+1  
 
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    number_of_fluxes_over_the_threshold = 0
    total_number_of_fluxes = 0
        
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print(cruisename)
    print(averaging_intervals_borders)
             
    #2 borders = 3 intervals
    oxygen_flux_BB_statistic = [None] * number_of_intervals
    oxygen_flux_Shih_statistic = [None] * number_of_intervals
    oxygen_flux_Osborn_statistic = [None] * number_of_intervals

    bathymetrie_statistic = [None] * number_of_intervals
    median_dissipation_statistic = [None] * number_of_intervals
    mean_dissipation_statistic = [None] * number_of_intervals
    all_dissipation_statistic = np.asarray([])
    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".npz":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    #print(DATAFILENAMES)
    
    profile_count = 0
    
    for DATAFILENAME in DATAFILENAMES:
    
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        transect_name = DATAFILENAME[:-4]
    
        #skip the short "S206" transects
        if transect_name[0:4] == "S106":
            print(transect_name,"skipped")
            continue
            
        print("\n",transect_name)
            
        
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
        
        oxygen_flux_Osborn_grid = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_Shih_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0
        lon_without_outliers = []
        bathymetrie_without_outliers = []
        
        transect_oxygen_flux_statistic = []
        

        count_of_short_profiles = 0
        
        for profile in range(number_of_profiles):
        
            profile_count+=1
            from_index = int(list_of_BBL_range_indices[profile]) 
            to_index = int(list_of_bathymetrie_indices[profile])
            
            if flux_through_halocline == True:
                from_index =  np.argmin(abs(eps_pressure-67))     
                to_index = np.argmin(abs(eps_pressure-77))
            
            for index,border in enumerate(averaging_intervals_borders):
                #print("index",index)
                
                #Sort the prifle into the intervals
                if lon[profile]  <= border:
                    interval_number = index
                    break
                    
                #longitude of the profile was greater as the last border --> last interval
                elif lon[profile] > averaging_intervals_borders[-1]:
                    interval_number =  len(averaging_intervals_borders)             
                
            #print(np.round(lon[profile],3),interval_number),np.round(np.nanmedian(np.log10(eps_grid[profile,30:-30])),3),np.shape(np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))))
                       

            #if the profile contains only nan values, profile is skipped
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,:])): #from_index:to_index
                print("NaN profile")
                continue

                               
            #check for an outlier profile, ergo too high dissipation rates compared with the surrounding
            if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+2*spread_of_profile_medians):      
                #print("\toutlier")
                outlier_count += 1
                continue
          
                
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
                    print("removed a short profile")
                    continue
                
            #if the water colum portion contains only nan values, save only the bathymetrie then skip it
            #useful if the interval to average over is deper than the current bathymetrie
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,from_index:to_index])):
                #if the list is empty
                if np.any(bathymetrie_statistic[interval_number]) == None:
                    bathymetrie_statistic[interval_number] = [bathymetrie[profile]]            
                else:
                    #concatenate all further profiles to the ones already in the array
                    bathymetrie_statistic[interval_number].append(bathymetrie[profile])
                continue
   

            total_number_of_fluxes += oxygen_flux_BB_grid[profile,from_index:to_index].size

                            
            #Sort the profiles into the intervals
            #if the list is empty
            if np.any(oxygen_flux_BB_statistic[interval_number]) == None:
  
                #fill it with the fluxes from the interval
                oxygen_flux_Osborn_statistic[interval_number] = oxygen_flux_Osborn_grid[profile,from_index:to_index]
                oxygen_flux_Shih_statistic[interval_number] = oxygen_flux_Shih_grid[profile,from_index:to_index]
                oxygen_flux_BB_statistic[interval_number] = oxygen_flux_BB_grid[profile,from_index:to_index]
                 
                median_dissipation_statistic[interval_number] = [np.nanmedian(eps_grid[profile,from_index:to_index])]
                mean_dissipation_statistic[interval_number] = [np.nanmean(eps_grid[profile,from_index:to_index])]                
            else:
                #concatenate all further profiles to the ones already in the array
                oxygen_flux_Osborn_statistic[interval_number] = np.append(oxygen_flux_Osborn_statistic[interval_number],oxygen_flux_Osborn_grid[profile,from_index:to_index])
                oxygen_flux_Shih_statistic[interval_number] = np.append(oxygen_flux_Shih_statistic[interval_number],oxygen_flux_Shih_grid[profile,from_index:to_index])
                oxygen_flux_BB_statistic[interval_number] = np.append(oxygen_flux_BB_statistic[interval_number],oxygen_flux_BB_grid[profile,from_index:to_index])
                                                
                median_dissipation_statistic[interval_number].append(np.nanmedian(eps_grid[profile,from_index:to_index]))
                mean_dissipation_statistic[interval_number].append(np.nanmean(eps_grid[profile,from_index:to_index]))

            #sort the bathymetrie data into the intervals
            if np.any(bathymetrie_statistic[interval_number]) == None:
                bathymetrie_statistic[interval_number] = [bathymetrie[profile]]            
            else:
                #concatenate all further profiles to the ones already in the array
                bathymetrie_statistic[interval_number].append(bathymetrie[profile])
                    


        print("removed",outlier_count,"profiles as outliers")
        print("removed",count_of_short_profiles,"profiles as they did not reach the sea floor")
        
        
        
        
    ###############################################################################################
    ###############################################################################################
    ###############################################################################################
    ###############################################################################################
    #Plotting of the maximum flux values per transect
    ###############################################################################################
    ###############################################################################################
    ###############################################################################################
    ###############################################################################################
    f_hist, hist_axarr = plt.subplots(nrows = 3, ncols = 2, sharex = True)

    bins = np.arange(-200,5,2.5)
    hist_axarr[0,0].hist(oxygen_flux_Osborn_statistic[0],bins = bins,edgecolor='black', linewidth=1.2, log = True)
    hist_axarr[0,1].hist(oxygen_flux_Osborn_statistic[1],bins = bins,edgecolor='black', linewidth=1.2, log = True)
    hist_axarr[1,0].hist(oxygen_flux_Shih_statistic[0],bins = bins,edgecolor='black', linewidth=1.2, log = True)
    hist_axarr[1,1].hist(oxygen_flux_Shih_statistic[1],bins = bins,edgecolor='black', linewidth=1.2, log = True)
    hist_axarr[2,0].hist(oxygen_flux_BB_statistic[0],bins = bins,edgecolor='black', linewidth=1.2, log = True)
    hist_axarr[2,1].hist(oxygen_flux_BB_statistic[1],bins = bins,edgecolor='black', linewidth=1.2, log = True)                
    
    """
    hist_axarr[0,0].axvline(np.nanmean(oxygen_flux_Osborn_statistic[0]), linewidth=2.4, c = "tab:red")
    hist_axarr[0,1].axvline(np.nanmean(oxygen_flux_Osborn_statistic[1]), linewidth=2.4, c = "tab:red")
    hist_axarr[1,0].axvline(np.nanmean(oxygen_flux_Shih_statistic[0]), linewidth=2.4, c = "tab:red")
    hist_axarr[1,1].axvline(np.nanmean(oxygen_flux_Shih_statistic[1]), linewidth=2.4, c = "tab:red")
    hist_axarr[2,0].axvline(np.nanmean(oxygen_flux_BB_statistic[0]), linewidth=2.4, c = "tab:red")
    hist_axarr[2,1].axvline(np.nanmean(oxygen_flux_BB_statistic[1]), linewidth=2.4, c = "tab:red")
    """
        
    hist_axarr[2,0].set_xlabel(r"oxygen flux [mmol/(m$^2$*d]")
    hist_axarr[2,1].set_xlabel(r"oxygen flux [mmol/(m$^2$*d]")
    hist_axarr[0,0].set_ylabel("frequency [#]")
    hist_axarr[1,0].set_ylabel("frequency [#]")
    hist_axarr[2,0].set_ylabel("frequency [#]")
    hist_axarr[0,0].set_title("Interior, Osborn flux")
    hist_axarr[0,1].set_title("Slope, Osborn flux")
    hist_axarr[1,0].set_title("Interior, Shih flux")
    hist_axarr[1,1].set_title("Slope, Shih flux")
    hist_axarr[2,0].set_title("Interior, BB flux")
    hist_axarr[2,1].set_title("Slope, BB flux")
    
    print("interior,osborn",oxygen_flux_Osborn_statistic[0].size)
    print("slope,osborn",oxygen_flux_Osborn_statistic[1].size)
    print("interior,Shih",oxygen_flux_Shih_statistic[0].size)
    print("slope,Shih",oxygen_flux_Shih_statistic[1].size)
    print("interior,BB",oxygen_flux_BB_statistic[0].size)
    print("slope,BB",oxygen_flux_BB_statistic[1].size)
                                            
    f_hist.suptitle(cruisename+": flux distribution around the halocline (67-77dbar)")
    f_hist.set_size_inches(18,10.5)
    f_hist.tight_layout() 
    f_hist.subplots_adjust(top=0.9)

    plt.show()
