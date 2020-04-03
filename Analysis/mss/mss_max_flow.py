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
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
    
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169","/home/ole/share-windows/processed_mss/emb177"]

#averaging_intervals_borders = [20.55,20.62]
averaging_intervals_borders = np.linspace(20.48,20.7,44)
height_above_ground = 10
maximum_reasonable_flux = 150
acceptable_slope = 20 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 
flux_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_flux_percentile = 97.72
dissip_percentile = 84.13 #percentile which is displayed as the error bar (variable spread)
second_dissip_percentile = 97.72

number_of_intervals = len(averaging_intervals_borders)+1   
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
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
    oxygen_flux_statistic = [None] * number_of_intervals

    bathymetrie_statistic = [None] * number_of_intervals
    dissipation_statistic = [None] * number_of_intervals
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
            
            
        #boundary_check_grid = np.zeros(np.shape(ozmidov_scale_grid))
        #check if ozimidov scale is bigger than the distance from the ground
        #for i in range(number_of_profiles):
        #density_grid[ozmidov_scale_grid,:])))
        
        Gamma_Osborn_eps_grid = thesis.Osborn(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity 
        turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid<0] = np.nan
        oxygen_flux_osborn_grid = turbulent_diffusivity_Osborn_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_osborn_grid = oxygen_flux_osborn_grid*86400*(1000/eps_density_grid)        
        
        Gamma_BB_eps_grid = thesis.BB(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_BB_grid = Gamma_BB_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity    
        turbulent_diffusivity_BB_grid[turbulent_diffusivity_BB_grid<0] = np.nan
        oxygen_flux_BB_grid = turbulent_diffusivity_BB_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_BB_grid = oxygen_flux_BB_grid*86400*(1000/eps_density_grid)
        
        Gamma_Skif_eps_grid = thesis.Skif(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity
        turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid<0] = np.nan
        oxygen_flux_Skif_grid = turbulent_diffusivity_Skif_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400*(1000/eps_density_grid)
        
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0
        lon_without_outliers = []
        bathymetrie_without_outliers = []
        
        transect_oxygen_flux_statistic = []
        
        was_the_last_profile_removed = False
        count_of_short_profiles = 0
        
        for profile in range(number_of_profiles):
        
            profile_count+=1
            from_index = np.argmin(abs(eps_pressure-67))        #int(list_of_BBL_range_indices[profile])     
            to_index = np.argmin(abs(eps_pressure-77))          #int(list_of_bathymetrie_indices[profile])
            
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

                               
            #check for an outlier profile 
            if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+2*spread_of_profile_medians):      
                #print("\toutlier")
                outlier_count += 1
                continue
           
            #check if the profile was stopped too early by comparing it to the predecessor and succesor. If yes, skipt it
            try:
                slope = (bathymetrie[profile]-bathymetrie[profile+1])    
                next_slope = (bathymetrie[profile]-bathymetrie[profile+2])
            except IndexError:
                slope = (bathymetrie[profile]-bathymetrie[profile-1])
                next_slope = (bathymetrie[profile]-bathymetrie[profile-2])
                           
            #only remove a profile if the slope to the next and to the overnext point is too high and the last profile was not removed
            if abs(slope)>acceptable_slope and abs(next_slope)>acceptable_slope and not was_the_last_profile_removed:
                was_the_last_profile_removed = True
                count_of_short_profiles +=1
                print("removed a short profile")
                continue
            else:
                was_the_last_profile_removed = False
                

            #if the water colum portion contains only nan values, save only the bathymetrie then skip it
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,from_index:to_index])):
                #if the list is empty
                if np.any(bathymetrie_statistic[interval_number]) == None:
                    bathymetrie_statistic[interval_number] = [bathymetrie[profile]]            
                else:
                    #concatenate all further profiles to the ones already in the array
                    bathymetrie_statistic[interval_number].append(bathymetrie[profile])
                continue
   
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
               
            lon_without_outliers.append(lon[profile])
            bathymetrie_without_outliers.append(bathymetrie[profile])
            
            all_dissipation_statistic = np.append(all_dissipation_statistic,eps_grid[profile,from_index:to_index])
                            
            #Sort the profiles into the intervals
            #if the list is empty
            if np.any(oxygen_flux_statistic[interval_number]) == None:
  
                #fill it with a reshaped profile
                oxygen_flux_statistic[interval_number] = min_max_array
                dissipation_statistic[interval_number] = [np.nanmean(eps_grid[profile,from_index:to_index])]
                                
            else:
                #concatenate all further profiles to the ones already in the array
                oxygen_flux_statistic[interval_number] = np.concatenate((oxygen_flux_statistic[interval_number],min_max_array),axis=0)
                dissipation_statistic[interval_number].append(np.nanmean(eps_grid[profile,from_index:to_index]))

            #sort the bathymetrie data into the intervals
            if np.any(bathymetrie_statistic[interval_number]) == None:
                bathymetrie_statistic[interval_number] = [bathymetrie[profile]]            
            else:
                #concatenate all further profiles to the ones already in the array
                bathymetrie_statistic[interval_number].append(bathymetrie[profile])
                    


        print("removed",outlier_count,"profiles as outliers")
        print("removed",count_of_short_profiles,"profiles as they did not reach the sea floor")

        ###############################################################################################
        #Plotting of the maximum flux values per transect
        ###############################################################################################
        f1, axarr1 = plt.subplots(nrows = 1, ncols = 1, sharex = True)
        
        bathymetrie_axes = axarr1.twinx()
        bathymetrie_axes.set_ylim((min(bathymetrie_without_outliers)-5,max(bathymetrie_without_outliers)))
        bathymetrie_axes.invert_yaxis()
        bathymetrie_axes.fill_between(lon_without_outliers,bathymetrie_without_outliers, np.ones(len(lon_without_outliers))*max(bathymetrie_without_outliers),color = "lightgrey", alpha = 0.8)
    
        axarr1.set_xlabel("longitude [degree]")    
        axarr1.set_ylabel(r"BB oxygen flux [mmol/(m$^2$*d]")
        bathymetrie_axes.set_ylabel("pressure [dbar]")
        #plot maximum flux
        axarr1.plot(lon_without_outliers,transect_oxygen_flux_statistic[:,1])
        
        #plot minimum flux
        axarr1.plot(lon_without_outliers,transect_oxygen_flux_statistic[:,0])
        axarr1.set_xlabel("longitude")        
                   
        f1.set_size_inches(9,5)
        f1.tight_layout()
        f1.subplots_adjust(top=0.95)
        f1.suptitle("max oxygen flux "+cruisename+" "+transect_name+" "+str(len(lon_without_outliers))+" profiles")
             
        f1.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/max_flux_transects/"+cruisename+"_"+transect_name+"_halocline", DPI = 300)         
        #plt.show()      
        plt.close(fig = "all")
        
        
    ###########################################################################################################################
    
    print("\n\nSorting result:")
    for index,interval in enumerate(dissipation_statistic):
        if index == 0:
            title = "-"+str(np.round(averaging_intervals_borders[0],2))
        elif index == number_of_intervals-1:
            title = str(np.round(averaging_intervals_borders[-1],2))+"-"
        else:
            title = str(np.round(averaging_intervals_borders[index-1],2))+"-"+str(np.round(averaging_intervals_borders[index],2))  
            
        print(title,np.shape(interval),"profiles")
            
    
    #compute mean and std over the saved intervals
    mean_max_flux = [None] * number_of_intervals
    median_max_flux = [None] * number_of_intervals
    upper_percentile_max_flux = [None] * number_of_intervals
    lower_percentile_max_flux = [None] * number_of_intervals  
      
    mean_min_flux = [None] * number_of_intervals
    median_min_flux = [None] * number_of_intervals
    upper_percentile_min_flux = [None] * number_of_intervals
    lower_percentile_min_flux = [None] * number_of_intervals  
    second_upper_percentile_max_flux = [None] * number_of_intervals
    second_lower_percentile_min_flux = [None] * number_of_intervals  
        
    bathymetrie_mean = [None] * number_of_intervals
    
    mean_dissipation = [None] * number_of_intervals
    median_dissipation = [None] * number_of_intervals
    std_dissipation = [None] * number_of_intervals
    lower_percentile_dissip = [None] * number_of_intervals
    upper_percentile_dissip = [None] * number_of_intervals

    second_lower_percentile_dissip = [None] * number_of_intervals
    second_upper_percentile_dissip = [None] * number_of_intervals

    for index in range(number_of_intervals):
        #print(np.shape(bathymetrie_statistic[index]))
        try:
            bathymetrie_mean[index] = np.nanmean(bathymetrie_statistic[index])
        except (TypeError,AttributeError):
            #print("error")
            bathymetrie_mean[index] = np.nan
                       
                       
                           
    for index in range(number_of_intervals):
        try:
            mean_max_flux[index] = np.nanmean(oxygen_flux_statistic[index][:,1],axis=0)
            median_max_flux[index] = np.nanmedian(oxygen_flux_statistic[index][:,1],axis=0)
            
            upper_percentile_max_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,1], flux_percentile)
            lower_percentile_max_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,1], 100-flux_percentile)
            
            mean_min_flux[index] = np.nanmean(oxygen_flux_statistic[index][:,0],axis=0)
            median_min_flux[index] = np.nanmedian(oxygen_flux_statistic[index][:,0],axis=0)

            upper_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], flux_percentile)
            lower_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], 100-flux_percentile)
            
            second_upper_percentile_max_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,1], second_flux_percentile)
            second_lower_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], 100-second_flux_percentile)
                    
            #bathymetrie_mean[index] = np.nanmean(bathymetrie_statistic[index])
            
            mean_dissipation[index] = np.nanmean(np.log10(dissipation_statistic[index]),axis=0)
            median_dissipation[index] = np.nanmedian(np.log10(dissipation_statistic[index]),axis=0)
            std_dissipation[index] = np.nanstd(np.log10(dissipation_statistic[index]),axis=0) 
            upper_percentile_dissip[index] = np.nanpercentile(np.log10(dissipation_statistic[index]), dissip_percentile)
            lower_percentile_dissip[index] = np.nanpercentile(np.log10(dissipation_statistic[index]), 100-dissip_percentile)
            
            second_upper_percentile_dissip[index] = np.nanpercentile(np.log10(dissipation_statistic[index]), second_dissip_percentile)
            second_lower_percentile_dissip[index] = np.nanpercentile(np.log10(dissipation_statistic[index]), 100-second_dissip_percentile)
        except (TypeError,AttributeError):
            mean_max_flux[index] = np.nan
            median_max_flux[index] = np.nan
            
            upper_percentile_max_flux[index] = np.nan
            lower_percentile_max_flux[index] = np.nan
            
            mean_min_flux[index] = np.nan
            median_min_flux[index] = np.nan

            upper_percentile_min_flux[index] = np.nan
            lower_percentile_min_flux[index] = np.nan
            
            second_upper_percentile_max_flux[index] = np.nan
            second_lower_percentile_min_flux[index] = np.nan
                    
            #bathymetrie_mean[index] = np.nan
            
            mean_dissipation[index] = np.nan
            median_dissipation[index] = np.nan
            std_dissipation[index] = np.nan
            upper_percentile_dissip[index] = np.nan
            lower_percentile_dissip[index] = np.nan
            
            second_upper_percentile_dissip[index] = np.nan
            second_lower_percentile_dissip[index] = np.nan
                       
    mean_max_flux = np.asarray(mean_max_flux)
    median_max_flux = np.asarray(median_max_flux)  
    
    mean_dissipation = np.asarray(mean_dissipation)
    std_dissipation = np.asarray(std_dissipation)
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
        
        
       
    f_flux,flux_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
    bathymetrie_axes = flux_axarr.twinx()
    
    #define flux_axarr as the foreground
    flux_axarr.set_zorder(10)
    flux_axarr.patch.set_visible(False)
        
    #append one extra border behind the last border in the mean distance of the borders 
    plot_longitude = np.append(averaging_intervals_borders,averaging_intervals_borders[-1]+np.mean(np.diff(averaging_intervals_borders)))
    #shift all plot points by half the border distance  
    plot_longitude - np.mean(np.diff(averaging_intervals_borders))/2
                      
    flux_axarr.plot(plot_longitude,mean_max_flux, zorder = 3,c = "k")#"tab:blue")
    flux_axarr.plot(plot_longitude,median_max_flux,ls = "--",c = "k")#"tab:blue")
    
    flux_axarr.fill_between(plot_longitude,upper_percentile_max_flux,lower_percentile_max_flux, color = "tab:blue", zorder = 2, alpha = 0.7)
    flux_axarr.fill_between(plot_longitude,upper_percentile_max_flux,second_upper_percentile_max_flux, color = "tab:blue", zorder = 2, alpha = 0.4)
    
    flux_axarr.plot(plot_longitude,mean_min_flux, color = "k")#"tab:green", zorder = 3)
    flux_axarr.plot(plot_longitude,median_min_flux,ls = "--",c = "k")#"tab:green")
    flux_axarr.fill_between(plot_longitude,upper_percentile_min_flux,lower_percentile_min_flux, color = "tab:green", zorder = 2, alpha = 0.6)   
    flux_axarr.fill_between(plot_longitude,second_lower_percentile_min_flux,lower_percentile_min_flux, color = "tab:green", zorder = 2, alpha = 0.4)  
    
    bathymetrie_axes.set_ylim((np.nanmin(bathymetrie_mean)-5,np.nanmax(bathymetrie_mean)))
    bathymetrie_axes.invert_yaxis()
    bathymetrie_axes.set_ylabel("pressure [dbar]")
        
    bathymetrie_axes.fill_between(plot_longitude,bathymetrie_mean, np.ones(len(bathymetrie_mean))*max(bathymetrie_mean),color = "lightgrey", zorder = 1, alpha = 0.8)
    
    mean_label = mlines.Line2D([], [], color='k', label='mean')
    median_label = mlines.Line2D([], [], ls = "--", color='k', label='median')
    down_label = mlines.Line2D([], [], color='tab:blue', label='maximum downwards BB flux')
    up_label = mlines.Line2D([], [], color='tab:green', label='maximum upwards BB flux')
    
    max_flux_label =  mpatches.Patch(color='tab:blue', alpha = 0.6,label='downwards flux '+str(flux_percentile)+"% percentile")
    min_flux_label =  mpatches.Patch(color='tab:green', alpha = 0.6, label='upwards flux '+str(flux_percentile)+"% percentile")
    second_max_flux_label =  mpatches.Patch(color='tab:blue', alpha = 0.4,label='downwards flux '+str(second_flux_percentile)+"% percentile")
    second_min_flux_label =  mpatches.Patch(color='tab:green', alpha = 0.4, label='upwards flux '+str(second_flux_percentile)+"% percentile")
       
    bathymetrie_label =  mpatches.Patch(color='lightgrey', label='bathymetrie')
    flux_axarr.legend(handles=[mean_label,median_label,down_label,up_label,max_flux_label, second_max_flux_label, min_flux_label, second_min_flux_label,bathymetrie_label]) #loc=8
    
    #flux_axarr.set_ylim((-55,55))
    flux_axarr.set_xlabel("longitude [degree]")    
    flux_axarr.set_ylabel(r"oxygen flux [mmol/(m$^2$*d]")
    #f_flux.suptitle(cruisename+": maximum up- and downwards BB oxygen flux in the lowermost "+str(height_above_ground)+" meters above ground")
    f_flux.suptitle(cruisename+": maximum up- and downwards BB oxygen flux around the halocline (67-77dbar) ("+str(number_of_intervals)+" intervals)")
    
    f_flux.set_size_inches(18,10.5)
    f_flux.tight_layout() 
    f_flux.subplots_adjust(top=0.95)
    #f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"_intervals_max_oxygen_flux")
    f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"_intervals_halocline_flux")
        
    ###############################################################################################################
    
    f_dissip,dissip_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
     #define dissip_axarr as the foreground
    dissip_axarr.set_zorder(10)
    dissip_axarr.patch.set_visible(False)
    
    bathymetrie_axes2 = dissip_axarr.twinx()
    bathymetrie_axes2.set_ylim((np.nanmin(bathymetrie_mean)-5,np.nanmax(bathymetrie_mean)))
    bathymetrie_axes2.invert_yaxis()
    bathymetrie_axes2.set_ylabel("pressure [dbar]")
        
    bathymetrie_axes2.fill_between(plot_longitude,bathymetrie_mean, np.ones(len(bathymetrie_mean))*max(bathymetrie_mean),color = "lightgrey", zorder = 1, label = "bathymetrie")
    
    dissip_axarr.plot(plot_longitude,mean_dissipation, "k", label ="mean dissipation")
    dissip_axarr.plot(plot_longitude,median_dissipation, "k--", label ="median dissipation")
    #dissip_axarr.fill_between(plot_longitude,mean_dissipation-std_dissipation,mean_dissipation+std_dissipation, alpha = 0.5)
    dissip_axarr.fill_between(plot_longitude,upper_percentile_dissip,lower_percentile_dissip, color = "tab:blue", alpha = 0.7)  
    dissip_axarr.fill_between(plot_longitude,upper_percentile_dissip,second_upper_percentile_dissip, color = "tab:blue", alpha = 0.4) 
    dissip_axarr.fill_between(plot_longitude,lower_percentile_dissip,second_lower_percentile_dissip, color = "tab:blue", alpha = 0.4) 
                
    dissip_axarr.set_ylabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")   
    dissip_axarr.set_xlabel("longitude [degree]")             
    bathymetrie_label =  mpatches.Patch(color='lightgrey', label='bathymetrie')
    dissip_mean_label = mlines.Line2D([], [], color= "k", label='mean dissipation $\epsilon$') #tab:blue
    dissip_median_label = mlines.Line2D([], [], color = "k", ls = "--", label='median dissipation $\epsilon$')
    dissip_percent_label =  mpatches.Patch(color='tab:blue', label=str(dissip_percentile)+"% percentile")
    dissip_second_percent_label =  mpatches.Patch(color='tab:blue', alpha = 0.4, label=str(second_dissip_percentile)+"% percentile")           
    dissip_axarr.legend(handles=[dissip_mean_label,dissip_median_label,dissip_percent_label,dissip_second_percent_label,bathymetrie_label])
    
    #f_dissip.suptitle(cruisename+": mean dissipation in the lowermost "+str(height_above_ground)+" meters above ground")
    f_dissip.suptitle(cruisename+": mean dissipation around the halocline (67-77dbar) ("+str(number_of_intervals)+" intervals)")

    f_dissip.set_size_inches(18,10.5)
    f_dissip.tight_layout() 
    f_dissip.subplots_adjust(top=0.95)
    #f_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"_intervals_mean_dissipation")
    f_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"_intervals_halocline_mean_dissipation")
    
    """
    hist_dissip,hist_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
    
    filtered_list = list(filter(None, dissipation_statistic))
    flat_list = [item for sublist in filtered_list for item in sublist]
    flattened_array = np.asarray(flat_list).flatten()
    
    print("profile_count",profile_count)
    print(np.shape(all_dissipation_statistic))
    hist_axarr.hist(np.log10(all_dissipation_statistic),bins = 40,edgecolor='black', linewidth=1.2)
    hist_axarr.axvline(np.nanmean(np.log10(all_dissipation_statistic)), linewidth=2.4, c = "tab:red", label = "mean of logs")
    hist_axarr.axvline(np.log10(np.nanmean(all_dissipation_statistic)), linewidth=2.4, c = "tab:green", label = "log of means")
    hist_axarr.axvline(np.nanmedian(np.log10(all_dissipation_statistic)), linewidth=2.4, ls = "-", c = "tab:orange", label = "median of logs")
    hist_axarr.axvline(np.log10(np.nanmedian(all_dissipation_statistic)), linewidth=2.4, ls = ":", c = "k", label = "log of median")
    
    hist_axarr.legend()
    
    hist_axarr.set_xlabel(r"log10($\epsilon$) $[m^2 s^{-3}]$") 
    hist_axarr.set_ylabel("frequency [#]")
                
    hist_dissip.suptitle(cruisename+": dissipation distribution around the halocline (67-77dbar)")
    hist_dissip.set_size_inches(18,10.5)
    hist_dissip.tight_layout() 
    hist_dissip.subplots_adjust(top=0.95)
    hist_dissip.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_dissipation_hist_halocline")
    """
                
    plt.show()
    
    
    
    
    
