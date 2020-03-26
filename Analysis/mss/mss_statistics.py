############################################################
#this program loads all profiles from a cruise, removes outliers
#and compute mean and standard deviation of profiles
#in choosable longitude intervals

##############################################################
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import warnings
warnings.filterwarnings('ignore') #ignoring warnings, specially numpy warnings

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)
    
    
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169","/home/ole/share-windows/processed_mss/emb177"]

averaging_intervals_borders = [20.55,20.62]
#averaging_intervals_borders = np.linspace(20.48,20.7,9)
number_of_intervals = len(averaging_intervals_borders)+1

for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print(cruisename)    
    #2 borders = 3 intervals
    print(averaging_intervals_borders)

    salinity_statistic = [None] * number_of_intervals
    temperature_statistic = [None] * number_of_intervals
    dissipation_statistic = [None] * number_of_intervals
    oxygen_sat_statistic = [None] * number_of_intervals
    oxygen_flux_statistic = [None] * number_of_intervals

        
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
    
        #define the pictures
        #f1, axarr1 = plt.subplots(nrows = 1, ncols = 7, sharey = True)
        #f2, axarr2 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
    
    

        
        
        data = np.load(datafile_path)
        
        number_of_profiles = data["number_of_profiles"] #
        lat = data["lat"] #Latitude of the profiles
        lon = data["lon"] #Longitude of the profiles
        distance = data["distance"] #distance from the starting profile (monotonically increasing)
        
        bathymetrie = data["bathymetrie"]
        list_of_bathymetrie_indices = data["list_of_bathymetrie_indices"]
        BBL = data["BBL"]
        list_of_BBL_indices = data["list_of_BBL_indices"]
        BBL_range = data["BBL_range"]
        list_of_BBL_range_indices = data["list_of_BBL_range_indices"]
        
        interp_pressure = data["interp_pressure"]
        oxygen_grid = data["oxygen_grid"]
        #oxygen_sat_grid = data["oxygen_sat_grid"]
        salinity_grid = data["salinity_grid"]
        consv_temperature_grid = data["consv_temperature_grid"]
        density_grid = data["density_grid"]
        
        eps_pressure = data["eps_pressure"]
        eps_grid = data["eps_grid"]
        corrected_eps_grid = data["corrected_eps_grid"]
        corrected_eps_wiki_grid = data["corrected_eps_wiki_grid"]
        eps_consv_temperature_grid = data["eps_consv_temperature_grid"]
        eps_oxygen_grid = data["eps_oxygen_grid"] 
        eps_oxygen_sat_grid = data["eps_oxygen_sat_grid"]         
        eps_salinity_grid = data["eps_salinity_grid"]
        
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
        eps_salinity_grid               absolute salinity in g/kg as a grid (number_of_profiles x len(eps_pressure)) 
        eps_oxygen_grid                 oxygen concentration in micromol per litre as a grid (number_of_profiles x len(eps_pressure))
        eps_N_squared_grid              N^2, the Brunt-Vaisala frequency in 1/s^2 as a grid (number_of_profiles x len(eps_pressure))
        eps_density_grid                density in kg/m^3 as a grid (number_of_profiles x len(eps_pressure))
        
        eps_viscosity_grid
        eps_Reynolds_bouyancy_grid
        corrected_eps_Reynolds_bouyancy_grid 
        eps_wiki_Reynolds_bouyancy_grid
        corrected_eps_wiki_Reynolds_bouyancy_grid 
        
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
        """
        
        print("Number of profiles:",number_of_profiles)
        
        print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
        
        eps_N_grid = np.sqrt(eps_N_squared_grid)
        #ozmidov scale
        ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))
        
        #conversion from pressure coordinates to depth
        eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
        #bathymetrie_in_m = gsw.z_from_p(bathymetrie,np.mean(lat))
        
        #eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
        
        #distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetrie_in_m,(-1,1))
        #boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)
            
        
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
        transect_median = np.nanmean(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        count_outliers = 0
        
        for profile in range(number_of_profiles):
            for index,border in enumerate(averaging_intervals_borders):
                #print("index",index)
                if lon[profile]  <= border:
                    interval_number = index
                    break
                    
                #longitude of the profile was greater as the last border --> last interval
                elif lon[profile] > averaging_intervals_borders[-1]:
                    interval_number =  len(averaging_intervals_borders)             
                
            #print(np.round(lon[profile],3),interval_number),np.round(np.nanmedian(np.log10(eps_grid[profile,30:-30])),3),np.shape(np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))))
                       

            #if the profile contains only nan values, profile is skipped
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,:])):
                print("NaN profile")
                continue
                
            #check for an outlier profile 
            if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+1.5*spread_of_profile_medians):      
                #print("\toutlier")
                count_outliers += 1
                continue
           
                
            #if the list is empty
            if np.any(oxygen_flux_statistic[interval_number]) == None:
                #fill it with a reshaped profile
                oxygen_flux_statistic[interval_number] = np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))
                salinity_statistic[interval_number] = np.reshape(eps_salinity_grid[profile,:],(1,-1))
                temperature_statistic[interval_number] = np.reshape(eps_consv_temperature_grid[profile,:],(1,-1)) 
                dissipation_statistic[interval_number] = np.reshape(eps_grid[profile,:],(1,-1)) 
                oxygen_sat_statistic[interval_number] = np.reshape(eps_oxygen_sat_grid[profile,:],(1,-1)) 
                    
            else:
                #concatenate all further profiles to the ones already in the array
                oxygen_flux_statistic[interval_number] = np.concatenate((oxygen_flux_statistic[interval_number],np.reshape(oxygen_flux_BB_grid[profile,:],(1,-1))),axis=0)
                salinity_statistic[interval_number] = np.concatenate((salinity_statistic[interval_number],np.reshape(eps_salinity_grid[profile,:],(1,-1))),axis=0)
                temperature_statistic[interval_number] = np.concatenate((temperature_statistic[interval_number],np.reshape(eps_consv_temperature_grid[profile,:],(1,-1))),axis=0)
                dissipation_statistic[interval_number] = np.concatenate((dissipation_statistic[interval_number],np.reshape(eps_grid[profile,:],(1,-1))),axis=0)
                oxygen_sat_statistic[interval_number] = np.concatenate((oxygen_sat_statistic[interval_number],np.reshape(eps_oxygen_sat_grid[profile,:],(1,-1))),axis=0)


        print("removed",count_outliers,"profiles as outliers")
        
        
        ###############################################################################################
        #Plotting of the outlier check pictures
        ###############################################################################################
        """
        cmap_hot = plt.get_cmap('hot_r')
        cmap_hot.set_bad(color = 'lightgrey')
        
        f1, axarr1 = plt.subplots(nrows = 2, ncols = 1, sharex = True)
        axarr1[0].plot(lon,np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        axarr1[0].axhline(transect_median,)
        axarr1[0].axhline(transect_median+1.5*spread_of_profile_medians, ls = "--")
        axarr1[0].axhline(transect_median-1.5*spread_of_profile_medians, ls = "--")
        axarr1[0].set_xlabel("longitude")        
        axarr1[0].set_ylabel(r"median(log10($\epsilon$)")  
               
        plotmesh_longitude = np.append(lon,2*lon[-1]-lon[-2])
        img1_1 = axarr1[1].pcolormesh(plotmesh_longitude,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -9.5, cmap = cmap_hot)
        
        #Draw interval borders
        #for border in averaging_intervals_borders:
        #    axarr1[1].axvline(border)
            
        axarr1[1].invert_yaxis()
        axarr1[1].set_ylabel("pressure [dbar]")   
        axarr1[1].set_xlabel("longitude") 
        colorbar(img1_1).set_label(r"log10(dissipation) $[m^2 s^{-3}]$")     
        f1.set_size_inches(9,5)
        f1.tight_layout()
         
        f1.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/outlier_visual_check/"+cruisename+"_"+transect_name+"_outlier", DPI = 300)         
        #plt.show()      
        plt.close(fig = "all")
        """
        
    ###########################################################################################################################
    
    for interval in oxygen_flux_statistic:
        print(np.shape(interval))
            
    #compute mean and std over the saved intervals
    mean_oxygen_flux = [None] * number_of_intervals
    std_oxygen_flux = [None] * number_of_intervals
    mean_salinity = [None] * number_of_intervals
    std_salinity = [None] * number_of_intervals
    mean_temperature = [None] * number_of_intervals
    std_temperature = [None] * number_of_intervals
    log_mean_dissipation = [None] * number_of_intervals
    log_std_dissipation = [None] * number_of_intervals
    
    mean_dissipation = [None] * number_of_intervals
    std_dissipation = [None] * number_of_intervals
    upper_percentile_dissipation = [None] * number_of_intervals
    lower_percentile_dissipation = [None] * number_of_intervals
    second_upper_percentile_dissipation = [None] * number_of_intervals
    second_lower_percentile_dissipation = [None] * number_of_intervals
            
    mean_oxygen_sat = [None] * number_of_intervals
    std_oxygen_sat = [None] * number_of_intervals
        
    for index in range(number_of_intervals):
        mean_oxygen_flux[index] = np.nanmean(oxygen_flux_statistic[index],axis=0)
        std_oxygen_flux[index] = np.nanstd(oxygen_flux_statistic[index],axis=0)
        mean_salinity[index] = np.nanmean(salinity_statistic[index],axis=0)
        std_salinity[index] = np.nanstd(salinity_statistic[index],axis=0)
        mean_temperature[index] = np.nanmean(temperature_statistic[index],axis=0)
        std_temperature[index] = np.nanstd(temperature_statistic[index],axis=0)
        log_mean_dissipation[index] = np.nanmean(np.log10(dissipation_statistic[index]),axis=0)
        log_std_dissipation[index] = np.nanstd(np.log10(dissipation_statistic[index]),axis=0) 
        
        mean_dissipation[index] = np.log10(np.nanmean(dissipation_statistic[index],axis=0))
        std_dissipation[index] = np.log10(np.nanstd(dissipation_statistic[index],axis=0)) 
        
        upper_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], 84.13, axis = 0))
        lower_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], 100-84.13, axis = 0))

        second_upper_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], 97.72, axis = 0))
        second_lower_percentile_dissipation[index] = np.log10(np.nanpercentile(dissipation_statistic[index], 100-97.72, axis = 0))
                
        mean_oxygen_sat[index] = np.nanmean(oxygen_sat_statistic[index],axis=0) 
        std_oxygen_sat[index] = np.nanstd(oxygen_sat_statistic[index],axis=0)
            
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
        
        
       
    f_flux,flux_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_salinity,sal_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_temperature,temp_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_dissipation,dissipation_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
    f_oxygen,oxygen_axarr = plt.subplots(nrows = 1, ncols = number_of_intervals, sharey = True, sharex = True) 
                    
    for index in range(number_of_intervals):
    
        if index == 0:
            title = "-"+str(np.round(averaging_intervals_borders[0],2))
        elif index == number_of_intervals-1:
            title = str(np.round(averaging_intervals_borders[-1],2))+"-"
        else:
            title = str(np.round(averaging_intervals_borders[index-1],2))+"-"+str(np.round(averaging_intervals_borders[index],2))   
            
        flux_axarr[index].plot(mean_oxygen_flux[index],eps_pressure)
        flux_axarr[index].fill_betweenx(eps_pressure,mean_oxygen_flux[index]-std_oxygen_flux[index],mean_oxygen_flux[index]+std_oxygen_flux[index], alpha = 0.5)
        flux_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        sal_axarr[index].plot(mean_salinity[index],eps_pressure)
        sal_axarr[index].fill_betweenx(eps_pressure,mean_salinity[index]-std_salinity[index],mean_salinity[index]+std_salinity[index], alpha = 0.5)
        sal_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        temp_axarr[index].plot(mean_temperature[index],eps_pressure)
        temp_axarr[index].fill_betweenx(eps_pressure,mean_temperature[index]-std_temperature[index],mean_temperature[index]+std_temperature[index], alpha = 0.5)
        temp_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        dissipation_axarr[index].plot(mean_dissipation[index],eps_pressure, label = "log of means")
        #dissipation_axarr[index].fill_betweenx(eps_pressure,mean_dissipation[index]-std_dissipation[index],mean_dissipation[index]+std_dissipation[index], alpha = 0.5)
        dissipation_axarr[index].fill_betweenx(eps_pressure,upper_percentile_dissipation[index],lower_percentile_dissipation[index], alpha = 0.6, label = "84.13% percentile")
        dissipation_axarr[index].fill_betweenx(eps_pressure,second_upper_percentile_dissipation[index],upper_percentile_dissipation[index], alpha = 0.4, color = "tab:blue",label = "97.72% percentile")
        dissipation_axarr[index].fill_betweenx(eps_pressure,second_lower_percentile_dissipation[index],lower_percentile_dissipation[index], alpha = 0.4, color = "tab:blue")
                        
        dissipation_axarr[index].plot(log_mean_dissipation[index],eps_pressure, label = "mean of logs")
        dissipation_axarr[index].fill_betweenx(eps_pressure,log_mean_dissipation[index]-log_std_dissipation[index],log_mean_dissipation[index]+log_std_dissipation[index], alpha = 0.4, label = "std")
        dissipation_axarr[index].legend()
        
        dissipation_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        oxygen_axarr[index].plot(mean_oxygen_sat[index],eps_pressure)
        oxygen_axarr[index].fill_betweenx(eps_pressure,mean_oxygen_sat[index]-std_oxygen_sat[index],mean_oxygen_sat[index]+std_oxygen_sat[index], alpha = 0.5)
        oxygen_axarr[index].set_title(str(np.shape(oxygen_flux_statistic[index])[0])+" profiles\n"+title)
        
        flux_axarr[index].set_xlabel(r"FO [mmol/(m$^2$*d]")
        oxygen_axarr[index].set_xlabel(r"O$_2$ saturation [$\%$]") 
        sal_axarr[index].set_xlabel("salinity [SA]") 
        temp_axarr[index].set_xlabel("consv_temperature [C]")
        dissipation_axarr[index].set_xlabel(r"log10($\epsilon$) $[m^2 s^{-3}]$")    
    
    flux_axarr[0].invert_yaxis()    
    flux_axarr[0].set_xlim((-10,10))
    flux_axarr[0].set_ylabel("pressure [dbar]")
    f_flux.suptitle("Oxygen flux")
       
    f_salinity.suptitle("Salinity")
    sal_axarr[0].invert_yaxis()
    sal_axarr[0].set_ylabel("pressure [dbar]")
        
    f_temperature.suptitle("temperature")
    temp_axarr[0].invert_yaxis()
    temp_axarr[0].set_ylabel("pressure [dbar]")
    
    f_dissipation.suptitle("dissipation")
    #dissipation_axarr[0].set_xlim((-10,-5))
    dissipation_axarr[0].invert_yaxis()
    dissipation_axarr[0].set_ylabel("pressure [dbar]")
    
    f_oxygen.suptitle("Oxygen saturation")
    oxygen_axarr[0].invert_yaxis()
    oxygen_axarr[0].set_ylabel("pressure [dbar]")    

    f_flux.set_size_inches(18,10.5)
    f_flux.tight_layout() 
    f_flux.subplots_adjust(top=0.91)
    f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_oxygen_flux")
    
    f_temperature.set_size_inches(18,10.5)
    f_temperature.tight_layout() 
    f_temperature.subplots_adjust(top=0.91)
    f_temperature.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_temperature")
    
    f_dissipation.set_size_inches(18,10.5)
    f_dissipation.tight_layout() 
    f_dissipation.subplots_adjust(top=0.91)
    f_dissipation.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_dissipation")
    
    f_oxygen.set_size_inches(18,10.5)
    f_oxygen.tight_layout() 
    f_oxygen.subplots_adjust(top=0.91)
    f_oxygen.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_oxygen")

    f_salinity.set_size_inches(18,10.5)
    f_salinity.tight_layout() 
    f_salinity.subplots_adjust(top=0.91)
    f_salinity.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_salinity")
    
    plt.show()
    
    
    
    
    
