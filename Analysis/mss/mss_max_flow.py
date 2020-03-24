############################################################
#this program loads all profiles from a cruise, removes outliers
#and retrieves the maximum oxyen flux in the lowermost meters 
#of the water column in choosable longitude intervals
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

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)
    
    
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169","/home/ole/share-windows/processed_mss/emb177"]

#averaging_intervals_borders = [20.55,20.62]
averaging_intervals_borders = np.linspace(20.48,20.7,9)
number_of_intervals = len(averaging_intervals_borders)+1

for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print(cruisename)    
    #2 borders = 3 intervals
    oxygen_flux_BB_statistic = [None] * number_of_intervals
    oxygen_flux_Shih_statistic = [None] * number_of_intervals
    oxygen_flux_Osborn_statistic = [None] * number_of_intervals
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
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02)
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
        print(averaging_intervals_borders)
        outlier_count = 0
        transect_oxygen_flux_statistic = []
        
        for profile in range(number_of_profiles):
            for index,border in enumerate(averaging_intervals_borders):
                #print("index",index)
                
                #Sort the prifle into tthe intervals
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
            if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+2*spread_of_profile_medians):      
                #print("\toutlier")
                outlier_count += 1
                continue
           

            from_index = int(list_of_BBL_range_indices[profile])     
            to_index = int(list_of_bathymetrie_indices[profile])
            #print(from_index,to_index)

            min_max_array = np.asarray([[np.nanmin(oxygen_flux_BB_grid[profile,from_index:to_index]),np.nanmax(oxygen_flux_BB_grid[profile,from_index:to_index])]])
            print(min_max_array)

            #if the list is empty
            if np.any(oxygen_flux_statistic[interval_number]) == None:
  
                #fill it with a reshaped profile
                oxygen_flux_statistic[interval_number] = min_max_array
           
            else:
                #concatenate all further profiles to the ones already in the array
                oxygen_flux_statistic[interval_number] = np.concatenate((oxygen_flux_statistic[interval_number],min_max_array),axis=0)
                
        transect_oxygen_flux_statistic = np.append(transect_oxygen_flux_statistic,min_max_array[0])

        print("removed",outlier_count,"profiles as outliers")
        for interval in oxygen_flux_statistic:
            print(np.shape(interval))
        print(np.shape(transect_oxygen_flux_statistic))
        print(transect_oxygen_flux_statistic)
        ###############################################################################################
        #Plotting of the maximum flux values per transect
        ###############################################################################################
        f1, axarr1 = plt.subplots(nrows = 1, ncols = 1, sharex = True)
        axarr1.plot(lon,transect_oxygen_flux_statistic[:,1])
        axarr1.set_xlabel("longitude")        
                   
        f1.set_size_inches(9,5)
        f1.tight_layout()
         
        f1.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/max_flux_transect/"+cruisename+"_"+transect_name+"_max_flux", DPI = 300)         
        plt.show()      
        plt.close(fig = "all")
        
        
    ###########################################################################################################################
    #compute mean and std over the saved intervals
    mean_max_flux = [None] * number_of_intervals
    std_max_flux = [None] * number_of_intervals

        
    for index in range(number_of_intervals):
        mean_oxygen_flux[index] = np.nanmean(oxygen_flux_statistic[index],axis=0)
        std_oxygen_flux[index] = np.nanstd(oxygen_flux_statistic[index],axis=0)

            
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    #####################################################PLOTTING#####################################################################
    ##################################################################################################################################
    ##################################################################################################################################
    ##################################################################################################################################     
        
        
       
    f_flux,flux_axarr = plt.subplots(nrows = 1, ncols = 1, sharey = True, sharex = True) 
    
    #append one extra border behind the last border in the mean distance of the borders 
    plot_longitude = np.append(averaging_intervals_borders,averaging_intervals_borders[-1]+np.mean(np.diff(averaging_intervals_borders)))
    #shift all plot points by half the border distance  
    plot_longitude - np.mean(np.diff(averaging_intervals_borders))/2
            
    flux_axarr.plot(plot_longitude,mean_max_flux)
    flux_axarr.fill_betweenx(plot_longitude,mean_max_flux-std_max_flux,mean_max_flux+std_max_flux, alpha = 0.5)

    flux_axarr[0].invert_yaxis()    
    flux_axarr[0].set_xlim((-20,20))
    flux_axarr[0].set_ylabel("pressure [dbar]")
    f_flux.suptitle("Oxygen flux")

    f_flux.set_size_inches(18,10.5)
    f_flux.tight_layout() 
    f_flux.subplots_adjust(top=0.91)
    f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruisename+"_"+str(number_of_intervals)+"intervals_oxygen_flux")
    

    
    plt.show()
    
    
    
    
    
