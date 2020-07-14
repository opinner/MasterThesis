#---------------------------------------------------------#
#Plots one mss transect 
#plus as an example one profile from that transect
#---------------------------------------------------------#

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw 
import mss_functions as thesis






def load_clean_interpolate_and_shift_data(datafile_path, shift_value):
    """
    loads a a .mat file of mss measurements
    a few files are cleaned through a hardcoded routine (emb217/TR1-8.mat, emb169/TS118, emb169/TRR109)
    
    
    
    datafile_path: absolute path of the mss measurements saved as a .mat file


    return values:
        number_of_profiles              number of profiles/casts in the transect
        lat                             latitude in degrees (as a float) of the casts
        lon                             longitude in degrees (as a float) of the casts
        distance                        distance in km from the starting point of the transect
        
        interp_pressure                 equidistant pressure array between the highest and the lowest measured pressure value
        pressure_grid
        oxygen_grid                     
        oxygen_sat_grid                 oxygen saturation in percent as a grid (number_of_profiles x len(interp_pressure))
        salinity_grid                   salinity in g/kg as a grid (number_of_profiles x len(interp_pressure)) 
        consv_temperature_grid          conservative temperature in degrees Celsius as a grid (number_of_profiles x len(interp_pressure))
        density_grid                    density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
        
        eps_pressure                    1D array of pressure values to the dissipation rate values (the pressure distance between points is bigger than in interp_pressure) 
        eps_pressure_grid               
        eps_grid                        measured dissipation rate values (number_of_profiles x len(eps_pressure))
        eps_salinity_grid               salinity in g/kg as a grid (number_of_profiles x len(eps_pressure)) 
        eps_consv_temperature_grid      conservative temperature as a grid (number_of_profiles x len(eps_pressure))
        eps_oxygen_grid                 
        eps_oxygen_sat_grid             oxygen saturation as a grid (number_of_profiles x len(eps_pressure))
        eps_N_squared_grid              N^2, the Brunt-Vaisala frequency in 1/s^2 as a grid (number_of_profiles x len(eps_pressure))
        eps_density_grid                density in kg/m^3 as a grid (number_of_profiles x len(eps_pressure))
        
        0                               error code, returned if distance is not monotonically increasing (ergo a loop in the data)

    """
    import scipy.io as sio
    import geopy.distance as geo
    import numpy as np
    import gsw 
    

    
    splitted_path = datafile_path.split("/")
    cruisename = splitted_path[4][0:6]
    DATAFILENAME = splitted_path[-1]
    
    transect_name = DATAFILENAME[:-4]
     
    print(cruisename+"_"+transect_name)
    #print("Filename:",sio.whosmat(FILENAME))

    data = sio.loadmat(datafile_path) 


    STA_substructure = data["STA"]
    DATA_substructure = data["DATA"]
    MIX_substructure = data["MIX"]
    CTD_substructure = data["CTD"]

    #print(STA_substructure.dtype)
    #print(DATA_substructure.dtype)
    #print(CTD_substructure.dtype)
    #print(MIX_substructure.dtype)

    lat = STA_substructure["LAT"][0]
    lon = STA_substructure["LON"][0]

    pressure = CTD_substructure["P"][0]
    absolute_salinity = CTD_substructure["SA"][0] #is this unit sufficient
    consv_temperature = CTD_substructure["CT"][0] #TODO better use conservative temperature?
    alpha = CTD_substructure["ALPHA"][0]
    beta = CTD_substructure["BETA"][0]
    
    if cruisename == "emb217":
            oxygen_sat = CTD_substructure["O2"][0]
            
    elif cruisename == "emb177":
            try:
                oxygen_data = sio.loadmat("/home/ole/windows/all_data/emb177/deployments/ship/mss/mss38/TODL/TODL_merged_oxy/" +transect_name+"_TODL_merged_oxy.mat")
            except OSError:
                print("##########################################")
                print(cruisename,transect_name,"is skipped!")
                print("No oxygen data for this file")
                print("##########################################")
                return 0    
            
            #the following code block is written by Peter Holtermann and Copy pasted here

            emb177_oxygen = []
            emb177_oxygen_pressure = []
            for icast in range(len(lon)):
                oxy_p_temp   = oxygen_data['TODL_MSS']['P'][0,:][icast][0,:]
                #mT  = oxygen_data['TODL_MSS']['mT'][0,:][icast][0,:]
                #C   = oxygen_data['TODL_MSS']['C'][0,:][icast][0,:]
                oxy  = oxygen_data['TODL_MSS']['oxy'][0,:][icast][:,4]
                emb177_oxygen.append(oxy)
                emb177_oxygen_pressure.append(oxy_p_temp)
                
            
    elif cruisename == "emb169":
            try:
                oxygen_data = sio.loadmat("/home/ole/windows/all_data/emb169/deployments/ship/mss/TODL/merged/" +transect_name+"_TODL_merged_oxy.mat")
            except OSError:
                print("##########################################")
                print(cruisename,transect_name,"is skipped!")
                print("No oxygen data for this file")
                print("##########################################")
                return 0    
            
            #the following code block is written by Peter Holtermann and Copy pasted here

            emb169_oxygen = []
            emb169_oxygen_pressure = []
            for icast in range(len(lon)):
                oxy_p_temp   = oxygen_data['TODL_MSS']['P'][0,:][icast][0,:]
                #mT  = oxygen_data['TODL_MSS']['mT'][0,:][icast][0,:]
                #C   = oxygen_data['TODL_MSS']['C'][0,:][icast][0,:]
                oxy  = oxygen_data['TODL_MSS']['oxy'][0,:][icast][:,4]
                emb169_oxygen.append(oxy)
                emb169_oxygen_pressure.append(oxy_p_temp)
    else:
        raise AssertionError
        
    eps = MIX_substructure["eps"][0]
    eps_pressure = MIX_substructure["P"][0]
       

    number_of_profiles = np.shape(pressure)[-1]

    latitude = []
    longitude = []

    #calculates the distance in km of the position from every profile to the position of the starting profile.  
    distance = np.zeros(number_of_profiles)
    origin = (float(lat[0][0][0]),float(lon[0][0][0])) #lots of brackets to get a number, not an array (workaround)
    for i in range(number_of_profiles):
        current_point = (float(lat[i][0][0]),float(lon[i][0][0]))
        latitude.append(float(lat[i][0][0]))
        longitude.append(float(lon[i][0][0]))
        distance[i] = geo.geodesic(origin,current_point).km #Distance in km, change to nautical miles?

    lat = np.asarray(latitude)
    lon = np.asarray(longitude)


    #Data Cleaning
    #-----------------------------------------------------------------------------------------------------------------
   
    
    #remove data from a file, with overlapping positional points
    if cruisename == "emb217" and  transect_name == "TR1-8": 
        lat = np.delete(lat,np.s_[33:47])
        lon = np.delete(lon,np.s_[33:47])
        distance = np.delete(distance,np.s_[33:47])
        
        pressure = np.delete(pressure,np.s_[33:47],axis=0)
        oxygen_sat = np.delete(oxygen_sat,np.s_[33:47],axis=0)
        absolute_salinity =  np.delete(absolute_salinity,np.s_[33:47],axis=0)
        consv_temperature = np.delete(consv_temperature,np.s_[33:47],axis=0)
        alpha = np.delete(alpha,np.s_[33:47],axis=0)
        beta = np.delete(beta,np.s_[33:47],axis=0)
        
        
        eps = np.delete(eps,np.s_[33:47],axis=0)
        eps_pressure = np.delete(eps_pressure,np.s_[33:47],axis=0)
        
        number_of_profiles = np.shape(pressure)[-1]
   
   
    if cruisename == "emb169" and  DATAFILENAME[:-4] == "TS118":
        lat = lat[:21]
        lon = lon[:21]
        distance = distance[:21]
        
        pressure = pressure[:21]
        oxygen_sat = oxygen_sat[:21]
        absolute_salinity =  absolute_salinity[:21]
        consv_temperature = consv_temperature[:21]
        alpha = alpha[:21]
        beta = beta[:21]
        
        
        eps = eps[:21]
        eps_pressure = eps_pressure[:21]
        number_of_profiles = np.shape(pressure)[-1]
    
    #removes the last data point, that dont seem to belong to the transect
    if cruisename == "emb169" and  DATAFILENAME[:-4] == "TRR109":
        lat = lat[:-1]
        lon = lon[:-1]
        distance = distance[:-1]
        
        pressure = pressure[:-1]
        oxygen_sat = oxygen_sat[:-1]
        absolute_salinity =  absolute_salinity[:-1]
        consv_temperature = consv_temperature[:-1]
        alpha = alpha[:-1]
        beta = beta[:-1]
        
        
        eps = eps[:-1]
        eps_pressure = eps_pressure[:-1]
        number_of_profiles = np.shape(pressure)[-1]
    #-----------------------------------------------------------------------------------------------------------------    

    #test if distance is monotonically increasing
    try:
        assert(np.all(np.diff(distance)>0))
    except AssertionError:  
        print("##########################################")
        print(cruisename,DATAFILENAME[:-4],"is skipped!")
        print("Distance is not monotonically increasing. Mabye due to a loop in the transect?")
        print("##########################################")
        return 0
        #continue #jump to the next datafile


    #TODO point density? 

    #initial values 
    min_pressure = 10
    max_pressure = 60
    max_size = 1000
    min_size = 3000

    #select the min and max values for the equidistant pressure array (later used for the grid)
    for i in range(number_of_profiles):
        if np.nanmin(pressure[i]) < min_pressure:
            min_pressure = np.nanmin(pressure[i])
        if np.nanmax(pressure[i]) > max_pressure:
            max_pressure = np.nanmax(pressure[i])
        if pressure[i].size > max_size:
            max_size = pressure[i].size       
        if pressure[i].size < min_size:
            min_size = pressure[i].size  


    #print("High resolution",min_pressure,max_pressure,min_size,max_size)

    #check if that worked correctly
    assert(max_size>= min_size)


    #creates a pressure axis between the maximum and minimum pressure    
    interp_pressure = np.linspace(min_pressure,max_pressure,min_size)

    #creates pressure grid, where every column is equal to interp_pressure
    pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],min_size))

    #create grids that changes in distance on x and in depth on y-axis

    salinity_grid = np.zeros((np.shape(pressure)[-1],min_size))
    consv_temperature_grid = np.copy(salinity_grid)
    alpha_grid = np.copy(salinity_grid)
    beta_grid = np.copy(salinity_grid)

    if cruisename == "emb217":
        oxygen_sat_grid =  np.copy(salinity_grid)
    elif cruisename == "emb177" or cruisename == "emb169": 
        oxygen_grid =  np.copy(salinity_grid)
    
    #check if the pressure points for every eps profile are the same
    for i in range(number_of_profiles):  
        assert np.all(eps_pressure[i].flatten() == eps_pressure[0].flatten())
    #if yes the pressure can be coverted to a 1D array instead of a 2D array    
    eps_pressure = eps_pressure[0].flatten()
    
    #print(eps[i].flatten().size,eps_pressure.size, np.arange(1,160.5,0.5).size, np.arange(1,160.5,0.5).size - eps_pressure.size)
    #print(np.shape(eps),np.shape(eps[0]))
    
    
    desired_pressure_array = np.arange(1,160.5,0.5) #stems from the orginal matlab script that process the raw MSS data
    last_element = eps_pressure[-1]
    old_size = eps_pressure.size
    if last_element < 160:
        for index,element in enumerate(eps):
            eps[index] = np.pad(eps[index], ((0,desired_pressure_array.size - eps_pressure.size),(0,0)), 'constant', constant_values= np.nan)
        eps_pressure = np.append(eps_pressure,np.arange(last_element+0.5,160.5,0.5))
    assert(eps[0].flatten().size == eps_pressure.size)    
            
    #averaged of approx 5 depth bins (???)
    eps_grid = np.zeros((number_of_profiles,eps_pressure.size))

    #needed for the interpolation of S and T to the same grid as the eps
    eps_salinity_grid = np.ones((np.shape(eps_grid)))
    eps_consv_temperature_grid = np.ones(np.shape(eps_grid))

    if cruisename == "emb217":
        eps_oxygen_sat_grid = np.copy(eps_salinity_grid)
    elif cruisename == "emb177" or cruisename == "emb169": 
        eps_oxygen_grid = np.copy(eps_salinity_grid)

    #vector times matrix multiplication to get a 2D array, where every column is equal to eps_pressure
    eps_pressure_grid = np.reshape(eps_pressure,(1,-1))*np.ones(np.shape(eps_grid))

    #create a pressure axis where every point is shifted by half the distance to the next one
    shifted_pressure = eps_pressure + np.mean(np.diff(eps_pressure))/2

    #prepend a point at the beginning to be a point longer than the pressure axis we want from the Nsquared function
    shifted_pressure = np.append(eps_pressure[0]-np.mean(np.diff(eps_pressure))/2, shifted_pressure)

    #from that create a grid, where we have just n-times the pressure axis with n the number of profiles 
    shifted_pressure_grid = np.reshape(shifted_pressure,(1,-1))*np.ones((number_of_profiles,shifted_pressure.size))
    shifted_salinity_grid = np.ones((np.shape(shifted_pressure_grid)))
    shifted_consv_temperature_grid = np.ones((np.shape(shifted_pressure_grid)))



    #Grid interpolation (loops over all profiles)
    for i in range(number_of_profiles):      
             
        #interpolation to a common fine grid    
        salinity_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
        consv_temperature_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
        alpha_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),alpha[i].flatten(), left = np.nan, right = np.nan)
        beta_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),beta[i].flatten(), left = np.nan, right = np.nan)
            
        #interpolation to a shifted grid for the Nsquared function
        shifted_salinity_grid[i] = np.interp(shifted_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
        shifted_consv_temperature_grid[i] = np.interp(shifted_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
        
        #just changes the format of eps slightly
        assert(eps[i].flatten().size == eps_pressure.size)
        eps_grid[i] = eps[i].flatten()

        #interpolate S and T to the same grid as eps
        eps_salinity_grid[i] = np.interp(eps_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
        eps_consv_temperature_grid[i] = np.interp(eps_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
        
        #interpolation of the oxygen depends on which source it is
        if cruisename == "emb217":
            oxygen_sat_grid[i] = np.interp(interp_pressure,pressure[i].flatten()+shift_value,oxygen_sat[i].flatten(), left = np.nan, right = np.nan)
            eps_oxygen_sat_grid[i] = np.interp(eps_pressure,pressure[i].flatten()+shift_value,oxygen_sat[i].flatten(), left = np.nan, right = np.nan)       
        elif cruisename == "emb177":
            oxygen_grid[i] = np.interp(interp_pressure,emb177_oxygen_pressure[i]+shift_value,emb177_oxygen[i],left=np.nan,right=np.nan)
            eps_oxygen_grid[i] = np.interp(eps_pressure,emb177_oxygen_pressure[i]+shift_value,emb177_oxygen[i].flatten(), left = np.nan, right = np.nan)  
        elif cruisename == "emb169": 
            oxygen_grid[i] = np.interp(interp_pressure,emb169_oxygen_pressure[i]+shift_value,emb169_oxygen[i],left=np.nan,right=np.nan)
            eps_oxygen_grid[i] = np.interp(eps_pressure,emb169_oxygen_pressure[i]+shift_value,emb169_oxygen[i].flatten(), left = np.nan, right = np.nan)  
            
        
    if cruisename == "emb217":
        assert(np.shape(oxygen_sat_grid) == np.shape(salinity_grid))
        assert(np.shape(eps_oxygen_sat_grid) == np.shape(eps_salinity_grid))
        
    if cruisename == "emb177" or cruisename == "emb169":
        assert(np.shape(oxygen_grid) == np.shape(salinity_grid))
        assert(np.shape(eps_oxygen_grid) == np.shape(eps_salinity_grid))
                    
    #fine density grid
    density_grid = gsw.rho(salinity_grid,consv_temperature_grid,pressure_grid)
    pot_density_grid = 1000+gsw.density.sigma0(salinity_grid,consv_temperature_grid)

    #density grid on the same points as the eps grid
    eps_density_grid = gsw.rho(eps_salinity_grid,eps_consv_temperature_grid,eps_pressure_grid)
    eps_pot_density_grid = 1000+gsw.density.sigma0(eps_salinity_grid,eps_consv_temperature_grid)
    
    #TODO compare with density_grid?
    #density_grid_check = (1 - alpha_grid * (consv_temperature_grid) + beta_grid * (salinity_grid))*rho_0
    #difference = density_grid-density_grid_check

    #calculate N^2 with the gsw toolbox, by using the shifted grid we should get a N^2 grid that is defined at the same points as eps_grid
    eps_N_squared_grid, crosscheck_pressure_grid = gsw.Nsquared(shifted_salinity_grid,shifted_consv_temperature_grid,shifted_pressure_grid, lat = np.mean(lat), axis = 1)
    crosscheck_pressure = np.mean(crosscheck_pressure_grid, axis = 0)

    #test if we really calculated N^2 for the same points as pressure points from the dissipation measurement
    assert(np.all(crosscheck_pressure == eps_pressure))
    
    #create grids that have the latitude/longitude values for every depth (size: number_of_profiles x len(interp_pressure))
    lat_grid = np.reshape(lat,(-1,1)) * np.ones((number_of_profiles,max_size))
    lon_grid = np.reshape(lon,(-1,1)) * np.ones((number_of_profiles,max_size))
    
    #check if lon grid was correctly created
    assert(np.all(lon_grid[:,0] == lon))
    
    #create grids that have the latitude/longitude values for every depth (size: number_of_profiles x len(eps_pressure))
    eps_lat_grid = np.reshape(lat,(-1,1)) * np.ones((number_of_profiles,eps_pressure.size))
    eps_lon_grid = np.reshape(lat,(-1,1)) * np.ones((number_of_profiles,eps_pressure.size))
    
    if cruisename == "emb217":
        #convert oxygen saturation to oxygen concentration (functions are selfwritten but use the gsw toolbox, for more informations see the function (also in this file))
        oxygen_grid = thesis.oxygen_saturation_to_concentration(oxygen_sat_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
        eps_oxygen_grid = thesis.oxygen_saturation_to_concentration(eps_oxygen_sat_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid) 
              
    elif cruisename == "emb177":
        #scale the oxygen to be at 100% at the surface
        maximum_concentration_grid = gsw.O2sol(salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
        
        correction_factor = np.ones((number_of_profiles,1))
        
        for profile in range(number_of_profiles):
            for i in range(100):   
                if np.isnan(oxygen_grid[profile,i]) or np.isnan(maximum_concentration_grid[profile,i]):
                    continue
                else:
                    correction_factor[profile] = maximum_concentration_grid[profile,i]/oxygen_grid[profile,i]
                    break
                    
        oxygen_grid = oxygen_grid * correction_factor
        eps_oxygen_grid = eps_oxygen_grid * correction_factor
        
        #convert oxygen concentration to oxygen saturation (functions are selfwritten but use the gsw toolbox, for more informations see the function (also in this file))
        oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(oxygen_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
        eps_oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(eps_oxygen_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid)  
    
    elif cruisename == "emb169":
        #convert oxygen concentration to oxygen saturation (functions are selfwritten but use the gsw toolbox, for more informations see the function (also in this file))
        oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(oxygen_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
        eps_oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(eps_oxygen_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid)  
    
    else:
        raise AssertionError
        
        
          
    return [[number_of_profiles,lat,lon,distance],[interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid,pot_density_grid],[eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid,eps_pot_density_grid]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#########################################################################################################################################
            


datafile_path = "/home/ole/windows/emb217_mss_data/TR1-8.mat"

splitted_filename = datafile_path.split("/")
cruisename = splitted_filename[4][0:6]
DATAFILENAME = splitted_filename[-1]
print("cruisename",cruisename)  



    
print("Filename:",sio.whosmat(datafile_path))
shift_values = [-0.5,-0.25,0.25,0.5,0]


f,axarr = plt.subplots(1,2, sharey = True)
f2,axarr2 = plt.subplots(2,1, sharex = True)

for shift_index,shift_value in enumerate(shift_values):


    results = load_clean_interpolate_and_shift_data(datafile_path, shift_value)
    
    try:
        number_of_profiles,lat,lon,distance = results[0]
        interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid, pot_density_grid = results[1]
        eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid, eps_pot_density_grid = results[2]
    except TypeError:
        #print(cruisename,transect_name,"is skipped!")
        continue


    #calculate the viscosity (2 different formula)
    eps_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_salinity_grid,eps_density_grid, "default")


    #calculate the Reynolds bouyancy number defined on the pressure values defined in eps_pressure
    eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared_grid)

    #A negative Reynolds bouyancy has no physical meaning, therefore gets replaced by NaNs
    #These negative values occur through a negative squared bouyancy frequency
    eps_Reynolds_bouyancy_grid[eps_Reynolds_bouyancy_grid < 0] = np.nan         
                
    #apply the correction on the dissipation rate following Garanaik2018           
    #corrected_eps_grid = thesis.correct_dissipation(eps_grid,eps_Reynolds_bouyancy_grid)

    #from that calculate again the now corrected Bouyancy Reynolds number (in 3D)
    #corrected_eps_Reynolds_bouyancy_grid = corrected_eps_grid/(eps_viscosity_grid*eps_N_squared_grid)
    
    #conversion from pressure coordinates to depth
    eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west   
    eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid)) #shape to the same from as eps_grid
        
    oxygen_flux_Osborn_grid = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
    oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
    oxygen_flux_Skif_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        
    
    profile_index = np.argmin(np.abs(lon-20.6)) 
    a = np.argmin(np.abs(eps_pressure-65))   
    b = np.argmin(np.abs(eps_pressure-75)) 

    #axarr[0].plot(np.nanmean(eps_oxygen_sat_grid, axis = 0),eps_pressure,label = str(shift_value))
    #axarr[1].plot(np.nanmean(oxygen_flux_Skif_grid, axis = 0),eps_pressure,label = str(shift_value))
    axarr[0].plot(eps_oxygen_sat_grid[profile_index],eps_pressure,label = str(shift_value))
    axarr[1].plot(oxygen_flux_Skif_grid[profile_index],eps_pressure,label = str(shift_value))
    
    axarr2[0].plot(lon,np.nanmax(eps_oxygen_sat_grid[:,a:b], axis = 1),label = str(shift_value))
    #axarr2[1].plot(lon,np.nanmax(oxygen_flux_Skif_grid[:,a:b], axis = 1),label = str(shift_value))
    axarr2[0].plot(lon,np.nanmin(eps_oxygen_sat_grid[:,a:b], axis = 1),"--",label = str(shift_value))
    axarr2[1].plot(lon,np.nanmean(oxygen_flux_Skif_grid[:,a:b], axis = 1),"--",label = str(shift_value))    


axarr[0].invert_yaxis()
axarr[0].legend()
axarr[1].legend()

axarr2[0].legend()
axarr2[1].legend()

plt.show()



