#---------------------------------------------------------#

#---------------------------------------------------------#

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw 
import mss_functions as thesis
import warnings
warnings.filterwarnings('ignore')



def load_data(datafile_path):
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
    
           
    if cruisename == "emb177":
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
            
    return [lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, emb177_oxygen,emb177_oxygen_pressure, number_of_profiles]
    
    
    
    
def interpolation1(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, emb177_oxygen, emb177_oxygen_pressure, number_of_profiles):
    """
    02 to eps, then d02/dz
    """

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

    oxygen_grid =  np.copy(salinity_grid)
    
    


    #averaged of approx 5 depth bins (???)
    eps_grid = np.zeros((number_of_profiles,eps_pressure.size))

    #needed for the interpolation of S and T to the same grid as the eps
    eps_salinity_grid = np.ones((np.shape(eps_grid)))
    eps_consv_temperature_grid = np.ones(np.shape(eps_grid))



    #vector times matrix multiplication to get a 2D array, where every column is equal to eps_pressure
    eps_pressure_grid = np.reshape(eps_pressure,(1,-1))*np.ones(np.shape(eps_grid))

    eps_oxygen_grid = np.copy(eps_salinity_grid)



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
        
        oxygen_grid[i] = np.interp(interp_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i],left=np.nan,right=np.nan)
        eps_oxygen_grid[i] = np.interp(eps_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i].flatten(), left = np.nan, right = np.nan)  

            
        
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

    
    
    
    return [[number_of_profiles,lat,lon,distance],[interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid,pot_density_grid],[eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid,eps_pot_density_grid]]
    
    
def interpolation2(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, emb177_oxygen, emb177_oxygen_pressure, number_of_profiles):
    """
    02, d02/dz then to eps
    """
    
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

    #oxygen_gradient_grid =  np.copy(salinity_grid)
  
    


    #averaged of approx 5 depth bins (???)
    eps_grid = np.zeros((number_of_profiles,eps_pressure.size))

    #needed for the interpolation of S and T to the same grid as the eps
    eps_salinity_grid = np.ones((np.shape(eps_grid)))
    eps_consv_temperature_grid = np.ones(np.shape(eps_grid))



    #vector times matrix multiplication to get a 2D array, where every column is equal to eps_pressure
    eps_pressure_grid = np.reshape(eps_pressure,(1,-1))*np.ones(np.shape(eps_grid))

    eps_oxygen_gradient_grid = np.copy(eps_salinity_grid)
    eps_oxygen_grid = np.copy(eps_salinity_grid)


    #create a pressure axis where every point is shifted by half the distance to the next one
    shifted_pressure = eps_pressure + np.mean(np.diff(eps_pressure))/2

    #prepend a point at the beginning to be a point longer than the pressure axis we want from the Nsquared function
    shifted_pressure = np.append(eps_pressure[0]-np.mean(np.diff(eps_pressure))/2, shifted_pressure)

    #from that create a grid, where we have just n-times the pressure axis with n the number of profiles 
    shifted_pressure_grid = np.reshape(shifted_pressure,(1,-1))*np.ones((number_of_profiles,shifted_pressure.size))
    shifted_salinity_grid = np.ones((np.shape(shifted_pressure_grid)))
    shifted_consv_temperature_grid = np.ones((np.shape(shifted_pressure_grid)))

    #print("shape emb177_oxygen_pressure",np.shape(emb177_oxygen_pressure),np.shape(emb177_oxygen_pressure[34]))
    # #mean lat should be sufficient, because the transect is east-west   
    #oxygen_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid)) #shape to the same from as eps_grid


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
        
        
     
    #create grids that have the latitude/longitude values for every depth (size: number_of_profiles x len(interp_pressure))
    lat_grid = np.reshape(lat,(-1,1)) * np.ones((number_of_profiles,max_size))
    lon_grid = np.reshape(lon,(-1,1)) * np.ones((number_of_profiles,max_size))
    maximum_concentration_grid = gsw.O2sol(salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
              
    #Grid interpolation (loops over all profiles)
    for i in range(number_of_profiles):        

        #Assumes a well mixed upper layer to get a constant correction factor (this is a little bit cheating)
        correction_factor = np.nanmean(maximum_concentration_grid[i,0:20])/np.nanmean(emb177_oxygen[i][0:20])
        emb177_oxygen[i] = emb177_oxygen[i] * correction_factor
        oxygen_depth = gsw.z_from_p(emb177_oxygen_pressure[i],np.mean(lat))
        emb177_oxygen_gradient =  thesis.central_differences(emb177_oxygen[i])/thesis.central_differences(oxygen_depth)
        eps_oxygen_gradient_grid[i] = np.interp(eps_pressure,emb177_oxygen_pressure[i],emb177_oxygen_gradient.flatten(), left = np.nan, right = np.nan)  
        eps_oxygen_grid[i] = np.interp(eps_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i],left=np.nan,right=np.nan)
            
        
    if cruisename == "emb177" or cruisename == "emb169":
        assert(np.shape(eps_oxygen_grid) == np.shape(eps_salinity_grid))
        assert(np.shape(eps_oxygen_gradient_grid) == np.shape(eps_salinity_grid))
                    
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
    

    
    #convert oxygen concentration to oxygen saturation (functions are selfwritten but use the gsw toolbox, for more informations see the function (also in this file))
    #oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(oxygen_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
    #eps_oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(eps_oxygen_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid)  
    

    
    #conversion from pressure coordinates to depth

    
    
    
    
    return [[number_of_profiles,lat,lon,distance],[interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid,pot_density_grid],[eps_pressure,eps_oxygen_grid,eps_oxygen_gradient_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid,eps_pot_density_grid]]    

    
    
    
        
    
    
def interpolation3(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, emb177_oxygen, emb177_oxygen_pressure, number_of_profiles):
    """
    eps to 02 and d02/dz
    """
    
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

    oxygen_grid =  np.copy(salinity_grid)
    
    


    #averaged of approx 5 depth bins (???)
    eps_grid = np.zeros((number_of_profiles,interp_pressure.size))

    #needed for the interpolation of S and T to the same grid as the eps
    #eps_salinity_grid = np.ones((np.shape(eps_grid)))
    #eps_consv_temperature_grid = np.ones(np.shape(eps_grid))



    #vector times matrix multiplication to get a 2D array, where every column is equal to eps_pressure
    #eps_pressure_grid = np.reshape(eps_pressure,(1,-1))*np.ones(np.shape(eps_grid))

    #eps_oxygen_grid = np.copy(eps_salinity_grid)



    #create a pressure axis where every point is shifted by half the distance to the next one
    shifted_pressure = interp_pressure + np.mean(np.diff(interp_pressure))/2

    #prepend a point at the beginning to be a point longer than the pressure axis we want from the Nsquared function
    shifted_pressure = np.append(interp_pressure[0]-np.mean(np.diff(interp_pressure))/2, shifted_pressure)

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
        eps_grid[i] = np.interp(interp_pressure,eps_pressure,eps[i].flatten(), left = np.nan, right = np.nan)

        #interpolate S and T to the same grid as eps
        #eps_salinity_grid[i] = np.interp(eps_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
        #eps_consv_temperature_grid[i] = np.interp(eps_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
        
        oxygen_grid[i] = np.interp(interp_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i],left=np.nan,right=np.nan)
        #eps_oxygen_grid[i] = np.interp(eps_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i].flatten(), left = np.nan, right = np.nan)  

            
        
    if cruisename == "emb217":
        assert(np.shape(oxygen_sat_grid) == np.shape(salinity_grid))
        #assert(np.shape(eps_oxygen_sat_grid) == np.shape(eps_salinity_grid))
        
    if cruisename == "emb177" or cruisename == "emb169":
        assert(np.shape(oxygen_grid) == np.shape(salinity_grid))
        #assert(np.shape(eps_oxygen_grid) == np.shape(eps_salinity_grid))
                    
    #fine density grid
    density_grid = gsw.rho(salinity_grid,consv_temperature_grid,pressure_grid)
    pot_density_grid = 1000+gsw.density.sigma0(salinity_grid,consv_temperature_grid)

    #density grid on the same points as the eps grid
    #eps_density_grid = gsw.rho(eps_salinity_grid,eps_consv_temperature_grid,eps_pressure_grid)
    #eps_pot_density_grid = 1000+gsw.density.sigma0(eps_salinity_grid,eps_consv_temperature_grid)
    
    #TODO compare with density_grid?
    #density_grid_check = (1 - alpha_grid * (consv_temperature_grid) + beta_grid * (salinity_grid))*rho_0
    #difference = density_grid-density_grid_check

    #calculate N^2 with the gsw toolbox, by using the shifted grid we should get a N^2 grid that is defined at the same points as eps_grid
    N_squared_grid, crosscheck_pressure_grid = gsw.Nsquared(shifted_salinity_grid,shifted_consv_temperature_grid,shifted_pressure_grid, lat = np.mean(lat), axis = 1)
    crosscheck_pressure = np.mean(crosscheck_pressure_grid, axis = 0)

    #for a,b in zip(crosscheck_pressure,interp_pressure):
    #    print(a,b)
    #test if we really calculated N^2 for the same points as pressure points from the dissipation measurement
    #assert(np.all(crosscheck_pressure == interp_pressure))
    np.testing.assert_allclose(crosscheck_pressure, interp_pressure, rtol=1e-5, atol=0)
    
    #create grids that have the latitude/longitude values for every depth (size: number_of_profiles x len(interp_pressure))
    lat_grid = np.reshape(lat,(-1,1)) * np.ones((number_of_profiles,max_size))
    lon_grid = np.reshape(lon,(-1,1)) * np.ones((number_of_profiles,max_size))
    
    #check if lon grid was correctly created
    assert(np.all(lon_grid[:,0] == lon))
    
    #create grids that have the latitude/longitude values for every depth (size: number_of_profiles x len(eps_pressure))
    #eps_lat_grid = np.reshape(lat,(-1,1)) * np.ones((number_of_profiles,interp_pressure.size))
    #eps_lon_grid = np.reshape(lat,(-1,1)) * np.ones((number_of_profiles,interp_pressure.size))
    
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
    #eps_oxygen_grid = eps_oxygen_grid * correction_factor
    
    #convert oxygen concentration to oxygen saturation (functions are selfwritten but use the gsw toolbox, for more informations see the function (also in this file))
    #oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(oxygen_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
    #eps_oxygen_sat_grid = thesis.oxygen_concentration_to_saturation(eps_oxygen_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid)  

    
    
    
    return [[number_of_profiles,lat,lon,distance],[interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid,pot_density_grid,eps_grid,N_squared_grid]]
    
    
        
    
def interpolation4(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, emb177_oxygen, emb177_oxygen_pressure, number_of_profiles):
    """
    average in bins to eps
    """

    import scipy.stats as scs

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

    oxygen_grid =  np.copy(salinity_grid)
    
    


    #averaged of approx 5 depth bins (???)
    eps_grid = np.zeros((number_of_profiles,eps_pressure.size))

    #needed for the interpolation of S and T to the same grid as the eps
    eps_salinity_grid = np.ones((np.shape(eps_grid)))
    eps_consv_temperature_grid = np.ones(np.shape(eps_grid))



    #vector times matrix multiplication to get a 2D array, where every column is equal to eps_pressure
    eps_pressure_grid = np.reshape(eps_pressure,(1,-1))*np.ones(np.shape(eps_grid))

    eps_oxygen_grid = np.copy(eps_salinity_grid)



    #create a pressure axis where every point is shifted by half the distance to the next one
    shifted_eps_pressure = eps_pressure + np.mean(np.diff(eps_pressure))/2

    #prepend a point at the beginning to be a point longer than the pressure axis we want from the Nsquared function
    shifted_eps_pressure = np.append(eps_pressure[0]-np.mean(np.diff(eps_pressure))/2, shifted_eps_pressure)
    print("eps_pressure",len(eps_pressure),eps_pressure[0],eps_pressure[-1])
    print("shifted_eps_pressure",len(shifted_eps_pressure),shifted_eps_pressure[0],shifted_eps_pressure[-1])
    print("appended eps pressure", len(np.append(shifted_eps_pressure,160.5))) 
    #from that create a grid, where we have just n-times the pressure axis with n the number of profiles 
    shifted_eps_pressure_grid = np.reshape(shifted_eps_pressure,(1,-1))*np.ones((number_of_profiles,shifted_eps_pressure.size))
    shifted_salinity_grid = np.ones((np.shape(shifted_eps_pressure_grid)))
    shifted_consv_temperature_grid = np.ones((np.shape(shifted_eps_pressure_grid)))



    #Grid interpolation (loops over all profiles)
    for i in range(number_of_profiles):      
             
        #interpolation to a common fine grid    
        salinity_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
        consv_temperature_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
        alpha_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),alpha[i].flatten(), left = np.nan, right = np.nan)
        beta_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),beta[i].flatten(), left = np.nan, right = np.nan)
         
        
        #bin edges of eps_pressure for the Nsquared function, so that the bins are centered around shifted_eps_pressure
        shift_sal_pressure = pressure[i].flatten()[~np.isnan(absolute_salinity[i].flatten())]
        shift_sal = absolute_salinity[i].flatten()[~np.isnan(absolute_salinity[i].flatten())]
        #print(np.shape(shift_sal_pressure),np.shape(shift_sal),np.shape(eps_pressure))
        shifted_salinity_grid[i] = scs.binned_statistic(x = shift_sal_pressure, values = shift_sal, statistic = "mean", bins = np.append([0.5],np.append(eps_pressure,160.5)))[0]
        
        shifted_consv_temperature_grid[i] = scs.binned_statistic(pressure[i].flatten()[~np.isnan(consv_temperature[i].flatten())],consv_temperature[i].flatten()[~np.isnan(consv_temperature[i].flatten())], statistic = "mean", bins = np.append([0.5],np.append(eps_pressure,160.5)))[0]
        
        """
        shifted_bin_salinity_grid[i] = scs.binned_statistic(x = shift_sal_pressure, values = shift_sal, statistic = "mean", bins = np.append([0.5],np.append(eps_pressure,160.5)))[0]
        shifted_bin_consv_temperature_grid[i] = scs.binned_statistic(pressure[i].flatten()[~np.isnan(consv_temperature[i].flatten())],consv_temperature[i].flatten()[~np.isnan(consv_temperature[i].flatten())], statistic = "mean", bins = np.append([0.5],np.append(eps_pressure,160.5)))[0]
        """
        
        #just changes the format of eps slightly
        assert(eps[i].flatten().size == eps_pressure.size)
        eps_grid[i] = eps[i].flatten()

        #bin S and T with bin_edges of shifted_eps_pressure, so that the bins are centered around eps_pressure (used for calculating the density)
        eps_salinity_grid[i] = scs.binned_statistic(pressure[i].flatten()[~np.isnan(absolute_salinity[i].flatten())],absolute_salinity[i].flatten()[~np.isnan(absolute_salinity[i].flatten())], statistic = "mean", bins = shifted_eps_pressure)[0]
        eps_consv_temperature_grid[i] = scs.binned_statistic(pressure[i].flatten()[~np.isnan(consv_temperature[i].flatten())],consv_temperature[i].flatten()[~np.isnan(consv_temperature[i].flatten())], statistic = "mean", bins = shifted_eps_pressure)[0]
        
        oxygen_grid[i] = np.interp(interp_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i],left=np.nan,right=np.nan)
        
        oxygen_wo_nan_pressure = emb177_oxygen_pressure[i][~np.isnan(emb177_oxygen[i].flatten())]
        oxygen_wo_nan = emb177_oxygen[i].flatten()[~np.isnan(emb177_oxygen[i].flatten())]
        #print(len(oxygen_wo_nan_pressure),len(oxygen_wo_nan))
        if len(oxygen_wo_nan) == 0:
            eps_oxygen_grid[i] = np.ones(len(eps_oxygen_grid[i]))*np.nan
        else:
            eps_oxygen_grid[i] = scs.binned_statistic(oxygen_wo_nan_pressure, oxygen_wo_nan, statistic = "mean", bins = shifted_eps_pressure)[0]

            
        
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
    eps_N_squared_grid, crosscheck_pressure_grid = gsw.Nsquared(shifted_salinity_grid,shifted_consv_temperature_grid,shifted_eps_pressure_grid, lat = np.mean(lat), axis = 1)
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

    
    
    
    return [[number_of_profiles,lat,lon,distance],[interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid,pot_density_grid],[eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid,eps_pot_density_grid]]
    
        
#########################################################################################################################################    
    

    
    
#########################################################################################################################################
            
f,axis = plt.subplots(ncols = 3, sharey = True)

colors = ["tab:blue","tab:green","tab:red","k"]

#datafile_path = "/home/ole/windows/emb217_mss_data/TR1-8.mat"
datafile_path = "/home/ole/windows/emb177_mss_data/TS1_8.mat"


splitted_filename = datafile_path.split("/")
cruisename = splitted_filename[4][0:6]
DATAFILENAME = splitted_filename[-1]
print("cruisename",cruisename)  

assert cruisename == "emb177"

upper_bound_halocline_as_density = 1006.9 #1006.9
lower_bound_halocline_as_density = 1008.2 #1007.9   
    
#print("Filename:",sio.whosmat(datafile_path))



data = load_data(datafile_path)
lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, oxygen, oxygen_pressure, number_of_profiles = data

profile_index = np.argmin(np.abs(lon-20.643))
#profile_index = np.argmin(np.abs(lon-20.546))
        

print(cruisename,DATAFILENAME[:-4],lon[profile_index])

#-------------------------------------------------------------------
#Interpolation 1 O2 to eps grid, then calculation of d02/dz
#-------------------------------------------------------------------
results = interpolation1(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, oxygen, oxygen_pressure, number_of_profiles)
    

number_of_profiles,lat,lon,distance = results[0]
interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid, pot_density_grid = results[1]
eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid, eps_pot_density_grid = results[2]

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
    



axis[0].plot(np.log10(eps_grid[profile_index,:]),eps_pressure, colors[0], label = "O2 to eps, then do2/dz")
axis[1].plot(thesis.central_differences(eps_oxygen_grid[profile_index,:])/thesis.central_differences(eps_depth),eps_pressure, colors[0], label = "O2 to eps, then do2/dz1")
axis[1].plot(eps_oxygen_grid[profile_index,:],eps_pressure, "--", c = colors[0])#, label = "O2 to eps, then do2/dz")
axis[2].plot(oxygen_flux_Skif_grid[profile_index,:],eps_pressure, colors[0], label = "O2 to eps, then do2/dz")


#-------------------------------------------------------------------
#Interpolation 2: O2, calculation of d02/dz, then to eps grid
#-------------------------------------------------------------------
results = interpolation2(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, oxygen, oxygen_pressure, number_of_profiles)
    
number_of_profiles,lat,lon,distance = results[0]
#interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid, pot_density_grid = results[1]
eps_pressure,eps_oxygen_grid,eps_oxygen_gradient_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid, eps_pot_density_grid = results[2]


#calculate the viscosity (2 different formula)
eps_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_salinity_grid,eps_density_grid, "default")


#calculate the Reynolds bouyancy number defined on the pressure values defined in eps_pressure
eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared_grid)

#A negative Reynolds bouyancy has no physical meaning, therefore gets replaced by NaNs
#These negative values occur through a negative squared bouyancy frequency
eps_Reynolds_bouyancy_grid[eps_Reynolds_bouyancy_grid < 0] = np.nan         
            
def oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_gradient_grid,eps_density_grid):         
    Gamma_Skif_eps_grid = thesis.Skif(eps_Reynolds_bouyancy_grid)
    turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
    #remove negative diffusivity
    turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid<0] = np.nan
    
    #print(np.nanmean(turbulent_diffusivity_Skif_grid),np.nanstd(turbulent_diffusivity_Skif_grid))

    oxygen_flux_Skif_grid = - turbulent_diffusivity_Skif_grid[:,:] * eps_oxygen_gradient_grid
    #convert from m*micromol/(kg*s) to mmol/(m^2*d)
    oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400*(1000/eps_density_grid)
    return oxygen_flux_Skif_grid   


oxygen_flux_Skif_grid = oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_gradient_grid,eps_density_grid)
    


axis[0].plot(np.log10(eps_grid[profile_index,:]),eps_pressure, colors[1], label = "02, d02/dz then to eps")
axis[1].plot(eps_oxygen_gradient_grid[profile_index,:],eps_pressure, colors[1], label = "02, d02/dz then to eps")
axis[1].plot(eps_oxygen_grid[profile_index,:],eps_pressure,"--", c=colors[1])#, label = "02, d02/dz then to eps")
axis[2].plot(oxygen_flux_Skif_grid[profile_index,:],eps_pressure, colors[1], label = "02, d02/dz then to eps")

#-------------------------------------------------------------------
#Interpolation 3 eps to fine grid, then flux calculation on the fine grid
#-------------------------------------------------------------------
results = interpolation3(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, oxygen, oxygen_pressure, number_of_profiles)
    

number_of_profiles,lat,lon,distance = results[0]
interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid,pot_density_grid,eps_grid,N_squared_grid = results[1]

#calculate the viscosity (2 different formula)
viscosity_grid = thesis.get_viscosity(consv_temperature_grid,salinity_grid,density_grid, "default")


#calculate the Reynolds bouyancy number defined on the pressure values defined in eps_pressure
Reynolds_bouyancy_grid = eps_grid/(viscosity_grid*N_squared_grid)

#A negative Reynolds bouyancy has no physical meaning, therefore gets replaced by NaNs
#These negative values occur through a negative squared bouyancy frequency
Reynolds_bouyancy_grid[Reynolds_bouyancy_grid < 0] = np.nan         
            
#apply the correction on the dissipation rate following Garanaik2018           
#corrected_eps_grid = thesis.correct_dissipation(eps_grid,eps_Reynolds_bouyancy_grid)

#from that calculate again the now corrected Bouyancy Reynolds number (in 3D)
#corrected_eps_Reynolds_bouyancy_grid = corrected_eps_grid/(eps_viscosity_grid*eps_N_squared_grid)

#conversion from pressure coordinates to depth
interp_depth = gsw.z_from_p(interp_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west   
#interp_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid)) #shape to the same from as eps_grid
    
oxygen_flux_Osborn_grid = thesis.get_oxygen_flux_osborn(Reynolds_bouyancy_grid,eps_grid,N_squared_grid,oxygen_grid,interp_depth,density_grid)
oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(Reynolds_bouyancy_grid,eps_grid,N_squared_grid,oxygen_grid,interp_depth,density_grid)
oxygen_flux_Skif_grid = thesis.get_oxygen_flux_skif(Reynolds_bouyancy_grid,eps_grid,N_squared_grid,oxygen_grid,interp_depth,density_grid)
    



#axis[0].plot(np.log10(eps_grid[profile_index,:]),interp_pressure, colors[2], label = "eps to 02 and d02/dz")
#axis[1].plot(thesis.central_differences(oxygen_grid[profile_index,:])/thesis.central_differences(interp_depth),interp_pressure, colors[2],label = "eps to 02 and d02/dz")
#axis[1].plot(oxygen_grid[profile_index,:],interp_pressure, "--", c= colors[2])#, label = "eps to 02 and d02/dz")
#axis[2].plot(oxygen_flux_Skif_grid[profile_index,:],interp_pressure, colors[2], label = "eps to 02 and d02/dz")

#-------------------------------------------------------------------
#Interpolation 4: average everything in bins centered around the values of eps_pressure
#-------------------------------------------------------------------
results = interpolation4(lat,lon,distance,pressure,absolute_salinity,consv_temperature, alpha, beta, eps, eps_pressure, oxygen, oxygen_pressure, number_of_profiles)
    

number_of_profiles,lat,lon,distance = results[0]
interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid, pot_density_grid = results[1]
eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid, eps_pot_density_grid = results[2]

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
    



axis[0].plot(np.log10(eps_grid[profile_index,:]),eps_pressure, ":", c = colors[3], label = "bins around eps")
axis[1].plot(thesis.central_differences(eps_oxygen_grid[profile_index,:])/thesis.central_differences(eps_depth),eps_pressure, ":", c = colors[3], label = "bins around eps")
axis[1].plot(eps_oxygen_grid[profile_index,:],eps_pressure, "--", c = colors[3])#, label = "O2 to eps, then do2/dz")
axis[2].plot(oxygen_flux_Skif_grid[profile_index,:],eps_pressure, "--", c =  colors[3], label = "bins around eps")


######################################################################################
axis[0].set_ylabel("depth [dbar]")
axis[0].set_xlabel("log10(eps)")
axis[1].set_xlabel("oxygen")
axis[2].set_xlabel("Shih flux")

axis[2].set_xlim(-200,5)
axis[0].legend(loc = "upper center")
axis[1].legend(loc = "upper center")
axis[2].legend(loc = "upper center")

axis[0].invert_yaxis()
f.suptitle("Impact of the interpolation on oxygen fluxes after Shih")
import pickle
pickle.dump(f, open('interpolation_impact.fig.pickle', 'wb'))

plt.show()



