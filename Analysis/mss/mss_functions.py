#######################################################################
#functions
#TODO better documentation
#######################################################################
   
   
import numpy as np   
   
def load_clean_and_interpolate_data(datafile_path):
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
    
    #checks if the shape of the data is as expected, if its smaller append nan
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
            oxygen_sat_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),oxygen_sat[i].flatten(), left = np.nan, right = np.nan)
            eps_oxygen_sat_grid[i] = np.interp(eps_pressure,pressure[i].flatten(),oxygen_sat[i].flatten(), left = np.nan, right = np.nan)       
        elif cruisename == "emb177":
            oxygen_grid[i] = np.interp(interp_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i],left=np.nan,right=np.nan)
            eps_oxygen_grid[i] = np.interp(eps_pressure,emb177_oxygen_pressure[i],emb177_oxygen[i].flatten(), left = np.nan, right = np.nan)  
        elif cruisename == "emb169": 
            oxygen_grid[i] = np.interp(interp_pressure,emb169_oxygen_pressure[i],emb169_oxygen[i],left=np.nan,right=np.nan)
            eps_oxygen_grid[i] = np.interp(eps_pressure,emb169_oxygen_pressure[i],emb169_oxygen[i].flatten(), left = np.nan, right = np.nan)  
            
        
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
        oxygen_grid = oxygen_saturation_to_concentration(oxygen_sat_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
        eps_oxygen_grid = oxygen_saturation_to_concentration(eps_oxygen_sat_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid) 
              
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
        oxygen_sat_grid = oxygen_concentration_to_saturation(oxygen_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
        eps_oxygen_sat_grid = oxygen_concentration_to_saturation(eps_oxygen_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid)  
    
    elif cruisename == "emb169":
        #convert oxygen concentration to oxygen saturation (functions are selfwritten but use the gsw toolbox, for more informations see the function (also in this file))
        oxygen_sat_grid = oxygen_concentration_to_saturation(oxygen_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)
        eps_oxygen_sat_grid = oxygen_concentration_to_saturation(eps_oxygen_grid,eps_salinity_grid, eps_consv_temperature_grid, eps_pressure_grid, eps_lat_grid, eps_lon_grid)  
    
    else:
        raise AssertionError
        
        
          
    return [[number_of_profiles,lat,lon,distance],[interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid,pot_density_grid],[eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid,eps_pot_density_grid]]
        

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def oxygen_saturation_to_concentration(oxygen_sat_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat, lon):
    """
    tranforms oxygen saturation in % to oxygen concentration in micro-moles per kg
    """

    import gsw 
    
    """
    try:
        assert not np.any(oxygen_sat_grid<0)
    except AssertionError:
        pass #TODO
        #print(oxygen_sat_grid[oxygen_sat_grid<0])
        #raise AssertionError 
    """
        
    maximum_concentration = gsw.O2sol(salinity_grid, consv_temperature_grid, pressure_grid, lat, lon)


    #print(maximum_concentration[maximum_concentration<0])
    #assert(not np.any(maximum_concentration<0))
    
    
    oxygen_concentration = (oxygen_sat_grid / 100 ) * maximum_concentration
    
    
    #assert(not np.any(oxygen_sat_grid<0))
    #assert(np.all(oxygen_concentration>0))

    return oxygen_concentration


def oxygen_concentration_to_saturation(oxygen_grid,salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid):
    """
    tranforms oxygen concentration in micro-moles per kg to oxygen saturation
    """

    import gsw 
    maximum_concentration_grid = gsw.O2sol(salinity_grid, consv_temperature_grid, pressure_grid, lat_grid, lon_grid)

    oxygen_sat_grid = 100 * oxygen_grid/maximum_concentration_grid

    return oxygen_sat_grid

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def find_bottom_and_bottom_currents(number_of_profiles,interp_pressure,density_grid,oxygen_sat_grid,height_above_ground = 10,minimal_density_difference = 0.02):
    import numpy as np

    """
    procedural function    
    
    searches for the highermost nanvalues as the ground
    calculates the position of the bottom boundary layer    
    
    input:
       number_of_profiles               number of profiles/casts in the transect
       interp_pressure
       
       density_grid                     density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
       oxygen_sat_grid                      oxygen concentration in percent as a grid (number_of_profiles x len(interp_pressure))
       
       height_above_ground              Default value 10m
       minimal_density_difference       Default value 0.02 kg/m^3
    
    
    return values:
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
    
    
    """    
    #search for bottom currents
    ###########################################################################################################################################################
    bathymetrie = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_bathymetrie_indices = np.zeros(number_of_profiles)

    halocline = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_halocline_indices = np.zeros(number_of_profiles)

    BBL = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_BBL_indices = np.zeros(number_of_profiles)

    BBL_range = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_BBL_range_indices = np.zeros(number_of_profiles)

    for i in range(number_of_profiles):

        #------------------search for bathymetrie values starts from below:-------------------------------
        #search is done in the fine grid

        #returns the pressure of the last nan value in a continuous row starting from high pressure (TODO:is it better to use the last index with data?)
        nan_index =  -np.argmax(np.flip(~np.isnan(density_grid[i,:]))) #at the moment the index is negative
        nan_index = density_grid[i,:].size + nan_index #now defined as positive index
        
       
        if nan_index == density_grid[i,:].size:
            if not np.isnan(density_grid[i,-1]): #if there are no NAN values towards the bottom
                nan_index = len(interp_pressure)-1 #set the last index as the index of the bottom
                
        list_of_bathymetrie_indices[i] = nan_index 
        bathymetrie[i] = interp_pressure[nan_index]
        
      
        #TODO
        #------------------search for halocline values starts from above:-------------------------------
        


        #------------------search for BBL values starts from below:-------------------------------  
        
        #index of maximal distance bottom plus 15m 
        
        BBL_boundary_index = np.argmax(interp_pressure >= (bathymetrie[i]-height_above_ground))
        assert(interp_pressure[BBL_boundary_index]<bathymetrie[i]) #tests if the point 15m above the ground is really above
        
        #TODO get the index (and from that the pressure) where the density difference is bigger than 0.01
        #BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:density_grid[i,:].size + nan_index])>0.01))
        
        #get the index (and from that the pressure) where the density difference is at maximum (in the lowermost 15 m)
        assert(nan_index>=0)
        BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:nan_index]))) -1 
        
        density_jump = np.diff(density_grid[i,:])[BBL_index-1]
        #print(density_jump)
        #print(np.max(np.diff(density_grid[i,BBL_boundary_index:nan_index])))
        #check if density jump is really the biggest density step
        #assert(density_jump == np.max(np.diff(density_grid[i,BBL_boundary_index:nan_index])))
        
        #check if the maximum is at the edge of the intervall or 
        if (BBL_index == BBL_boundary_index) or (BBL_index == (BBL_boundary_index+1)):
            #print(nan_index,BBL_index, "at the upper edge")
            BBL_index = nan_index

        #check if the maximum is too small
        elif (density_jump < minimal_density_difference):
            #print(nan_index,BBL_index, "maximum to small")
            BBL_index = nan_index
    
                    
        #print(BBL_index,nan_index)
        #print("BBL",interp_pressure[BBL_index])
        #print(bathymetrie[i])
       
        list_of_BBL_indices[i] = BBL_index 
        BBL[i] = interp_pressure[BBL_index]
        
        list_of_BBL_range_indices[i] = BBL_boundary_index 
        BBL_range[i] = interp_pressure[BBL_boundary_index]
        

                
    return [[bathymetrie,list_of_bathymetrie_indices],[BBL,list_of_BBL_indices],[BBL_range,list_of_BBL_range_indices]]

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def get_list_of_short_profiles(number_of_profiles,bathymetrie,acceptable_slope):
    list_of_short_profiles = []
    count_of_short_profiles = 0
            
    for profile in range(number_of_profiles):
        #check if the profile was stopped too early by comparing it to the predecessor and succesor. If yes, skip it
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
                
    return list_of_short_profiles          
                    
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def central_differences(array):
    dimensions = len(np.shape(array))
    diff_grid = np.zeros(np.shape(array))
     
    #compute the central differences in the 1D array and set the first and last value of the derivative array to NaN    
    if dimensions == 1:
        vertical_points = np.asarray(array).size
        for i in range(1,vertical_points-1):
            diff_grid[i] = (array[i+1] - array[i-1])/2

        diff_grid[0] = np.nan 
        diff_grid[-1] = np.nan
     
    #compute the central differences in every row of the 2D array and set the first and last value of each row of the derivative array to NaN       
    elif dimensions == 2:
        number_of_profiles = np.shape(array)[0]
        vertical_points = np.shape(array)[1]
        
        for profile in range(number_of_profiles):
            for i in range(1,vertical_points-1):
                diff_grid[profile,i] = (array[profile,i+1] - array[profile,i-1])/2

        diff_grid[:,0] = np.nan 
        diff_grid[:,-1] = np.nan
           
    else:
        raise ValueError  
        
    return diff_grid       
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def get_viscosity(T, S, eps_density_grid, formula = "default"):
    """
    returns temperature and density dependent viscosity 
    depending on the keyword it use different formula
    
    input:
        temperature in C
        salinity in g/kg
    
    output:
        dynamic viscosity kg/(m s)
    """
    #dynamic viscosity:
    if formula == "default":
        """
        source: sw_viscosity.m from Lars Umlauf, based on:
    
        %=========================================================================        
        %   2009-12-18: Mostafa H. Sharqawy (mhamed@mit.edu), MIT
        %               - Initial version
        %   2012-06-06: Karan H. Mistry (mistry@mit.edu), MIT
        %               - Allow T,S input in various units
        %               - Allow T,S to be matrices of any size
        %
        % DISCLAIMER:
        %   This software is provided "as is" without warranty of any kind.
        %   See the file sw_copy.m for conditions of use and licence.
        %
        % REFERENCES:
        %   [1] M. H. Sharqawy, J. H. Lienhard V, and S. M. Zubair, Desalination
        %       and Water Treatment, 16, 354-380, 2010. (http://web.mit.edu/seawater/)
        %   [2] B. M. Fabuss, A. Korosi, and D. F. Othmer, J., Chem. Eng. Data 14(2), 192, 1969.
        %   [3] J. D. Isdale, C. M. Spence, and J. S. Tudhope, Desalination, 10(4), 319 - 328, 1972
        %   [4] F. J. Millero, The Sea, Vol. 5, 3  80, John Wiley, New York, 1974
        %   [5] IAPWS release on the viscosity of ordinary water substance 2008
        %=========================================================================
        """
        
        
        
        #check the function range
        assert(np.nanmax(T) < 180 and np.nanmin(T) > 0)
        assert(np.nanmax(S) < 150 and np.nanmin(S) > 0)
        
        S = S/1000;

        a = [1.5700386464E-01,6.4992620050E+01,-9.1296496657E+01,4.2844324477E-05,1.5409136040E+00,1.9981117208E-02,-9.5203865864E-05,7.9739318223E+00,-7.5614568881E-02,4.7237011074E-04]

        mu_w = a[3] + 1./(a[0]*(T+a[1])**2+a[2]);


        A  = a[4] + a[5] * T + a[6] * T**2;
        B  = a[7] + a[8] * T + a[9] * T**2;
        mu = mu_w*(1 + A*S + B*S**2);
    
        kinematic_viscosity = mu/eps_density_grid
        
        assert(np.shape(kinematic_viscosity) == np.shape(T))
        kinematic_viscosity_without_nans = kinematic_viscosity[~np.isnan(kinematic_viscosity)]
        #print(mu_without_nans)
        assert(np.all(np.log10(kinematic_viscosity_without_nans) < -5) and np.all(np.log10(kinematic_viscosity_without_nans) > -7))
    
        return kinematic_viscosity 
    
    if formula == "Ilker":
        """
        Viscosity function I got from Peter (probable source: Ilker)
        """
        #% T is temperature in degC; vis in m2/s
        #% vis=(1.792747-(.05126103*T)+(0.0005918645*T*T))*1e-6;
        #% Ilker
        return (1.792747-(0.05126103*T)+(0.0005918645*T*T))*1e-6



    if formula == "Wikipedia":
        """
        Viscosity function I copied from Wikipedia
        """

        A = 29.39*1e-3
        B = 507.88 
        C = 149.3
        
        return (A * np.exp(B/(T-C)))/eps_density_grid

    
   
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#TURBULENCE PARAMETRIZATIONS




#vectorized version (acts on every element of an array)
def BB(Reb):
    """
    returns vectorized version of the Bouffard and Boegman (2013) turbulence parametrizaion
    """

    #basic version
    def basic_BB(Reb):
        Pr = 700

        if Reb < 0.18:
            return 0
            
        elif ((Reb >= 0.18) and (Reb <= 96.56)):
            return 0.1 * Pr**(-0.5) * Reb**(0.5)

        else:
            return 2*Reb**(-0.5)



    vBB = np.vectorize(basic_BB,otypes=[float])
    return vBB(Reb)





#vectorized version (acts on every element of an array)
def Skif(Reb):
    """
    vectorized version of the Shih(2005) turbulence parametrizaion
    """

    def basic_Skif(Reb):
        if Reb < 7:
            return 0
            
        elif ((Reb >= 7) and (Reb<=100)):
            return 0.2
            
        else:
            return 2*Reb**(-0.5)


    vSkif = np.vectorize(basic_Skif, otypes=[float]) 
    return vSkif(Reb)
     
        

        
#vectorized version (acts on every element of an array)
def Osborn(Reb):
    """
    vectorized version of the Osborn(1980) turbulence parametrizaion
    """
    
    def basic_Osborn(Reb):
        return 0.2
    
    vOsborn = np.vectorize(basic_Osborn, otypes=[float])    
    return vOsborn(Reb)
    
    
    
    
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
"""
Dissipitation correction  
"""

def correct_dissipation(dissipation_1D,Reb_1D, a = 1.3):
    """
    Correction of Dissipation rate correction following Garanaik2018
    """
    #Reb_1D[Reb_1D < 0] = np.nan
    
    correction_factor = (1-np.exp(-a*np.log10(Reb_1D)))
    #assert((np.nanmax(correction_factor) <= 1) and (np.nanmin(correction_factor)>= 0))
    
    dissipation_3D =  dissipation_1D * (1-np.exp(-a*np.log10(Reb_1D)))
    
    return dissipation_3D
    
    
    
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    


def get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid):
    Gamma_Osborn_eps_grid = Osborn(eps_Reynolds_bouyancy_grid)
    turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)

    #remove negative diffusivity 
    turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid<0] = np.nan
    
    #print(np.nanmean(turbulent_diffusivity_Osborn_grid),np.nanstd(turbulent_diffusivity_Osborn_grid))
        
    oxygen_flux_osborn_grid = - turbulent_diffusivity_Osborn_grid[:,:] * central_differences(eps_oxygen_grid)/central_differences(eps_depth)
    #convert from m*micromol/(kg*s) to mmol/(m^2*d)
    oxygen_flux_osborn_grid = oxygen_flux_osborn_grid*86400*(1000/eps_density_grid)        
    return oxygen_flux_osborn_grid

def get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid):        
    Gamma_BB_eps_grid = BB(eps_Reynolds_bouyancy_grid)
    turbulent_diffusivity_BB_grid = Gamma_BB_eps_grid * eps_grid / (eps_N_squared_grid)
    #remove negative diffusivity    
    turbulent_diffusivity_BB_grid[turbulent_diffusivity_BB_grid<0] = np.nan
    oxygen_flux_BB_grid = - turbulent_diffusivity_BB_grid[:,:] * central_differences(eps_oxygen_grid)/central_differences(eps_depth)
    #convert from m*micromol/(kg*s) to mmol/(m^2*d)
    oxygen_flux_BB_grid = oxygen_flux_BB_grid*86400*(1000/eps_density_grid)
    return oxygen_flux_BB_grid
    
def get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid):         
    Gamma_Skif_eps_grid = Skif(eps_Reynolds_bouyancy_grid)
    turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
    #remove negative diffusivity
    turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid<0] = np.nan
    
    #print(np.nanmean(turbulent_diffusivity_Skif_grid),np.nanstd(turbulent_diffusivity_Skif_grid))

    oxygen_flux_Skif_grid = - turbulent_diffusivity_Skif_grid[:,:] * central_differences(eps_oxygen_grid)/central_differences(eps_depth)
    #convert from m*micromol/(kg*s) to mmol/(m^2*d)
    oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400*(1000/eps_density_grid)
    return oxygen_flux_Skif_grid   
        

#####################################################################################################
def colorbar(mappable, ax = None):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    if ax == None:
        ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)
