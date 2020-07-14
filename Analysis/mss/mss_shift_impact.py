############################################################
#compares the impact of shifting the measured raw oxygen data by a small amount on the computed oxygen fluxes

#EMB169 is not yet possible as the raw data is saved in two foldes instead of one. (No case handling of that) 
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

import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')


#################################################
#Load raw data and shift the oxygen the defined values
#################################################

FOLDERNAME_raw = "/home/ole/windows/emb217_mss_data" 
FOLDERNAME_shifted =  "/home/ole/windows/shifted_mss_data"
shift_values = [-0.5,-0.25,0.25,0.5,0]
#shift_values = [0.5,0.5,0]



f_flux,flux_axarr = plt.subplots(1)
#flux_axarr.invert_yaxis()


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
    

##################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!main function starts!!!!!!!!!!!!!!!!!!!!!
##################################################################


############################
#loop over the raw data and prepare it with shifted oxygen data points
############################
for shift_index,shift_value in enumerate(shift_values):

    print("-----------------------------------------------")
    print("-----------------------------------------------")
    print("compute flux with a shift of "+str(shift_value)+" m")
    print("-----------------------------------------------")
    print("-----------------------------------------------")

    path = pathlib.Path(FOLDERNAME_raw)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME_raw.split("/")
    
    cruisename = splitted_foldername[4][0:6]
    
    print("load, clean and shift data")
    print(cruisename)    
   
                
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".mat":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME_raw+"/"+DATAFILENAME
        
        #skip this filename, not what thsi program expects
        if DATAFILENAME == "TS11_TODL_merged.mat":
            continue
        if DATAFILENAME[0] == "S":
            continue   

        transect_name = DATAFILENAME[:-4]

        
        results = load_clean_interpolate_and_shift_data(datafile_path, shift_value)
        
        try:
            number_of_profiles,lat,lon,distance = results[0]
            interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid, pot_density_grid = results[1]
            eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid, eps_pot_density_grid = results[2]
        except TypeError:
            print(cruisename,transect_name,"is skipped!")
            continue
 

        #calculate the viscosity (2 different formula)
        eps_wiki_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_salinity_grid,eps_density_grid,"Wikipedia")
        eps_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_salinity_grid,eps_density_grid, "default")

        #print("eps_viscosity_grid\n",np.shape(eps_viscosity_grid))

        #calculate the Reynolds bouyancy number defined on the pressure values defined in eps_pressure
        eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared_grid)
        eps_wiki_Reynolds_bouyancy_grid = eps_grid/(eps_wiki_viscosity_grid*eps_N_squared_grid)


        #A negative Reynolds bouyancy has no physical meaning, therefore gets replaced by NaNs
        #These negative values occur through a negative squared bouyancy frequency
        eps_Reynolds_bouyancy_grid[eps_Reynolds_bouyancy_grid < 0] = np.nan 
        eps_wiki_Reynolds_bouyancy_grid[eps_wiki_Reynolds_bouyancy_grid < 0] = np.nan            
                    
        #apply the correction on the dissipation rate following Garanaik2018           
        corrected_eps_grid = thesis.correct_dissipation(eps_grid,eps_Reynolds_bouyancy_grid)
        corrected_eps_wiki_grid = thesis.correct_dissipation(eps_grid,eps_wiki_Reynolds_bouyancy_grid)

        #from that calculate again the now corrected Bouyancy Reynolds number (in 3D)
        corrected_eps_Reynolds_bouyancy_grid = corrected_eps_grid/(eps_viscosity_grid*eps_N_squared_grid)
        corrected_eps_wiki_Reynolds_bouyancy_grid = corrected_eps_wiki_grid/(eps_wiki_viscosity_grid*eps_N_squared_grid)

        #use self written function to get BBL
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,interp_pressure,density_grid,oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02)
        bathymetrie,list_of_bathymetrie_indices = results[0]
        BBL,list_of_BBL_indices = results[1]
        BBL_range,list_of_BBL_range_indices = results[2]

          
          
 
          
        np.savez("/home/ole/windows/shifted_mss_data/"+cruisename+"_"+transect_name ,number_of_profiles = number_of_profiles, lat = lat,lon = lon,distance = distance, bathymetrie = bathymetrie,list_of_bathymetrie_indices = list_of_bathymetrie_indices,BBL = BBL,list_of_BBL_indices = list_of_BBL_indices,BBL_range = BBL_range,list_of_BBL_range_indices = list_of_BBL_range_indices, interp_pressure = interp_pressure,oxygen_grid = oxygen_grid, oxygen_sat_grid = oxygen_sat_grid, salinity_grid = salinity_grid,consv_temperature_grid = consv_temperature_grid, density_grid = density_grid, pot_density_grid = pot_density_grid, eps_pressure = eps_pressure,eps_grid = eps_grid, corrected_eps_wiki_grid = corrected_eps_wiki_grid, corrected_eps_grid = corrected_eps_grid, eps_salinity_grid = eps_salinity_grid, eps_consv_temperature_grid = eps_consv_temperature_grid, eps_oxygen_grid = eps_oxygen_grid, eps_oxygen_sat_grid = eps_oxygen_sat_grid, eps_N_squared_grid = eps_N_squared_grid,eps_density_grid = eps_density_grid, eps_pot_density_grid = eps_pot_density_grid, eps_Reynolds_bouyancy_grid =  eps_Reynolds_bouyancy_grid, corrected_eps_Reynolds_bouyancy_grid = corrected_eps_Reynolds_bouyancy_grid,  eps_wiki_Reynolds_bouyancy_grid = eps_wiki_Reynolds_bouyancy_grid, corrected_eps_wiki_Reynolds_bouyancy_grid = corrected_eps_wiki_Reynolds_bouyancy_grid)
        
        
               
    ########################################################################################################################################################## 
    rolling_window_size = 12

    maximum_reasonable_flux = 500 #float('Inf') #200 #Fluxes above this value will be discarded
    acceptable_slope = 2 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 

    print("compute rolling average")
         
    print("#######################################")

    number_of_fluxes_over_the_threshold = 0
    number_of_transects = 0
    total_number_of_fluxes = 0
    number_of_zero_flux = 0
    amount_of_missing_values = 0
    total_number_of_valid_profiles = 0    
    total_number_of_correct_profiles = 0
    
    path = pathlib.Path(FOLDERNAME_shifted)
    DATAFILENAMES = []

    #splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    #cruisename = splitted_foldername[-1]
    
    print("cruisename",cruisename)
    
    if cruisename == "emb217":
        upper_bound_halocline_as_density = 1006.4 #1005.75
        lower_bound_halocline_as_density = 1008.5 #1006.25
    elif cruisename == "emb177":
        upper_bound_halocline_as_density = 1006.9 #1006.9
        lower_bound_halocline_as_density = 1008.2 #1007.9   
    elif cruisename == "emb169":
        upper_bound_halocline_as_density = 1006.5 
        lower_bound_halocline_as_density = 1008.6    

    dissipation_list = []
    BB_flux_list = []
    Shih_flux_list = []
    Osborn_flux_list = []
    longitude_list = []
    bathymetry_list = []
    bathymetry_longitude_list = []
    interval_list = []
    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        #print(all_files_name[:6],cruisename)
        if all_files_name[-4:] == ".npz" and all_files_name[:6] == cruisename:
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    print(DATAFILENAMES)
    
    
    for DATAFILENAME in DATAFILENAMES:
    
        datafile_path = FOLDERNAME_shifted+"/"+DATAFILENAME
        
        transect_name = DATAFILENAME[:-4]
    
        #skip the short "S206" transects
        if transect_name[0:4] == "S106":
            print(transect_name,"skipped")
            continue
            
        #something is not correct with this measurement
        if cruisename == "emb169" and transect_name[0:4] == "TS13":
            print(transect_name,"skipped, measurement looks wrong")
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
        
        number_of_transects+=1 
            
        print("Number of profiles:",number_of_profiles)
        
        #print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
        
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
        oxygen_flux_Skif_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0
        lon_without_outliers = []
        bathymetrie_without_outliers = []
        
        transect_oxygen_flux_statistic = []
        

        #compare the pressure value of the lowest valid data point per profile th the neighbouring profiles to determine outliers
        list_of_short_profiles = thesis.get_list_of_short_profiles(number_of_profiles,bathymetrie,acceptable_slope)
        
        #hard coded to accout for one short profile, that my condtions doesn't recognize
        #if the boolean array contains at least one true
        if np.any(lon == 20.466176666666666):
            list_of_short_profiles.append(np.argmax(lon == 20.466176666666666))
        
        for profile in range(number_of_profiles):
        
            #if the current profile is too short, skip it
            if profile in list_of_short_profiles:
                print(str(lon[profile])+": short profile")
                continue
                
            
            if np.nanmean(eps_oxygen_sat_grid[profile]) < 0:
                print(cruisename,transect_name,"negative oxygen values")
                continue
                

                    
            from_index =  np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-upper_bound_halocline_as_density))     
            to_index = np.nanargmin(abs(np.asarray(eps_pot_density_grid[profile])-lower_bound_halocline_as_density))
                     
                
                
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
            
            total_number_of_correct_profiles+=1
                                    
            #if the water colum portion contains only nan values, save only the bathymetrie then skip it
            #useful if the interval to average over is deeper than the current bathymetrie
            if np.all(np.isnan(oxygen_flux_BB_grid[profile,from_index:to_index])): # or  to_index == list_of_bathymetrie_indices[profile]-1:

                #find the correct position in the sorted list
                for index,value in enumerate(bathymetry_longitude_list):
                    if value > lon[profile]:
                        list_position = index
                        break
                    elif index == len(bathymetry_longitude_list)-1:
                        list_position = index+1
                        break
                
                #if the list is empty                   
                if len(bathymetry_longitude_list) == 0:   
                    bathymetry_list.append(bathymetrie[profile])
                    bathymetry_longitude_list.append(lon[profile])
                    interval_list.append([np.nan,np.nan])
                else:
                    bathymetry_list.insert(list_position,bathymetrie[profile])
                    bathymetry_longitude_list.insert(list_position,lon[profile])
                    interval_list.insert(list_position,[np.nan,np.nan])
                continue
   
            
            
                
            #find the correct position in the sorted list
            for index,value in enumerate(bathymetry_longitude_list):
                if value > lon[profile]:
                    list_position = index
                    break
                elif index == len(bathymetry_longitude_list)-1:
                    list_position = len(bathymetry_longitude_list)
                    break
                                
            if len(bathymetry_longitude_list) == 0:   
                bathymetry_list.append(bathymetrie[profile])
                bathymetry_longitude_list.append(lon[profile])
                interval_list.append([eps_pressure[from_index],eps_pressure[to_index]])
            else:
                bathymetry_list.insert(list_position,bathymetrie[profile])
                bathymetry_longitude_list.insert(list_position,lon[profile])
                interval_list.insert(list_position,[eps_pressure[from_index],eps_pressure[to_index]])
                                            
            #find the correct position in the sorted list
            for index,value in enumerate(longitude_list):
                if value > lon[profile]:
                    list_position = index
                    break
                elif index == len(longitude_list)-1:
                    list_position = len(longitude_list)
                    break
            
                 
            if len(longitude_list) == 0:    
                dissipation_list.append(eps_grid[profile,from_index:to_index])
                #BB_flux_list.append(oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.append(oxygen_flux_Skif_grid[profile,from_index:to_index])
                Osborn_flux_list.append(oxygen_flux_Osborn_grid[profile,from_index:to_index])
                longitude_list.append(lon[profile])
            
            
            else:
                
                #Sort the current profile into the list            
                dissipation_list.insert(list_position,eps_grid[profile,from_index:to_index])
                #BB_flux_list.insert(list_position,oxygen_flux_BB_grid[profile,from_index:to_index])
                Shih_flux_list.insert(list_position,oxygen_flux_Skif_grid[profile,from_index:to_index])
                Osborn_flux_list.insert(list_position,oxygen_flux_Osborn_grid[profile,from_index:to_index])
                longitude_list.insert(list_position,lon[profile])

    
            assert(np.all(longitude_list == sorted(longitude_list)))

            total_number_of_valid_profiles+=1

        
    ###########################################################################################################################
    #print(len(longitude_list))
    assert(len(longitude_list) != 0)
    assert(np.all(longitude_list == sorted(longitude_list)))
    assert(np.all(bathymetry_longitude_list == sorted(bathymetry_longitude_list)))        
    interval_list = np.asarray(interval_list)
    
    #compute mean and std over the saved intervals
    mean_Osborn_flux = [None] * total_number_of_valid_profiles
    mean_Shih_flux = [None] * total_number_of_valid_profiles
    median_flux = [None] * total_number_of_valid_profiles
    #upper_percentile_flux = [None] * total_number_of_valid_profiles
    #lower_percentile_flux = [None] * total_number_of_valid_profiles
    #second_upper_percentile_flux = [None] * total_number_of_valid_profiles
    #second_lower_percentile_flux = [None] * total_number_of_valid_profiles
  
    #bathymetrie_mean = [None] * number_of_intervals
    #bathymetrie_percentile = [None] * number_of_intervals
    
    log_mean_dissipation = [None] * total_number_of_valid_profiles
    arith_mean_dissipation = [None] * total_number_of_valid_profiles
    median_dissipation = [None] * total_number_of_valid_profiles
    #lower_percentile_dissip = [None] * total_number_of_valid_profiles
    #upper_percentile_dissip = [None] * total_number_of_valid_profiles

    #second_lower_percentile_dissip = [None] * total_number_of_valid_profiles
    #second_upper_percentile_dissip = [None] * total_number_of_valid_profiles

    """
    mean_dissipation_med = [None] * number_of_profiles
    median_dissipation_med = [None] * number_of_profiles
    #std_dissipation_med = [None] * number_of_profiles
    lower_percentile_dissip_med = [None] * number_of_profiles
    upper_percentile_dissip_med = [None] * number_of_profiles

    second_lower_percentile_dissip_med = [None] * number_of_profiles
    second_upper_percentile_dissip_med = [None] * number_of_profiles
    """
    
    #compute statistical properties of the saved values
    for index in range(total_number_of_valid_profiles):
        temp_Shih_flux = Shih_flux_list[index]
        number_of_fluxes_over_the_threshold += np.sum(np.abs(temp_Shih_flux)>maximum_reasonable_flux)
        number_of_zero_flux += np.sum(np.abs(temp_Shih_flux)==0)
        amount_of_missing_values += np.sum(np.isnan(temp_Shih_flux))
        #count the number of flux data points        
        total_number_of_fluxes += temp_Shih_flux.size
        
        temp_Shih_flux[np.abs(temp_Shih_flux)>maximum_reasonable_flux] = np.nan
    
        mean_Shih_flux[index] = np.nanmean(temp_Shih_flux)
        median_flux[index] = np.nanmedian(temp_Shih_flux)
        
        temp_Osborn_flux = Osborn_flux_list[index]
        #replace fluxes above the threshold wit nan
        temp_Osborn_flux[np.abs(temp_Osborn_flux)>maximum_reasonable_flux] = np.nan  
        mean_Osborn_flux[index] = np.nanmean(temp_Osborn_flux)
        
        #upper_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, flux_percentile)
        #lower_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, 100-flux_percentile)
        #second_upper_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, second_flux_percentile)
        #second_lower_percentile_flux[index] = np.nanpercentile(temp_Shih_flux, 100-second_flux_percentile)
                        
        """        
        mean_min_flux[index] = np.nanmean(oxygen_flux_statistic[index][:,0],axis=0)
        median_min_flux[index] = np.nanmedian(oxygen_flux_statistic[index][:,0],axis=0)

        upper_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], flux_percentile)
        lower_percentile_min_flux[index] = np.nanpercentile(oxygen_flux_statistic[index][:,0], 100-flux_percentile)
        """        
        #bathymetrie_mean[index] = np.nanmean(bathymetrie_statistic[index])
        
        log_mean_dissipation[index] = np.nanmean(np.log10(dissipation_list[index]))
        arith_mean_dissipation[index] = np.log10(np.nanmean(dissipation_list[index]))
        median_dissipation[index] = np.log10(np.nanmedian(dissipation_list[index]))
        #upper_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], dissip_percentile))
        #lower_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], 100-dissip_percentile))
        #second_upper_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], second_dissip_percentile))
        #second_lower_percentile_dissip[index] = np.log10(np.nanpercentile(dissipation_list[index], 100-second_dissip_percentile))
                    
                       
                           
    rolling_mean_Shih_flux = [None] * total_number_of_valid_profiles
    rolling_mean_Osborn_flux = [None] * total_number_of_valid_profiles
    rolling_median_flux = [None] * total_number_of_valid_profiles
    #rolling_upper_percentile_flux = [None] * total_number_of_valid_profiles
    #rolling_lower_percentile_flux = [None] * total_number_of_valid_profiles
    #rolling_second_upper_percentile_flux = [None] * total_number_of_valid_profiles
    #rolling_second_lower_percentile_flux = [None] * total_number_of_valid_profiles
    
    rolling_log_mean_dissipation = [None] * total_number_of_valid_profiles
    rolling_arith_mean_dissipation = [None] * total_number_of_valid_profiles
    rolling_median_dissipation = [None] * total_number_of_valid_profiles
    #rolling_lower_percentile_dissip = [None] * total_number_of_valid_profiles
    #rolling_upper_percentile_dissip = [None] * total_number_of_valid_profiles
    #rolling_second_upper_percentile_dissip = [None] * total_number_of_valid_profiles
    #rolling_second_lower_percentile_dissip = [None] * total_number_of_valid_profiles
            
    max_longitude_gap_index = np.argmax(np.diff(longitude_list))
    max_longitude_gap = np.diff(longitude_list)[max_longitude_gap_index]
    print("LONGITUDE GAP",max_longitude_gap)
    #compute rolling average
    for index in range(total_number_of_valid_profiles):
    
        #controls that the mean is not computed over too distant points
        if max_longitude_gap > 0.001:
            if ((index+rolling_window_size//2) >= max_longitude_gap_index+1) and ((index-rolling_window_size//2)<=max_longitude_gap_index+1):
                
                rolling_mean_Shih_flux[index] = np.nan
                rolling_mean_Osborn_flux[index] = np.nan
                rolling_median_flux[index] = np.nan  
                #rolling_upper_percentile_flux[index] = np.nan
                #rolling_lower_percentile_flux[index] = np.nan
                #rolling_second_upper_percentile_flux[index] = np.nan
                #rolling_second_lower_percentile_flux[index] = np.nan
                                
                rolling_log_mean_dissipation[index] = np.nan
                rolling_arith_mean_dissipation[index] = np.nan
                rolling_median_dissipation[index] = np.nan   
                #rolling_lower_percentile_dissip[index] =  np.nan
                #rolling_upper_percentile_dissip[index] =  np.nan
                #rolling_second_upper_percentile_dissip[index] =  np.nan
                #rolling_second_lower_percentile_dissip[index] =  np.nan
                continue
    
        if rolling_window_size ==1:
            rolling_mean_Shih_flux[index] = mean_Shih_flux[index]
            rolling_mean_Osborn_flux[index] = mean_Osborn_flux[index]
            rolling_median_flux[index] = mean_flux[index]



            #print(index,longitude_list[index],np.round(mean_flux[index],3))
            
        else:
            try:
                
                rolling_mean_Shih_flux[index] = np.nanmean(mean_Shih_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_mean_Osborn_flux[index] = np.nanmean(mean_Osborn_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_median_flux[index] = np.nanmean(median_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_upper_percentile_flux[index] = np.nanmean(upper_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_lower_percentile_flux[index] = np.nanmean(lower_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_second_upper_percentile_flux[index] = np.nanmean(second_upper_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_second_lower_percentile_flux[index] = np.nanmean(second_lower_percentile_flux[index-(rolling_window_size//2):index+rolling_window_size//2])
                
                #print("test",np.nanmean(log_mean_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2]))
                rolling_log_mean_dissipation[index] = np.nanmean(log_mean_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_arith_mean_dissipation[index] = np.nanmean(arith_mean_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2])
                rolling_median_dissipation[index] = np.nanmean(median_dissipation[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_lower_percentile_dissip[index] =  np.nanmean(lower_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_upper_percentile_dissip[index] =  np.nanmean(upper_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_second_upper_percentile_dissip[index] = np.nanmean(second_upper_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
                #rolling_second_lower_percentile_dissip[index] = np.nanmean(second_lower_percentile_dissip[index-(rolling_window_size//2):index+rolling_window_size//2])
    
                #print(index,longitude_list[index],np.round(mean_flux[index-rolling_window_size//2:index+rolling_window_size//2],3))
    

            except (IndexError,ValueError):
                rolling_mean_Shih_flux[index] = np.nan
                rolling_mean_Osborn_flux[index] = np.nan
                rolling_median_flux[index] = np.nan  
                #rolling_upper_percentile_flux[index] = np.nan
                #rolling_lower_percentile_flux[index] = np.nan
                #rolling_second_upper_percentile_flux[index] = np.nan
                #rolling_second_lower_percentile_flux[index] = np.nan
                                
                rolling_log_mean_dissipation[index] = np.nan
                rolling_arith_mean_dissipation[index] = np.nan
                rolling_median_dissipation[index] = np.nan   
                #rolling_lower_percentile_dissip[index] =  np.nan
                #rolling_upper_percentile_dissip[index] =  np.nan
                #rolling_second_upper_percentile_dissip[index] =  np.nan
                #rolling_second_lower_percentile_dissip[index] =  np.nan 
                
    print(cruisename,"total number of transects =",number_of_transects)  
    print("total_number_of profiles",total_number_of_correct_profiles)        
    print("total_number_of_valid_profiles",total_number_of_valid_profiles)     
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
    

    color = ["tab:blue","tab:green","tab:red", "tab:orange","k"]
       
    if shift_value == 0:
        ls = ":"
    else:
        ls = "-"   
                  
    flux_axarr.plot(longitude_list,rolling_mean_Shih_flux, ls, lw = 2.5, zorder = 3, c = color[shift_index], label = "shifted by " + str(shift_value)+" m")#"tab:blue")
    #flux_axarr.plot(longitude_list,rolling_mean_Osborn_flux, ls = "-.", lw = 2.5, zorder = 3, c = color[n], label = str(n))#, label = label_name)
    
        
        
#flux_axarr.set_ylabel("pressure [dbar]")



flux_axarr.legend(loc = "lower left")

#flux_axarr.set_ylim((-30,1))    
         
flux_axarr.set_ylabel(r"oxygen flux [mmol/(m$^2$*d]")


f_flux.set_size_inches(16,10.5)
f_flux.tight_layout() 
f_flux.subplots_adjust(top=0.94)



f_flux.suptitle(cruisename+": Impact of shifting the raw oxygen measurements")

f_flux.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/isopycnal_averaging/"+cruisename+"_shift_comparison")

         
plt.show()
    
    
