import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw 
import pathlib
import mss_functions as thesis


LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/processed_mss/emb169","/home/ole/share-windows/processed_mss/emb177","/home/ole/share-windows/processed_mss/emb217"]


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
    
    print(DATAFILENAMES)
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
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
        
        interp_pressure                 equidistant pressure array between the highest and the lowest measured pressure value
        oxygen_grid                     oxygen saturation in percent as a grid (number_of_profiles x len(interp_pressure))
        salinity_grid                   salinity in g/kg as a grid (number_of_profiles x len(interp_pressure)) 
        consv_temperature_grid          conservative temperature in degrees Celsius as a grid (number_of_profiles x len(interp_pressure))
        density_grid                    density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
        
        eps_pressure                    pressure values to the dissipation rate values (the pressure distance between points is bigger than in interp_pressure) 
        eps_grid                        measured dissipation rate values (number_of_profiles x len(eps_pressure))
        eps_consv_temperature_grid      conservative temperature as a grid (number_of_profiles x len(eps_pressure))
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
        
        
        
        Gamma_Osborn_eps_grid = thesis.Osborn(corrected_eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Osborn = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)
        #TODO: Which differention scheme should I use? 
        #Here I remove the uppermost row of the diffusivity to get the same shape (diff removes one row)
        oxygen_flux_osborn = turbulent_diffusivity_Osborn[:,1:] * np.diff(eps_oxygen_grid)/np.diff(eps_pressure)
        
        
        Gamma_BB_eps_grid = thesis.BB(corrected_eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_BB = Gamma_BB_eps_grid * eps_grid / (eps_N_squared_grid)
        #TODO: Which differention scheme should I use? 
        #Here I remove the uppermost row of the diffusivity to get the same shape (diff removes one row)
        oxygen_flux_BB = turbulent_diffusivity_BB[:,1:] * np.diff(eps_oxygen_grid)/np.diff(eps_pressure)
        
        
        Gamma_Skif_eps_grid = thesis.Skif(corrected_eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Skif = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
        #TODO: Which differention scheme should I use? 
        #Here I remove the uppermost row of the diffusivity to get the same shape (diff removes one row)
        oxygen_flux_Skif = turbulent_diffusivity_Skif[:,1:] * np.diff(eps_oxygen_grid)/np.diff(eps_pressure)
        
        
        
        
        
        
        
        
        
