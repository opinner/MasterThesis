############################################################
#this program loads all profiles from a cruise, removes outliers
#and makes a scatter plot 

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
import gsw.conversions as gswc
import gsw 
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as plc
import warnings
warnings.filterwarnings('ignore')
    
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169",
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217","/home/ole/windows/processed_mss/emb177"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177"]

vmin = 20.47
vmax = 20.7

list_of_bad_profiles,reasons = np.loadtxt("./data/list_of_bad_profiles.txt", dtype=str, unpack=True)
 
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    f_density, density_axarr = plt.subplots(nrows = 1, ncols = 2, sharey = True) 
    f_oxygen, oxygen_axarr = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruise_name = splitted_foldername[-1]
    
    print(cruise_name)
           
    if cruise_name == "emb217":
        upper_boundary = 1006.4 #1006.10 #1005.75
        lower_boundary = 1008.5 #1006.6 #1006.25
    elif cruise_name == "emb177":
        upper_boundary = 1006.9 #1006.6 #1007.4
        lower_boundary = 1008.2 #1007.6 #1007.9 
    elif cruise_name == "emb169":
        upper_boundary = 1006.5 
        lower_boundary = 1008.6 
            
    density_axarr[0].axhspan(lower_boundary, upper_boundary, alpha=0.3, color='tab:red')
    density_axarr[1].axhspan(lower_boundary, upper_boundary, alpha=0.3, color='tab:red')  
    oxygen_axarr[0].axhspan(lower_boundary, upper_boundary, alpha=0.3, color='tab:red')  
    oxygen_axarr[1].axhspan(lower_boundary, upper_boundary, alpha=0.3, color='tab:red') 
    
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

    minima = 0
    maxima = len(DATAFILENAMES)

    norm = plc.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap = cm.viridis) #cmap=cm.Greys_r)

    profile_count = 0
    
    for transect_index,DATAFILENAME in enumerate(DATAFILENAMES):
    
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        transect_name = DATAFILENAME[:-4]
    
        if "_".join((cruise_name,transect_name)) in list_of_bad_profiles:
            print("_".join((cruise_name,transect_name)),"skipped")
            continue
        
        #if transect_name not in ["TR1-4","TR1-7","TR1-10"]:
        #    continue
        if transect_name not in ["TS1_4","TS1_7","TS1_10"]:
            continue
                                
        print("\n",transect_name)
        
        data = np.load(datafile_path)
        try:
            number_of_profiles = data["number_of_profiles"] #
            lat = data["lat"] #Latitude of the profiles
            lon = data["lon"] #Longitude of the profiles
            distance = data["distance"] #distance from the starting profile (monotonically increasing)
            
            #interp_pressure = data["interp_pressure"]
            #oxygen_grid = data["oxygen_grid"]
            #salinity_grid = data["salinity_grid"]
            #consv_temperature_grid = data["consv_temperature_grid"]
            #density_grid = data["density_grid"]
            
            eps_pressure = data["bin_pressure"]
            eps_grid = data["bin_eps_grid"]
            corrected_eps_grid = data["corrected_bin_eps_grid"]
            eps_consv_temperature_grid = data["bin_consv_temperature_grid"]
            eps_oxygen_sat_grid = data["bin_oxygen_sat_grid"]
            eps_oxygen_grid = data["bin_oxygen_grid"] 
            
            eps_N_squared_grid = data["bin_N_squared_grid"]
            eps_density_grid = data["bin_density_grid"]
            eps_pot_density_grid = data["bin_pot_density_grid"]
            #eps_viscosity_grid = data["eps_viscosity_grid"]
            eps_Reynolds_bouyancy_grid = data["bin_Reynolds_bouyancy_grid"]
            corrected_eps_Reynolds_bouyancy_grid = data["corrected_bin_Reynolds_bouyancy_grid"]
        
        except KeyError:
            print(transect_name," is skipped, Error during loading data")
            continue
            
        #print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
            
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
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid)
        bathymetry,list_of_bathymetrie_indices = results[0]
        
        #conversion from pressure coordinates to depth
        eps_depth = gswc.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west

        #get oxygen fluxes
        oxygen_flux_osborn_grid = thesis.get_oxygen_flux_osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_BB_grid = thesis.get_oxygen_flux_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
        oxygen_flux_Skif_grid = thesis.get_oxygen_flux_skif(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid,eps_oxygen_grid,eps_depth,eps_density_grid)
      
        example_profile_index = [np.argmin(np.abs(lon-20.55)),np.argmin(np.abs(lon-20.57)),np.argmin(np.abs(lon-20.61))]
        
        for profile in range(number_of_profiles):
        
            if "_".join((cruise_name,transect_name,str(profile))) in list_of_bad_profiles:
                print("_".join((cruise_name,transect_name,str(profile))),"skipped")
                continue
                                    
            
            if profile not in example_profile_index:
                pass #continue
            
            color = np.ones((np.shape(eps_pot_density_grid[profile])[0],4)) * mapper.to_rgba(transect_index) #np.ones(np.shape(eps_pot_density_grid[profile]))*lon[profile]
            

            density_axarr[0].scatter(np.log10(eps_grid[profile,:]),eps_pot_density_grid[profile,:], c = color , cmap = "viridis", vmin = vmin , vmax = vmax, marker = ",", s= 1.4)           
            image = density_axarr[1].scatter(oxygen_flux_Skif_grid[profile,:],eps_pot_density_grid[profile,:], c = color, cmap = "viridis", vmin = vmin , vmax = vmax, marker = ",", s= 1.4)
            oxygen_axarr[0].scatter(thesis.central_differences(eps_oxygen_sat_grid[profile,:]/thesis.central_differences(eps_depth)),eps_pot_density_grid[profile,:], c = color, cmap = "viridis", vmin = vmin , vmax = vmax, marker = ",", s= 1.4)
            image_oxy = oxygen_axarr[1].scatter(eps_pressure,eps_pot_density_grid[profile,:], c = color, cmap = "viridis", vmin = vmin , vmax = vmax, marker = ",", s= 1.4)



    def colorbar(mappable):
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        ax = mappable.axes
        fig = ax.figure
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.1)
        #cax.set_label('Temperature / Celsius')
        return fig.colorbar(mappable, cax=cax)
        
        
    print("\n\n\n Results: \n\n")

    oxygen_axarr[0].invert_yaxis()
    oxygen_axarr[0].set_ylabel(r"potential density $\sigma$ [kg/m$^3$]")

    
    oxygen_axarr[0].set_xlabel(r"dO$_2$ /dz [micromol/kg/m]")
    oxygen_axarr[1].set_xlabel("pressure [dbar]")
    colorbar(image_oxy).set_label(r"longitude $[\degree$E]")  
               
    f_oxygen.set_size_inches(18,10.5)
    f_oxygen.tight_layout()
    f_oxygen.subplots_adjust(top=0.95)
    f_oxygen.suptitle(cruise_name+": potential density vs Oxygen gradient and pressure")
    f_oxygen.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruise_name+"_pot_density_scatter_oxy")

    
    colorbar(image).set_label(r"longitude $[\degree$E]")        

    density_axarr[1].set_xlim(-10.5,-3)
    density_axarr[0].invert_yaxis()
    if cruise_name == "emb177":
        density_axarr[1].set_xlim(-300,100)
    else:  
        density_axarr[1].set_xlim(-200,100)
        
    density_axarr[0].set_xlabel(r"log10($\epsilon$) [m$^2$ s$^{-3}$]") 

    density_axarr[1].set_xlabel(r"Shih oxygen flux [mmol/(m$^2$d]")
    density_axarr[0].set_ylabel(r"potential density $\sigma$ [kg/m$^3$]")

               
    f_density.set_size_inches(18,10.5)
    f_density.tight_layout()
    f_density.subplots_adjust(top=0.95)
    f_density.suptitle(cruise_name+": potential density vs Turbulence and Fluxes")
    f_density.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/"+cruise_name+"_pot_density_scatter")
    
plt.show()
   
