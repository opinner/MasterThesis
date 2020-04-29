##################################################


#TODO check units of oxygen concentration: micromol per kg to micromol per m^3
####################################################
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import warnings
#warnings.filterwarnings('ignore')

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)
    
    
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/processed_mss/emb217"]#,"/home/ole/share-windows/processed_mss/emb169","/home/ole/share-windows/processed_mss/emb177"]


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
    
    #print(DATAFILENAMES)
    
    for DATAFILENAME in DATAFILENAMES:
    
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        transect_name = DATAFILENAME[:-4]
    
        #if transect_name != "TR1-10":
        #    continue
    
        #define the pictures
        f1, axarr1 = plt.subplots(nrows = 1, ncols = 7, sharey = True)
        f2, axarr2 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
    
    

        
        
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
        
        interp_pressure                 equidistant 1D pressure array between the highest and the lowest measured pressure value
        oxygen_grid                     oxygen concentration in in micromol per kg as a grid (number_of_profiles x len(interp_pressure))
        salinity_grid                   salinity in g/kg as a grid (number_of_profiles x len(interp_pressure)) 
        consv_temperature_grid          conservative temperature in degrees Celsius as a grid (number_of_profiles x len(interp_pressure))
        density_grid                    density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
        
        eps_pressure                    pressure values to the dissipation rate values (the pressure distance between points is bigger than in interp_pressure) 
        eps_grid                        measured dissipation rate values (number_of_profiles x len(eps_pressure))
        eps_consv_temperature_grid      conservative temperature as a grid (number_of_profiles x len(eps_pressure))
        eps_oxygen_grid                 oxygen concentration in micromol per kg as a grid (number_of_profiles x len(eps_pressure))
        eps_oxygen_sat_grid             oxygen saturation in percent as a grid (number_of_profiles x len(eps_pressure))
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
        
        
        """
        Gamma_Osborn_eps_grid = thesis.Osborn(corrected_eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * corrected_eps_grid / (eps_N_squared_grid)
        #TODO: Which differention scheme should I use? 
        #Here I remove the uppermost row of the diffusivity to get the same shape (diff removes one row)
        oxygen_flux_osborn_grid = turbulent_diffusivity_Osborn_grid[:,1:] * np.diff(eps_oxygen_grid)/np.diff(eps_depth)
        #convert from m*micromol/(l*s) to mmol/(m^2*d)
        oxygen_flux_osborn_grid = oxygen_flux_osborn_grid*86400/(1000**2)      
        
        Gamma_BB_eps_grid = thesis.BB(corrected_eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_BB_grid = Gamma_BB_eps_grid * corrected_eps_grid / (eps_N_squared_grid)
        #TODO: Which differention scheme should I use? 
        #Here I remove the uppermost row of the diffusivity to get the same shape (diff removes one row)
        oxygen_flux_BB_grid = turbulent_diffusivity_BB_grid[:,1:] * np.diff(eps_oxygen_grid)/np.diff(eps_depth)
        #convert from m*micromol/(l*s) to mmol/(m^2*d)
        oxygen_flux_BB_grid = oxygen_flux_BB_grid*86400/(1000**2)
        
        Gamma_Skif_eps_grid = thesis.Skif(corrected_eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * corrected_eps_grid / (eps_N_squared_grid)
        #TODO: Which differention scheme should I use? 
        #Here I remove the uppermost row of the diffusivity to get the same shape (diff removes one row)
        oxygen_flux_Skif_grid = turbulent_diffusivity_Skif_grid[:,1:] * np.diff(eps_oxygen_grid)/np.diff(eps_depth)
        #convert from m*micromol/(l*s) to mmol/(m^2*d)
        oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400/(1000**2)          
        """
        
        Gamma_Osborn_eps_grid = thesis.Osborn(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity 
        turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid<0] = np.nan
        oxygen_flux_osborn_grid = - turbulent_diffusivity_Osborn_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_osborn_grid = oxygen_flux_osborn_grid*86400*(1000/eps_density_grid)        
        
        Gamma_BB_eps_grid = thesis.BB(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_BB_grid = Gamma_BB_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity    
        turbulent_diffusivity_BB_grid[turbulent_diffusivity_BB_grid<0] = np.nan
        oxygen_flux_BB_grid = - turbulent_diffusivity_BB_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_BB_grid = oxygen_flux_BB_grid*86400*(1000/eps_density_grid)
        
        Gamma_Skif_eps_grid = thesis.Skif(eps_Reynolds_bouyancy_grid)
        turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
        #remove negative diffusivity
        turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid<0] = np.nan
        oxygen_flux_Skif_grid = - turbulent_diffusivity_Skif_grid[:,:] * thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
        #convert from m*micromol/(kg*s) to mmol/(m^2*d)
        oxygen_flux_Skif_grid = oxygen_flux_Skif_grid*86400*(1000/eps_density_grid)
        
        
        #plotting a example profile
        #axarr1[0].plot(eps_oxygen_grid[-5,:],eps_pressure) #distance_from_ground_grid[-5,:],eps_pressure)
        #axarr1[0].axvline(x=0)
        
        mask = ma.masked_where(boundary_check_grid[-5,:],distance_from_ground_grid[-5,:])
        #print(mask)
        
        axarr1[0].plot(np.log10(eps_N_squared_grid[-5,:]),eps_pressure)
        
        """
        for i in range(distance_from_ground_grid[-5,:].size):
            print(distance_from_ground_grid[-5,i],ozmidov_scale_grid[-5,i],boundary_check_grid[-5,i])
            
        print("Boundary effects:",np.any(boundary_check_grid[-5,:]))
        """    
            
        #axarr1[0].plot(mask,eps_pressure,"r")
        axarr1[1].plot(np.log10(eps_grid[-5,:]),eps_pressure, label = "dissipation")
        
        #axarr1[2].plot(ozmidov_scale_grid[-5,:],eps_pressure, label = "Ozmidov") #        axarr1[1].plot(eps_N_squared_grid[-5,:],eps_pressure)

        axarr1[2].plot(eps_Reynolds_bouyancy_grid[-5,:],eps_pressure, label = "raw")
        axarr1[2].plot(corrected_eps_Reynolds_bouyancy_grid[-5,:],eps_pressure, label = "corrected")
        
        axarr1[3].plot(Gamma_Osborn_eps_grid[-5,:],eps_pressure,label = "Osborn")
        axarr1[3].plot(Gamma_BB_eps_grid[-5,:],eps_pressure,label = "BB")
        axarr1[3].plot(Gamma_Skif_eps_grid[-5,:],eps_pressure,label = "Skif")
        
        axarr1[4].plot(np.log10(turbulent_diffusivity_Osborn_grid[-5,:]),eps_pressure, label = "Osborn")     
        axarr1[4].plot(np.log10(turbulent_diffusivity_BB_grid[-5,:]),eps_pressure, label = "BB")   
        axarr1[4].plot(np.log10(turbulent_diffusivity_Skif_grid[-5,:]),eps_pressure, label = "Skif")


        
        #axarr1[4].set_xlim((np.nanmin(turbulent_diffusivity_Skif_grid[-5,:]),np.nanmax(turbulent_diffusivity_Skif_grid[-5,:])))        
        axarr1[4].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        
        
        axarr1[5].plot(np.diff(eps_oxygen_grid[-5,:])/np.diff(eps_depth),eps_pressure[1:])
                                
        axarr1[6].plot(oxygen_flux_osborn_grid[-5,:],eps_pressure[:], label = "Osborn")
        axarr1[6].plot(oxygen_flux_BB_grid[-5,:],eps_pressure[:], label = "BB")
        axarr1[6].plot(oxygen_flux_Skif_grid[-5,:],eps_pressure[:], label = "Skif")        
        axarr1[6].set_xlim(-20,20)
                
        #axarr1[6].set_xlim((np.nanmin(oxygen_flux_Skif_grid[-5,:])-0.1*np.nanmax(oxygen_flux_Skif_grid[-5,:]),1.1*np.nanmax(oxygen_flux_Skif_grid[-5,:])))
        axarr1[6].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        
        axarr1[2].legend()
        axarr1[3].legend()
        axarr1[4].legend()        
        
        axarr1[0].set_xlabel(r"$N^2$")#oxygen [micromol/l]")
        axarr1[0].set_ylabel("pressure [dbar]")
        axarr1[1].set_xlabel("Dissipation")#(r"$N^2$ [$1/s^2$]")
        axarr1[2].set_xlabel(r"Reb")
        axarr1[3].set_xlabel(r"$\Gamma$")
        axarr1[4].set_xlabel("turbulent diffusitivity")
        axarr1[5].set_xlabel(r"d 0$_2$ / dz [micromol/(kg m)]")
        axarr1[6].set_xlabel(r"FO [mmol/(m$^2$*d]")

        #axarr1[3].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
        axarr1[0].set_ylim(0,125)
        
        #plotting the transects
        plotmesh_longitude = np.append(lon,2*lon[-1]-lon[-2])

        cmap_hot = plt.get_cmap('hot_r')
        cmap_hot.set_bad(color = 'lightgrey')

        cmap_RdBu = plt.get_cmap('RdBu_r')
        cmap_RdBu.set_bad(color = 'lightgrey')
        
        img2_0 = axarr2[0].pcolormesh(plotmesh_longitude,interp_pressure,oxygen_grid.T,cmap = cmap_RdBu)
        img2_1 = axarr2[1].pcolormesh(plotmesh_longitude,interp_pressure,density_grid.T,cmap = cmap_RdBu)
        img2_2 = axarr2[2].pcolormesh(plotmesh_longitude,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -9.5, cmap = cmap_hot)
        img2_3 = axarr2[3].pcolormesh(plotmesh_longitude,eps_pressure,oxygen_flux_Skif_grid.T, cmap = cmap_RdBu, vmin = -20, vmax= 20)  
        
        colorbar(img2_0).set_label("Oxygen [micromol/kg]") 
        colorbar(img2_1).set_label(r"Density [kg/m$^3$]") 
        colorbar(img2_2).set_label(r"log10(dissipation) $[m^2 s^{-3}]$")
        colorbar(img2_3).set_label(r"FO [mmol/(m$^2$*d]")
        
        print(transect_name)
        
        axarr1[0].invert_yaxis()
        axarr2[0].invert_yaxis()

        f1.suptitle(cruisename+" "+DATAFILENAME[:-4]+" exemplary profiles")
        f2.suptitle(cruisename+" "+DATAFILENAME[:-4]+" Transects")
        

        f1.set_size_inches(18,10.5)
        f2.set_size_inches(18,10.5)

        f1.tight_layout() 
        f2.tight_layout()    
        
        f1.savefig("/home/ole/Thesis/Analysis/mss/pictures/oxygen_fluxes/"+cruisename+"_"+transect_name+"_profiles", DPI = 300)
        f2.savefig("/home/ole/Thesis/Analysis/mss/pictures/oxygen_fluxes/"+cruisename+"_"+transect_name+"_transects", DPI = 300)
        

        plt.show()
        
        #close the pictures after saving
        plt.close(fig = "all")
