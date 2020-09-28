############################################################
#plots an example profile of a specified transect
##############################################################
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
#matplotlib preferences:
MINI_SIZE = 9
SMALL_SIZE = 10.95
MEDIUM_SIZE = 12
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE, titleweight = "bold")     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MINI_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight = "bold")  # fontsize of the figure title

import gsw
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore')
import mss_functions as thesis



#datafile_path = "/home/ole/windows/processed_mss/emb177/TS1_6.npz"
#datafile_path = "/home/ole/windows/processed_mss/emb217/TR1-4.npz"
paths = ["/home/ole/windows/processed_mss/emb169/TS11.npz","/home/ole/windows/processed_mss/emb177/TS1_6.npz","/home/ole/windows/processed_mss/emb217/TR1-4.npz"]

for datafile_path in paths:
    data = np.load(datafile_path)

    splitted_path = datafile_path.split("/")
    cruise_name = splitted_path[-2]
    transect_name = splitted_path[-1][:-4]

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
        assert np.all(eps_grid[~np.isnan(eps_grid)] > 0)
        corrected_eps_grid = data["corrected_bin_eps_grid"]
        eps_consv_temperature_grid = data["bin_consv_temperature_grid"]
        eps_salinity_grid = data["bin_salinity_grid"]
        eps_oxygen_sat_grid = data["bin_oxygen_sat_grid"]
        eps_oxygen_grid = data["bin_oxygen_grid"] 
        
        eps_N_squared_grid = data["bin_N_squared_grid"]
        eps_density_grid = data["bin_density_grid"]
        eps_pot_density_grid = data["bin_pot_density_grid"]
        #eps_viscosity_grid = data["eps_viscosity_grid"]
        eps_Reynolds_bouyancy_grid = data["bin_Reynolds_bouyancy_grid"]
        corrected_eps_Reynolds_bouyancy_grid = data["corrected_bin_Reynolds_bouyancy_grid"]

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
        raise AssertionError
    
    eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west 
    
    Gamma_Skif_eps_grid = thesis.Skif(eps_Reynolds_bouyancy_grid)
    turbulent_diffusivity_Skif_grid = Gamma_Skif_eps_grid * eps_grid / (eps_N_squared_grid)
    #remove negative diffusivity
    turbulent_diffusivity_Skif_grid[turbulent_diffusivity_Skif_grid==0] = 1e-7
    
    Gamma_Osborn_eps_grid = thesis.Osborn(eps_Reynolds_bouyancy_grid)
    turbulent_diffusivity_Osborn_grid = Gamma_Osborn_eps_grid * eps_grid / (eps_N_squared_grid)

    #remove negative diffusivity 
    turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid==0] = 1e-7
        
    profile_index = np.argmin(np.abs(lon-20.56)) 
    print(cruise_name,transect_name,lon[profile_index])
     
    """   
    def make_patch_spines_invisible(ax):
        ax.set_frame_on(True)
        ax.patch.set_visible(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
    
    fout, axis = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    axis[0].invert_yaxis()
    axis[0].plot(np.log10(turbulent_diffusivity_Osborn_grid[profile_index]),eps_pressure, "--",label = "Osborn")
    axis[0].plot(np.log10(turbulent_diffusivity_Skif_grid[profile_index]),eps_pressure, label = "Shih")
    axis[1].plot(np.log10(eps_Reynolds_bouyancy_grid[profile_index]),eps_pressure, label = "Reb")
    axis[1].plot(np.log10(eps_N_squared_grid[profile_index]),eps_pressure, label = r"N^2")
    
    fout.suptitle(cruise_name)
    axis[0].legend()
    axis[1].legend()
            
    figure, temperature_axis = plt.subplots(nrows = 1, ncols = 1)

    temperature_axis.invert_yaxis()

    salinity_axis = temperature_axis.twiny()
    oxygen_axis = temperature_axis.twiny()

    # Offset the right spine of oxygen_axis.  The ticks and label have already been
    # placed on the right by twinx above.
    oxygen_axis.spines["top"].set_position(("axes", 1.1))
    # Having been created by twinx, oxygen_axis has its frame off, so the line of its
    # detached spine is invisible.  First, activate the frame but make the patch
    # and spines invisible.
    make_patch_spines_invisible(oxygen_axis)
    # Second, show the right spine.
    oxygen_axis.spines["top"].set_visible(True)

    t, = temperature_axis.plot(eps_consv_temperature_grid[profile_index,:],eps_pressure, c = "tab:red", lw = 2, label = "conservative temperature")
    s, = salinity_axis.plot(eps_salinity_grid[profile_index,:],eps_pressure, c = "tab:green", lw = 2, label = "practical salinity")
    o, = oxygen_axis.plot(eps_oxygen_sat_grid[profile_index,:],eps_pressure, c = "tab:blue", lw = 2, label = "oxygen saturation")

      

    tick_spacing = 10
    temperature_axis.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing)) 
        
    temperature_axis.set_xlabel(r"conservative temperature [$\degree$C]")
    temperature_axis.set_ylabel("pressure [dbar]")
    salinity_axis.set_xlabel("absolute salinity [g/kg]")
    oxygen_axis.set_xlabel("oxygen saturation [%]")


    temperature_axis.xaxis.label.set_color(t.get_color())
    salinity_axis.xaxis.label.set_color(s.get_color())
    oxygen_axis.xaxis.label.set_color(o.get_color())

    tkw = dict(size=4, width=1.5)
    temperature_axis.tick_params(axis='x', colors=t.get_color())#, **tkw)
    salinity_axis.tick_params(axis='x', colors=s.get_color())#, **tkw)
    oxygen_axis.tick_params(axis='x', colors=o.get_color())#, **tkw)
    temperature_axis.tick_params(axis='y')#, **tkw)

    #lines = [s,t,o]
    #temperature_axis.legend(lines, [l.get_label() for l in lines])

    if cruise_name == "emb169":
        season = "spring cruise"
    elif cruise_name == "emb177":
        season = "winter cruise"
    elif cruise_name == "emb217":
        season = "summer cruise"
        
    cruise_str = season+": "+cruise_name
    lon_str = str(np.round(lon[profile_index],3))+r" $\degree$E"
    textstr = "\n".join((cruise_str,transect_name,lon_str))
    props = dict(boxstyle='square', facecolor = "white")
    temperature_axis.text(0.5, 0.96, textstr, transform=temperature_axis.transAxes, fontsize=14,va='top', ha = "center", bbox=props, multialignment = "center")
    
    figure.set_size_inches(9,10.5)
    figure.tight_layout()
    """

    fout2, axis2 = plt.subplots(nrows = 1, ncols = 4, sharey = True)
    t, = axis2[0].plot(eps_consv_temperature_grid[profile_index,:],eps_pressure, c = "tab:red", lw = 2, label = "conservative temperature")
    s, = axis2[1].plot(eps_salinity_grid[profile_index,:],eps_pressure, c = "tab:green", lw = 2, label = "practical salinity")
    o, = axis2[2].plot(eps_oxygen_sat_grid[profile_index,:],eps_pressure, c = "tab:blue", lw = 2, label = "oxygen saturation")
    e, = axis2[3].plot(np.log10(eps_grid[profile_index,:]),eps_pressure, c = "k", lw = 2, label = "dissipation rate")
    
    upper_boundary = np.argmin(np.abs(eps_pressure-52))
    lower_boundary = np.argmin(np.abs(eps_pressure-90))
    
    oxy_index = upper_boundary+np.nanargmin(thesis.central_differences(eps_oxygen_sat_grid[profile_index,upper_boundary:lower_boundary]))
    salt_index = upper_boundary+np.nanargmax(thesis.central_differences(eps_salinity_grid[profile_index,upper_boundary:lower_boundary]))
    temp_index = upper_boundary+np.nanargmax(thesis.central_differences(eps_consv_temperature_grid[profile_index,upper_boundary:lower_boundary]))
        
        
    axis2[0].axhline(eps_pressure[temp_index], c="k")
    axis2[1].axhline(eps_pressure[salt_index],c="k")
    axis2[2].axhline(eps_pressure[oxy_index],c="k")
        
    axis2[0].invert_yaxis()
    
    axis2[0].set_ylabel("pressure [dbar]")
    axis2[0].set_xlabel(r"$\theta$ [$\degree$C]")
    axis2[1].set_xlabel("SA[g/kg]")
    axis2[2].set_xlabel(r"O$_2$ saturation [%]")
    axis2[3].set_xlabel(r"$\epsilon$ [m²/s³]")

    width = 6.2012
    height = width / 1.618
    
    if cruise_name == "emb169":
        season = "spring cruise"
    elif cruise_name == "emb177":
        season = "winter cruise"
    elif cruise_name == "emb217":
        season = "summer cruise"
        
    cruise_str = season+" "+cruise_name
    lon_str = str(np.round(lon[profile_index],3))+r" $\degree$E"
    textstr = ", ".join((cruise_str,transect_name,lon_str))
    
    fout2.set_size_inches(width,height)
    fout2.suptitle(textstr)
    fout2.subplots_adjust(top=0.92,bottom=0.17,left=0.122,right=0.974,hspace=0.058,wspace=0.341)
    plt.savefig("/home/ole/Thesis/Analysis/mss/pictures/halocline_profile_"+cruise_name+"_"+transect_name, dpi = 400)
     
plt.show()




