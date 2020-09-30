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
        #eps_pot_density_grid = np.sort(eps_pot_density_grid, axis = 0)
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
    
    results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid)
    """
    bathymetry                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
    list_of_bathymetry_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
    BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
    list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
    BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
    list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
    """
    bathymetry,list_of_bathymetry_indices = results[0]
    BBL,list_of_BBL_indices = results[1] #not needed here
    #BBL_range,list_of_BBL_range_indices = results[2]
    
    #print(np.asarray(list_of_bathymetry_indices) - np.asarray(list_of_BBL_indices))
    
    #eps_N_grid = np.sqrt(eps_N_squared_grid)
    #ozmidov scale
    #ozmidov_scale_grid = np.sqrt(eps_grid/(eps_N_grid**3))
    
    #conversion from pressure coordinates to depth
    eps_depth = gsw.z_from_p(eps_pressure,np.mean(lat)) #mean lat should be sufficient, because the transect is east-west
    bathymetry_in_m = gsw.z_from_p(bathymetry,np.mean(lat))
    
    eps_depth_grid = np.reshape(eps_depth,(1,-1))*np.ones(np.shape(eps_grid))
    
    distance_from_ground_grid = eps_depth_grid - np.reshape(bathymetry_in_m,(-1,1))
    distance_from_ground_grid[distance_from_ground_grid < 0] = np.nan
    #boundary_check_grid = ~(distance_from_ground_grid < ozmidov_scale_grid)
    
    turbulent_diffusivity_Osborn_grid = thesis.get_turbulent_diffusivity_Osborn(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
    turbulent_diffusivity_Shih_grid = thesis.get_turbulent_diffusivity_Shih(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
    #turbulent_diffusivity_BB_grid = thesis.get_turbulent_diffusivity_BB(eps_Reynolds_bouyancy_grid,eps_grid,eps_N_squared_grid)
      
                 
    shear_velocity_grid = thesis.get_shear_velocity(eps_grid,distance_from_ground_grid)
    
    #compute the mean shear velocity in the BBL as theoretically it should be constant
    shear_velocity = np.zeros(number_of_profiles)
    for profile in range(number_of_profiles):
        BBL_from = int(list_of_BBL_indices[profile])
        BBL_to = int(list_of_bathymetry_indices[profile])
    
        #print(list_of_BBL_indices[profile],list_of_bathymetry_indices[profile])
        BBL_shear_velocity = np.nanmean(shear_velocity_grid[BBL_from:BBL_to])
        
        law_of_the_wall_turbulent_diffusivity = thesis.von_Karman_constant * BBL_shear_velocity * distance_from_ground_grid[profile,BBL_from:BBL_to]
        
        #replace the corresponding bins with the turbulent diffusivity from the law of the wall
        turbulent_diffusivity_Osborn_grid[profile,BBL_from:BBL_to] = law_of_the_wall_turbulent_diffusivity
        turbulent_diffusivity_Shih_grid[profile,BBL_from:BBL_to] = law_of_the_wall_turbulent_diffusivity
    
    oxygen_gradient_grid = thesis.central_differences(eps_oxygen_grid)/thesis.central_differences(eps_depth)
    unit_conversion_grid = 86400*(1000/eps_density_grid) #to convert from m*micromol/(kg*s) to mmol/(m^2*d)

    oxygen_flux_Osborn_grid = - turbulent_diffusivity_Osborn_grid * oxygen_gradient_grid * unit_conversion_grid
    oxygen_flux_Shih_grid = - turbulent_diffusivity_Shih_grid * oxygen_gradient_grid * unit_conversion_grid

    #remove negative diffusivity for the plot
    turbulent_diffusivity_Osborn_grid[turbulent_diffusivity_Osborn_grid==0] = 1e-7
        
    profile_index = np.argmin(np.abs(lon-20.57)) 
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
    #f, = axis2[4].plot(oxygen_flux_Shih_grid[profile_index,:],eps_pressure, c = "k", lw = 2, label = r"Shih O$_2$ flux")
    #f, = axis2[4].plot(oxygen_flux_Osborn_grid[profile_index,:],eps_pressure, ls = "--", c = "k", lw = 2, label = r"Osborn O$_2$ flux")
        
    upper_boundary = np.argmin(np.abs(eps_pressure-52))
    lower_boundary = np.argmin(np.abs(eps_pressure-90))
    
    oxy_index = upper_boundary+np.nanargmin(thesis.central_differences(eps_oxygen_sat_grid[profile_index,upper_boundary:lower_boundary]))
    salt_index = upper_boundary+np.nanargmax(thesis.central_differences(eps_salinity_grid[profile_index,upper_boundary:lower_boundary]))
    temp_index = upper_boundary+np.nanargmax(thesis.central_differences(eps_consv_temperature_grid[profile_index,upper_boundary:lower_boundary]))
        
    axis2[0].axhline(eps_pressure[temp_index], c="k")
    axis2[1].axhline(eps_pressure[salt_index],c="k")
    axis2[2].axhline(eps_pressure[oxy_index],c="k")
    
    halocline_index = int(np.median([oxy_index,salt_index,temp_index]))
    halocline_depth = eps_pressure[halocline_index]
    print(halocline_index)
    halocline_density = eps_pot_density_grid[profile_index,halocline_index] 
    #print(halocline_density)
    interval_start_index = np.nanargmin(np.abs(eps_pot_density_grid[profile_index,:] - (halocline_density - 0.75)))
    interval_stop_index = np.nanargmin(np.abs(eps_pot_density_grid[profile_index,:] - (halocline_density + 0.75)))
        
    #    for datapoint in eps_pressure[halocline_index-20:halocline_index+20]:
    #    print(datapoint 
    
    
    #check if the vertical interval is bigger than the maximum halocline thickness
    while True:
        if halocline_depth - eps_pressure[interval_start_index] <= 10:
            break
        elif interval_start_index >= halocline_index:
            break
        else:
            interval_start_index += 1
            
            
    while True:
        if eps_pressure[interval_stop_index] - halocline_depth <= 10:
            break
        elif interval_stop_index <= halocline_index:
            break
        else:
            interval_stop_index -= 1
                        
    
    #axis2[0].axhspan(eps_pressure[interval_start_index],eps_pressure[interval_stop_index], color = "r",  alpha = 0.5)
    #axis2[1].axhspan(eps_pressure[interval_start_index],eps_pressure[interval_stop_index], color = "r",  alpha = 0.5)
    #axis2[2].axhspan(eps_pressure[interval_start_index],eps_pressure[interval_stop_index], color = "r",  alpha = 0.5)
    #axis2[3].axhspan(eps_pressure[interval_start_index],eps_pressure[interval_stop_index], color = "r",  alpha = 0.5)
    #axis2[4].axhspan(eps_pressure[interval_start_index],eps_pressure[interval_stop_index], color = "r",  alpha = 0.5)
            
    axis2[0].invert_yaxis()
    #axis2[4].set_xlim(-12,1.5)
    
    axis2[0].set_ylabel("pressure [dbar]")
    axis2[0].set_xlabel(r"$\theta$ [$\degree$C]")
    axis2[1].set_xlabel("SA[g/kg]")
    axis2[2].set_xlabel(r"O$_2$ saturation [%]")
    axis2[3].set_xlabel(r"$\epsilon$ [m²/s³]")
    #axis2[4].set_xlabel("OF [mmol/m²d]")
    
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
    fout2.subplots_adjust(top=0.92,bottom=0.17,left=0.122,right=0.946,hspace=0.058,wspace=0.2)
    plt.savefig("/home/ole/Thesis/Analysis/mss/pictures/halocline_profile_"+cruise_name+"_"+transect_name, dpi = 400)
     
plt.show()




