#---------------------------------------------------------#
#Plots every mss measurement as a transect 
#plus as an example one profile from that transect
#---------------------------------------------------------#

#TODO Why all the runtime warnings? Workaround "turn off warnings"?
#TODO RuntimeWarning for nan values. Solve by switching to pandas?
#TODO change plots from distance to longitude
#TODO change linecolor when N^2 < 0
#TODO color Re_b after the three classifications
#---------------------------------------------------------#
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw 
import pathlib
import mss_functions as thesis
import warnings
warnings.filterwarnings('ignore') #ignoring warnings, specially numpy warnings

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)



#contains the MSS Data
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/emb217_mss_data","/home/ole/windows/emb177_mss_data/","/home/ole/windows/emb169_mss_data/MSS055/matlab","/home/ole/windows/emb169_mss_data/MSS038/matlab/"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/emb177_mss_data/","/home/ole/windows/emb217_mss_data"]
#LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/emb217_mss_data"]
 
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    cruisename = splitted_foldername[4][0:6]
    
    print(cruisename)    
    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".mat":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        if DATAFILENAME == "TS11_TODL_merged.mat":
            continue
            

        transect_name = DATAFILENAME[:-4]

        #define the pictures
        f1, axarr1 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
        f2, axarr2 = plt.subplots(nrows = 4, ncols = 1, sharex = True, sharey = True)
        f4, axarr4 = plt.subplots(nrows = 1, ncols = 7, sharey = True)#, sharex = True, 


        results = thesis.load_clean_and_interpolate_data(datafile_path)
        
        try:
            number_of_profiles,lat,lon,distance = results[0]
            interp_pressure,oxygen_sat_grid,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid = results[1]
            eps_pressure,eps_oxygen_sat_grid,eps_oxygen_grid,eps_grid,eps_salinity_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid = results[2]
        except TypeError:
            #print(cruisename,transect_name,"is skipped!")
            continue
 

        #calculate the viscosity (2 different formula)
        eps_wiki_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_density_grid,"Wikipedia")
        eps_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid,eps_density_grid)

        #print("eps_viscosity_grid\n",np.shape(eps_viscosity_grid))

        #calculate the Reynolds bouyancy number defined on the pressure values defined in eps_pressure
        eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared_grid)
        eps_wiki_Reynolds_bouyancy_grid = eps_grid/(eps_wiki_viscosity_grid*eps_N_squared_grid)


        #A negative Reynolds bouyancy has no physical meaning, therefore get replaced by NaNs
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

          
          
        ##########################################################################################################################################################  
          
        np.savez("/home/ole/windows/processed_mss/"+cruisename+"/"+transect_name ,number_of_profiles = number_of_profiles, lat = lat,lon = lon,distance = distance, bathymetrie = bathymetrie,list_of_bathymetrie_indices = list_of_bathymetrie_indices,BBL = BBL,list_of_BBL_indices = list_of_BBL_indices,BBL_range = BBL_range,list_of_BBL_range_indices = list_of_BBL_range_indices, interp_pressure = interp_pressure,oxygen_grid = oxygen_grid,salinity_grid = salinity_grid,consv_temperature_grid = consv_temperature_grid, density_grid = density_grid, eps_pressure = eps_pressure,eps_grid = eps_grid, corrected_eps_wiki_grid = corrected_eps_wiki_grid, corrected_eps_grid = corrected_eps_grid, eps_salinity_grid = eps_salinity_grid, eps_consv_temperature_grid = eps_consv_temperature_grid, eps_oxygen_grid = eps_oxygen_grid, eps_oxygen_sat_grid = eps_oxygen_sat_grid, eps_N_squared_grid = eps_N_squared_grid,eps_density_grid = eps_density_grid, eps_Reynolds_bouyancy_grid =  eps_Reynolds_bouyancy_grid, corrected_eps_Reynolds_bouyancy_grid = corrected_eps_Reynolds_bouyancy_grid,  eps_wiki_Reynolds_bouyancy_grid = eps_wiki_Reynolds_bouyancy_grid, corrected_eps_wiki_Reynolds_bouyancy_grid = corrected_eps_wiki_Reynolds_bouyancy_grid)
        
        
        ##########################################################################################################################################################
        #TODO outsource this?
        #Plotting


        #choose which profile get exemplarily plotted
        profile_index = -5
        
        print("Profile at Longitude",lon[profile_index])
        #print(lon)
        #print(np.all(np.diff(lon)>0))
        img4_0 = axarr4[0].plot(oxygen_sat_grid[profile_index,:],interp_pressure, label = "fine grid")
        img4_0a = axarr4[0].plot(eps_oxygen_sat_grid[profile_index,:],eps_pressure, label = "eps grid")
        img4_0b = axarr4[0].plot(0,BBL[profile_index],"Dr")
        img4_0c = axarr4[0].plot(0,bathymetrie[profile_index],"Dg")
        img4_0d = axarr4[0].plot(0,BBL_range[profile_index],"ok")
        
        
        print("Bottom",bathymetrie[profile_index],"BBL",BBL[profile_index],"max BBL",interp_pressure[profile_index],)


        #img4_1 = axarr4[1].plot(salinity_grid[profile_index,:],interp_pressure, label = "fine grid")
        img4_1 = axarr4[1].plot(density_grid[profile_index,:],interp_pressure, label = "fine grid")
        img4_1b = axarr4[1].plot(eps_density_grid[profile_index,:],eps_pressure, label = "eps grid")

        img4_2 = axarr4[2].plot(consv_temperature_grid[profile_index,:],interp_pressure, label = "fine grid")
        img4_2b = axarr4[2].plot(eps_consv_temperature_grid[profile_index,:],eps_pressure, label = "eps grid")

        #img4_3 = axarr4[3].plot(BV_freq_squared_grid_gsw[profile_index,:],mid_point_pressure, label = "fine grid")
        img4_3b = axarr4[3].plot(eps_N_squared_grid[profile_index,:],eps_pressure, label = "eps grid")


        img4_4 = axarr4[4].plot(eps_viscosity_grid[profile_index,:]*10**6,eps_pressure,label = "Ilker")
        img4_4b = axarr4[4].plot(eps_wiki_viscosity_grid[profile_index,:]*10**6,eps_pressure,"--",label = "Wikipedia")
        
        img4_5 = axarr4[5].plot(np.log10(eps_grid[profile_index,:]),eps_pressure, label = "eps grid")
        img4_5a = axarr4[5].plot(np.log10(corrected_eps_grid[profile_index,:]),eps_pressure, label = "corrected/Ilker")     
        img4_5b = axarr4[5].plot(np.log10(corrected_eps_wiki_grid[profile_index,:]),eps_pressure,label = "corrected/Wikipedia")
        
        img4_6a = axarr4[6].plot(np.log10(eps_Reynolds_bouyancy_grid[profile_index,:]),eps_pressure,label = "Ilker")
        img4_6b = axarr4[6].plot(np.log10(eps_wiki_Reynolds_bouyancy_grid[profile_index,:]),eps_pressure,label = "Wikipedia")

        
        axarr4[0].set_xlabel("oxygen [%]")
        axarr4[0].set_ylabel("pressure [dbar]")
        #axarr4[1].set_xlabel("SA")
        axarr4[1].set_xlabel(r"density [kg/m$^3$]")
        axarr4[2].set_xlabel("CT")
        axarr4[3].set_xlabel("N^2")
        axarr4[4].set_xlabel("viscosity / 10^(-6)")
        axarr4[5].set_xlabel("log10(eps)")
        axarr4[6].set_xlabel("Reynolds_bouyancy")
        axarr4[0].legend()
        axarr4[1].legend()
        axarr4[2].legend()
        axarr4[3].legend()
        axarr4[5].legend()
        axarr4[4].legend()
        axarr4[6].legend()

        #axarr4[6].set_xlim([,400])
        axarr4[0].set_title("Measurement 02")
        axarr4[1].set_title(r"Measurement $\rho$") #"Measurement SA")
        axarr4[2].set_title("Measurement CT")
        axarr4[3].set_title("Calculation N^2")
        axarr4[4].set_title("Calculation nu")
        axarr4[5].set_title("Measurement eps")
        axarr4[6].set_title("Calculation Re_b")
        
       
         
        #Plot the data  
        #append the last distance plus the last difference (for plotting all the n profiles we need a distance array of size n+1 
        plotmesh_distance = np.append(distance,2*distance[-1]-distance[-2])
        plotmesh_longitude = np.append(lon,2*lon[-1]-lon[-2])

        cmap_hot = plt.get_cmap('hot_r')
        cmap_hot.set_bad(color = 'lightgrey')

        cmap_RdBu = plt.get_cmap('RdBu_r')
        cmap_RdBu.set_bad(color = 'lightgrey')
 
        img1_0 = axarr1[0].pcolormesh(plotmesh_longitude,interp_pressure,oxygen_grid.T,cmap = cmap_RdBu)
        img1_1 = axarr1[1].pcolormesh(plotmesh_longitude,interp_pressure,salinity_grid.T,cmap = cmap_RdBu)
        img1_2 = axarr1[2].pcolormesh(plotmesh_longitude,interp_pressure,consv_temperature_grid.T,cmap = cmap_RdBu)
        img1_3 = axarr1[3].pcolormesh(plotmesh_longitude,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -10,cmap = cmap_hot)


        img2_0 = axarr2[0].pcolormesh(plotmesh_longitude,interp_pressure,density_grid.T,cmap = cmap_RdBu)
        img2_1 = axarr2[1].pcolormesh(plotmesh_longitude,eps_pressure,eps_N_squared_grid.T,vmin = 0, vmax = 0.015,cmap = cmap_RdBu)
        img2_2 = axarr2[2].pcolormesh(plotmesh_longitude,eps_pressure,np.log10(eps_grid.T), vmax = -7, vmin = -9,cmap = cmap_hot)
        img2_3 = axarr2[3].pcolormesh(plotmesh_longitude,eps_pressure,eps_Reynolds_bouyancy_grid.T, vmin = 0, vmax = 100,cmap = cmap_hot)
        
        
        #draw the calculated layers in the plot    
        axarr1[0].plot(lon,bathymetrie)
        axarr1[0].plot(lon,BBL,"g")
           
        axarr2[0].plot(lon,bathymetrie)
        axarr2[0].plot(lon,BBL,"g")


        axarr1[3].set_xlabel("Longitude")#("distance [km]")
        axarr2[3].set_xlabel("Longitude")
            
        f1.set_size_inches(18,10.5)
        f2.set_size_inches(18,10.5)
        f4.set_size_inches(18,10.5)

        colorbar(img1_0).set_label(r"Oxygen [$\mu mol/kg$]") 
        colorbar(img1_1).set_label("salinity [SA]") 
        colorbar(img1_2).set_label("consv_temperature [C]")
        colorbar(img1_3).set_label(r"log10(dissipation) $[m^2 s^{-3}]$")

        colorbar(img2_0).set_label(r"density [kg/$m^3$]")
        colorbar(img2_1).set_label(r"$N^2$ $[1/s^2]$")
        colorbar(img2_2).set_label(r"log10($\epsilon$) $[m^2 s^{-3}]$")  
        colorbar(img2_3).set_label(r"$Re_b$")

        for i in range(4):
            if cruisename == "emb169":
                axarr1[i].set_ylim((0,160))           
                axarr2[i].set_ylim((0,160))
                assert(np.min(interp_pressure)>0)
                assert(np.max(interp_pressure)<160)
                 
            if cruisename == "emb177":
                axarr1[i].set_ylim((0,135))           
                axarr2[i].set_ylim((0,135))   
                assert(np.min(interp_pressure)>0)
                assert(np.max(interp_pressure)<135)
                                
            if cruisename == "emb217":
                if DATAFILENAME[:4] == "S106":
                    axarr1[i].set_ylim((0,90))           
                    axarr2[i].set_ylim((0,90)) 
                    assert(np.min(interp_pressure)>0)
                    assert(np.max(interp_pressure)<90)                    
                else:
                    axarr1[i].set_ylim((0,160))           
                    axarr2[i].set_ylim((0,160))
                    assert(np.min(interp_pressure)>0)
                    assert(np.max(interp_pressure)<160)                       



        axarr1[0].invert_yaxis()
        axarr2[0].invert_yaxis()
        axarr4[0].invert_yaxis() 

        f1.suptitle(cruisename+" "+DATAFILENAME[:-4]+" Measurements")
        f2.suptitle(cruisename+" "+DATAFILENAME[:-4]+" Calculations")
        f4.suptitle(cruisename+" "+DATAFILENAME[:-4]+" Profile at Longitude "+str(lon[profile_index]))

        f1.tight_layout() 
        f2.tight_layout() 
        f4.tight_layout()      
        
        f4.subplots_adjust(top=0.92) 
        
        
        f1.savefig("./pictures/"+cruisename+"/"+cruisename+"_"+DATAFILENAME[:-4]+"_Measurements")
        f2.savefig("./pictures/"+cruisename+"/"+cruisename+"_"+DATAFILENAME[:-4]+"_calculations")
        f4.savefig("./pictures/"+cruisename+"/"+cruisename+"_"+DATAFILENAME[:-4]+"_profiles")

        #close the pictures after saving
        plt.close(fig = "all")
    #plt.show()


