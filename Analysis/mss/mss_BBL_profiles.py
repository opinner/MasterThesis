#---------------------------------------------------------------------#
#Plots 5 succesive profiles with markers where the BBL was determined
#additionally shows where on the transect these 5 profiles lay
#TODO add new grid
#---------------------------------------------------------------------#

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw 
import pathlib

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)

    
#Constants
rho_0 = 1000 #kg/m³
g = 9.81 #m/s² #could be replace by gsw.grav(lat,p)

#contains the MSS Data
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/emb217_mss_data"]#,"/home/ole/share-windows/emb177_mss_data/","/home/ole/share-windows/emb169_mss_data/MSS055/matlab/","/home/ole/share-windows/emb169_mss_data/MSS038/matlab/"]

 
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
          
        if DATAFILENAME[:2] !=  "TR":
            continue

        #define the pictures
        f1, axarr1 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)
        f2, axarr2 = plt.subplots(nrows = 1, ncols = 10, sharey = True)
        #f4, axarr4 = plt.subplots(nrows = 1, ncols = 7, sharey = True)#, sharex = True, 

    
        results = thesis.load_clean_and_interpolate_data(datafile_path)
        
        try:
            number_of_profiles,lat,lon,distance = results[0]
            interp_pressure,oxygen_grid,salinity_grid,consv_temperature_grid,density_grid = results[1]
            eps_pressure,eps_grid,eps_consv_temperature_grid,eps_N_squared_grid,eps_density_grid = results[2]
        except TypeError:
            print(cruisename,DATAFILENAME[:-4],"is skipped!")
            continue


        #calculate the viscosity (2 different formula)
        eps_wiki_viscosity_grid = thesis.wiki_viscosity(eps_consv_temperature_grid)/eps_density_grid
        eps_viscosity_grid = thesis.get_viscosity(eps_consv_temperature_grid)

        #calculate the Reynolds bouyancy number defined on eps_pressure
        eps_Reynolds_bouyancy_grid = eps_grid/(eps_viscosity_grid*eps_N_squared_grid)
        eps_wiki_Reynolds_bouyancy_grid = eps_grid/(eps_wiki_viscosity_grid*eps_N_squared_grid)

        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,interp_pressure,density_grid,oxygen_grid,height_above_ground = 10,minimal_density_difference = 0.02)
        bathymetrie,list_of_bathymetrie_indices = results[0]
        BBL,list_of_BBL_indices = results[1]
        BBL_range,list_of_BBL_range_indices = results[2]


        """"
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
            nan_index =  -np.argmax(np.flip(~np.isnan(salinity_grid[i,:]))) #at the moment the index is negative
            nan_index = density_grid[i,:].size + nan_index #now defined as positive index
            
           
            if nan_index == density_grid[i,:].size:
                if not np.isnan(salinity_grid[i,-1]): #if there are no NAN values towards the bottom
                    nan_index = len(interp_pressure)-1
                    
            list_of_bathymetrie_indices[i] = nan_index 
            bathymetrie[i] = interp_pressure[nan_index]
            
            #print(nan_index)
            #print("bathymetrie check\n")
            #print(interp_pressure[nan_index-1],salinity_grid[i,nan_index-1])
            #print(interp_pressure[nan_index],salinity_grid[i,nan_index])
            #print(interp_pressure[nan_index+1],salinity_grid[i,nan_index+1])
            #while interp_coarse_pressure[j]
            #if 
          
            #TODO
            #------------------search for halocline values starts from above:-------------------------------
            


            #------------------search for BBL values starts from below:-------------------------------  
            
            #index of maximal distance bottom plus 15m 
            height_above_ground = 10
            BBL_boundary_index = np.argmax(interp_pressure >= (bathymetrie[i]-height_above_ground))
            assert(interp_pressure[BBL_boundary_index]<bathymetrie[i]) #tests if the point 15m above the ground is really above
            
            #TODO get the index (and from that the pressure) where the density difference is bigger than 0.01
            #BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:density_grid[i,:].size + nan_index])>0.01))
            
            #get the index (and from that the pressure) where the density difference is at maximum (in the lowermost 15 m)
            assert(nan_index>=0)
            BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:nan_index]))) -1 
            
            #check if the maximum is at the edge of the intervall or if the maximum is too small
            minimal_density_difference = 0.02
            if (BBL_index == BBL_boundary_index) or (BBL_index == (BBL_boundary_index+1)) or ((density_grid[i,BBL_index]-density_grid[i,BBL_index-1]) < minimal_density_difference):
                BBL_index = nan_index #equivalent to a BBL thickness of 0
            
            #print(BBL_index,nan_index)
            #print("BBL",interp_pressure[BBL_index])
            #print(bathymetrie[i])
           
            list_of_BBL_indices[i] = BBL_index 
            BBL[i] = interp_pressure[BBL_index]
            
            list_of_BBL_range_indices[i] = BBL_boundary_index 
            BBL_range[i] = interp_pressure[BBL_boundary_index]

        """


        ##########################################################################################################################################################   



        #print(np.max(interp_pressure),np.min(interp_pressure))
        BBL_thickness = bathymetrie - BBL
        
        #searches for the index of the profile with the biggest BBL thickness
        index_of_maximum_thickness = np.argmax(BBL_thickness)
        
        
        first_profile = index_of_maximum_thickness-2
        last_profile = index_of_maximum_thickness + 2


        #handels edge cases, if the choosen profile is to near to the start or end of the transect
        if (last_profile >= (number_of_profiles-1)):
            last_profile = number_of_profiles - 1 
            first_profile = last_profile-4
            index_of_maximum_thickness = first_profile+2

        if (first_profile < 0):
            first_profile = 0
            last_profile = 4 
            index_of_maximum_thickness = 2


        #Plotting
        count = 0
        for transect_index in range(first_profile,last_profile+1): #in total 5 profiles
        
            #print(transect_index,index_of_maximum_thickness,number_of_profiles)
            img2_0 = axarr2[count].plot(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:],interp_pressure[int(list_of_BBL_range_indices[transect_index])-5:],"g", label = "fine grid")
            
            #TODO plot oxygen on a second axis
            #img4_0 = axarr4[count].plot(oxygen_grid[transect_index,:],interp_pressure,"g", label = "fine grid")
            
        
            img2_0b = axarr2[count].plot(np.nanmean(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),BBL[transect_index],"Dr")
            img2_0c = axarr2[count].plot(np.nanmean(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),bathymetrie[transect_index],"Dg")
            img2_0d = axarr2[count].plot(np.nanmean(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),BBL_range[transect_index],"ok")
        
            axarr2[count].set_xlabel("density")
        
        
            img2_0 = axarr2[count+1].plot(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),interp_pressure[int(list_of_BBL_range_indices[transect_index])+1-5:])
        

            img2_1b = axarr2[count+1].plot(np.nanmean(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:])),bathymetrie[transect_index],"Dg")
            img2_1c = axarr2[count+1].plot(np.nanmean(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:])),BBL_range[transect_index],"ok")
            img2_1d = axarr2[count+1].plot(np.nanmean(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:])),BBL[transect_index],"Dr")
        
            count+=2
        
                 
         
        #Plot the data  
        #append the last distance plus the last difference (for plotting all the n profiles we need a distance array of size n+1 
        plotmesh_distance = np.append(distance,2*distance[-1]-distance[-2])
        
        cmap = plt.get_cmap('viridis')
        cmap.set_bad(color = 'lightgrey')
 
        img1_0 = axarr1[0].pcolormesh(plotmesh_distance,interp_pressure,oxygen_grid.T,cmap = cmap)
        img1_1 = axarr1[1].pcolormesh(plotmesh_distance,interp_pressure,density_grid.T,cmap = cmap)

        
        
        #draw the calculated layers in the plot    
        axarr1[0].plot(distance,bathymetrie)
        axarr1[0].plot(distance,BBL)
           
        axarr1[1].plot(distance,bathymetrie)
        axarr1[1].plot(distance,BBL)

        
        axarr1[0].plot(distance[first_profile:last_profile+1],10*np.ones(5),"rD")
        axarr1[1].plot(distance[first_profile:last_profile+1],10*np.ones(5),"rD")
        
            
            
        f1.set_size_inches(18,10.5)
        f2.set_size_inches(18,10.5)

        colorbar(img1_0).set_label(r"density [kg/$m^3$]")
        colorbar(img1_1).set_label(r"oxygen $[%]$")


        for i in range(2):
            if cruisename == "emb169":        
                axarr1[i].set_ylim((0,160))
                assert(np.min(interp_pressure)>0)
                assert(np.max(interp_pressure)<160)
                 
            if cruisename == "emb177":
                axarr1[i].set_ylim((0,135))            
                assert(np.min(interp_pressure)>0)
                assert(np.max(interp_pressure)<135)
                                
            if cruisename == "emb217":
                if DATAFILENAME[:4] == "S106":
                    axarr1[i].set_ylim((0,90))           
                    assert(np.min(interp_pressure)>0)
                    assert(np.max(interp_pressure)<90)                    
                else:
                    axarr1[i].set_ylim((0,160))           
                    assert(np.min(interp_pressure)>0)
                    assert(np.max(interp_pressure)<160)                       



        axarr1[0].invert_yaxis()
        axarr2[0].invert_yaxis()

        f1.suptitle(cruisename+" "+DATAFILENAME[:-4]+" transect")
        f2.suptitle(cruisename+" "+DATAFILENAME[:-4]+ " centered around index "+str(index_of_maximum_thickness))

        

        f1.tight_layout() 
        f2.tight_layout()     
        
        f1.savefig("./BBL_profiles/"+cruisename+"/"+cruisename+"_"+DATAFILENAME[:-4]+"_transect", dpi=450)
        f2.savefig("./BBL_profiles/"+cruisename+"/"+cruisename+"_"+DATAFILENAME[:-4]+"_profiles", dpi=450)

        plt.close(fig = "all")
plt.show()


