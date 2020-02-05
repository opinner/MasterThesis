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

def get_viscosity(T):
    #% T is temperature in degC; vis in m2/s
    #% vis=(1.792747-(.05126103*T)+(0.0005918645*T*T))*1e-6;
    #% Ilker
    return (1.792747-(0.05126103*T)+(0.0005918645*T*T))*1e-6

def wiki_viscosity(T):
    A = 29.39*1e-3
    B = 507.88 
    C = 149.3
    
    return A * np.exp(B/(T-C))
    
#Constants
rho_0 = 1000 #kg/m³
g = 9.81 #m/s² #could be replace by gsw.grav(lat,p)

#contains the MSS Data
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/emb217_mss_data","/home/ole/share-windows/emb177_mss_data/","/home/ole/share-windows/emb169_mss_data/MSS055/matlab/","/home/ole/share-windows/emb169_mss_data/MSS038/matlab/"]

 
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
            

        #define the pictures
        f1, axarr1 = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = True)
        f2, axarr2 = plt.subplots(nrows = 1, ncols = 10, sharey = True)
        #f4, axarr4 = plt.subplots(nrows = 1, ncols = 7, sharey = True)#, sharex = True, 

        print(cruisename+"_"+DATAFILENAME[:-4])
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

        #print(lat)

        pressure = CTD_substructure["P"][0]
        oxygen = CTD_substructure["O2"][0]
        absolute_salinity = CTD_substructure["SA"][0] #is this unit sufficient
        consv_temperature = CTD_substructure["CT"][0] #TODO better use conservative temperature?
        alpha = CTD_substructure["ALPHA"][0]
        beta = CTD_substructure["BETA"][0]

        eps = MIX_substructure["eps"][0]
        eps_pressure = MIX_substructure["P"][0]

        number_of_transects = np.shape(pressure)[-1]

        latitude = []
        longitude = []

        distance = np.zeros(number_of_transects)
        origin = (float(lat[0][0][0]),float(lon[0][0][0])) #lots of brackets to get a number, not an array (workaround)
        for i in range(number_of_transects):
            current_point = (float(lat[i][0][0]),float(lon[i][0][0]))
            latitude.append(float(lat[i][0][0]))
            longitude.append(float(lon[i][0][0]))
            distance[i] = geo.geodesic(origin,current_point).km #Distance in km, change to nautical miles?

        lat = np.asarray(latitude)
        lon = np.asarray(longitude)

        #remove data from a file, with overlapping positional points
        if (datafile_path == "/home/ole/share-windows/emb217_mss_data/TR1-8.mat"):
            lat = np.delete(lat,np.s_[33:47])
            lon = np.delete(lon,np.s_[33:47])
            distance = np.delete(distance,np.s_[33:47])
            
            pressure = np.delete(pressure,np.s_[33:47],axis=0)
            oxygen = np.delete(oxygen,np.s_[33:47],axis=0)
            absolute_salinity =  np.delete(absolute_salinity,np.s_[33:47],axis=0)
            consv_temperature = np.delete(consv_temperature,np.s_[33:47],axis=0)
            alpha = np.delete(alpha,np.s_[33:47],axis=0)
            beta = np.delete(beta,np.s_[33:47],axis=0)
            
            
            eps = np.delete(eps,np.s_[33:47],axis=0)
            eps_pressure = np.delete(eps_pressure,np.s_[33:47],axis=0)
            
            number_of_transects = np.shape(pressure)[-1]
       
       
        if cruisename == "emb169" and  DATAFILENAME[:-4] == "TS118":
            lat = lat[:21]
            lon = lon[:21]
            distance = distance[:21]
            
            pressure = pressure[:21]
            oxygen = oxygen[:21]
            absolute_salinity =  absolute_salinity[:21]
            consv_temperature = consv_temperature[:21]
            alpha = alpha[:21]
            beta = beta[:21]
            
            
            eps = eps[:21]
            eps_pressure = eps_pressure[:21]
            number_of_transects = np.shape(pressure)[-1]
        
        #removes the last data point, that dont seem to belong to the transect
        if cruisename == "emb169" and  DATAFILENAME[:-4] == "TRR109":
            lat = lat[:-1]
            lon = lon[:-1]
            distance = distance[:-1]
            
            pressure = pressure[:-1]
            oxygen = oxygen[:-1]
            absolute_salinity =  absolute_salinity[:-1]
            consv_temperature = consv_temperature[:-1]
            alpha = alpha[:-1]
            beta = beta[:-1]
            
            
            eps = eps[:-1]
            eps_pressure = eps_pressure[:-1]
            number_of_transects = np.shape(pressure)[-1]
    
    
        #test if distance is monotonically increasing
        try:
            assert(np.all(np.diff(distance)>0))
        except AssertionError:  
            print("##########################################")
            print(cruisename,DATAFILENAME[:-4],"is skipped!")
            print("##########################################")
            continue #jump to the next datafile

        #initial values 
        min_pressure = 10
        max_pressure = 60
        max_size = 1000
        min_size = 3000

        #eps profile has a coarser resolution
        min_eps_pressure = 10
        max_eps_pressure = 60
        max_eps_size = 100
        min_eps_size = 400

        #select the start and end point for the
        for i in range(number_of_transects):
            if np.nanmin(pressure[i]) < min_pressure:
                min_pressure = np.nanmin(pressure[i])
            if np.nanmax(pressure[i]) > max_pressure:
                max_pressure = np.nanmax(pressure[i])
            if pressure[i].size > max_size:
                max_size = pressure[i].size       
            if pressure[i].size < min_size:
                min_size = pressure[i].size  

            if np.nanmin(eps_pressure[i]) < min_eps_pressure:
                min_eps_pressure = np.nanmin(eps_pressure[i])
            if np.nanmax(eps_pressure[i]) > max_eps_pressure:
                max_eps_pressure = np.nanmax(eps_pressure[i])
            if eps_pressure[i].size > max_eps_size:
                max_eps_size = eps_pressure[i].size       
            if eps_pressure[i].size < min_eps_size:
                min_eps_size = eps_pressure[i].size  

        #print("High resolution",min_pressure,max_pressure,min_size,max_size)
        #print("Coarse resolution",min_eps_pressure,max_eps_pressure,min_eps_size,max_eps_size)

        #check if that worked correctly
        assert(max_size>= min_size)
        assert(max_eps_size>= min_eps_size)


            
        interp_pressure = np.linspace(min_pressure,max_pressure,min_size)
        #print("max_pressure",max_pressure)


        #creates pressure grid, where every column is equal to interp_pressure
        pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],min_size))

        interp_coarse_pressure = np.linspace(min_eps_pressure,max_eps_pressure,min_eps_size)


        #create grids with distance on x and depth on y-axis
        oxygen_grid = np.zeros((np.shape(pressure)[-1],min_size))
        salinity_grid = np.copy(oxygen_grid)
        consv_temperature_grid = np.copy(oxygen_grid)
        alpha_grid = np.copy(consv_temperature_grid)
        beta_grid = np.copy(salinity_grid)

        #averaged of approx 5 depth bins (???)
        eps_grid = np.zeros((np.shape(pressure)[-1],min_eps_size))


        for i in range(number_of_transects): 
            #print(np.shape(oxygen_grid),np.shape(interp_pressure),np.shape(pressure[i]),np.shape(oxygen[i]),np.shape(eps[i].flatten()))
            oxygen_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),oxygen[i].flatten(), left = np.nan, right = np.nan)
            salinity_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
            consv_temperature_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
            alpha_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),alpha[i].flatten(), left = np.nan, right = np.nan)
            beta_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),beta[i].flatten(), left = np.nan, right = np.nan)
            
            eps_grid[i] = np.interp(interp_coarse_pressure,eps_pressure[i].flatten(),eps[i].flatten(), left = np.nan, right = np.nan)
            
        #calculating the density using teh gsw toolbox
        density_grid = gsw.rho(salinity_grid,consv_temperature_grid,pressure_grid)

        #creates pressure grid, where every column is equal to interp_pressure
        pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],min_size))

        #calculate N^2 with the gsw toolbox 
        BN_freq_squared_grid_gsw, midpoint_pressure_grid = gsw.Nsquared(salinity_grid,consv_temperature_grid,pressure_grid, lat = np.mean(lat), axis = 1)

        #create 1D-array of the midpoint pressure values
        mid_point_pressure = np.mean(midpoint_pressure_grid, axis = 0)
        #print(np.mean(np.std(midpoint_pressure_grid, axis = 0)))
        assert(np.all(np.std(midpoint_pressure_grid, axis = 0) <10**(-10)))

        """
        BN_freq_squared_cleaned_grid = np.copy(BN_freq_squared_grid)
        BN_freq_squared_cleaned_grid_check = np.copy(BN_freq_squared_grid_check)

        #replace all negative values with 0 (is that good?) to compute the frequency
        BN_freq_squared_cleaned_grid[BN_freq_squared_grid < 0] = np.nan #replace all negative values with 0 (is that good?)
        BN_freq_squared_cleaned_grid_check[BN_freq_squared_grid_check < 0]
        #draw the square root
        BN_freq_cleaned_grid = np.sqrt(BN_freq_squared_cleaned_grid)
        BN_freq_cleaned_grid_check = np.sqrt(BN_freq_squared_cleaned_grid_check)
        """

        #interpolate the temperature and to the eps pressure grid
        coarse_consv_temperature_grid = np.zeros(np.shape(eps_grid)) #TODO DO I need the in-situ temperature here?
        coarse_N_squared_grid = np.zeros(np.shape(eps_grid))
        coarse_density_grid = np.copy(coarse_N_squared_grid)

        for i in range(number_of_transects):
             
            #midpoint grid to coarse grid 
            coarse_N_squared_grid[i] = np.interp(interp_coarse_pressure,mid_point_pressure,BN_freq_squared_grid_gsw[i].flatten(), left = np.nan, right = np.nan)
             
            #fine grid to coarse grid
            coarse_consv_temperature_grid[i] = np.interp(interp_coarse_pressure,interp_pressure,consv_temperature_grid[i].flatten(), left = np.nan, right = np.nan)
            coarse_density_grid[i] = np.interp(interp_coarse_pressure,interp_pressure,density_grid[i].flatten(), left = np.nan, right = np.nan)

        coarse_viscosity_T_grid = get_viscosity(coarse_consv_temperature_grid)

        coarse_viscosity_TS_grid = wiki_viscosity(coarse_consv_temperature_grid)/coarse_density_grid


        #all used grids should be defined on the coarse pressure grid
        Reynolds_bouyancy_grid = eps_grid/(coarse_viscosity_T_grid*coarse_N_squared_grid) #from Maffioli 2014
        Reynolds_bouyancy_TS_grid = eps_grid/(coarse_viscosity_TS_grid*coarse_N_squared_grid)

        #search for bottom currents
        bathymetrie = np.zeros(number_of_transects)-99 #fill value (or error value) of -99
        list_of_bathymetrie_indices = np.zeros(number_of_transects)

        halocline = np.zeros(number_of_transects)-99 #fill value (or error value) of -99
        list_of_halocline_indices = np.zeros(number_of_transects)

        BBL = np.zeros(number_of_transects)-99 #fill value (or error value) of -99
        list_of_BBL_indices = np.zeros(number_of_transects)

        BBL_range = np.zeros(number_of_transects)-99 #fill value (or error value) of -99
        list_of_BBL_range_indices = np.zeros(number_of_transects)

        for i in range(number_of_transects):

            #------------------search for bathymetrie values starts from below:-------------------------------

            #returns the pressure of the last nan value in a continuous row starting from high pressure (TODO:is it better to use the last index with data?)
            nan_index =  -np.argmax(np.flip(~np.isnan(salinity_grid[i,:]))) #at the moment the index is negative
            nan_index = density_grid[i,:].size + nan_index #now defined as positive index
            assert(nan_index>=0)
           
            if nan_index == density_grid[i,:].size:
                if not np.isnan(salinity_grid[i,-1]): #if there are no NAN values towards the bottom
                    nan_index = len(interp_pressure)-1
                    
            list_of_bathymetrie_indices[i] = nan_index 
            bathymetrie[i] = interp_pressure[nan_index]
         
            #------------------search for halocline values starts from above:-------------------------------
            #to be added TODO


            #------------------search for BBL values starts from below:-------------------------------  
            
            #index of bottom plus 15m
            height_above_ground = 10 
            BBL_boundary_index = np.argmax(interp_pressure >= (bathymetrie[i]-height_above_ground))
            assert(interp_pressure[BBL_boundary_index]<bathymetrie[i]) #tests if the point 15m above the ground is really above

            
            #get the index (and from that the pressure) where the density difference is at maximum (in the lowermost 15 m)
            BBL_index =  nan_index - np.argmax(np.flip(np.diff(density_grid[i,BBL_boundary_index:nan_index]))) -1 
            
            #check if the maximum is at the edge of the intervall or if the maximum is too small
            minimal_density_difference = 0.02
            if (BBL_index == BBL_boundary_index) or (BBL_index == (BBL_boundary_index+1)) or ((density_grid[i,BBL_index]-density_grid[i,BBL_index-1]) < minimal_density_difference):
                BBL_index = nan_index #equivalent to a BBL thickness of 0
           
            list_of_BBL_indices[i] = BBL_index 
            BBL[i] = interp_pressure[BBL_index]
            
            list_of_BBL_range_indices[i] = BBL_boundary_index 
            BBL_range[i] = interp_pressure[BBL_boundary_index]



        #print(np.max(interp_pressure),np.min(interp_pressure))
        BBL_thickness = bathymetrie - BBL
        
        #searches for the index of the profile with the biggest BBL thickness
        index_of_maximum_thickness = np.argmax(BBL_thickness)
        
        
        first_profile = index_of_maximum_thickness-2
        last_profile = index_of_maximum_thickness + 2

        if (last_profile >= (number_of_transects-1)):
            print("test1")
            last_profile = number_of_transects - 1 
            first_profile = last_profile-4
            index_of_maximum_thickness = first_profile+2

        if (first_profile < 0):
            print("test2")
            first_profile = 0
            last_profile = 4 
            index_of_maximum_thickness = 2


        #Plotting
        count = 0
        for transect_index in range(first_profile,last_profile+1): #in total 5 profiles
        
            print(transect_index,index_of_maximum_thickness,number_of_transects)
            img2_0 = axarr2[count].plot(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:],interp_pressure[int(list_of_BBL_range_indices[transect_index])-5:],"g", label = "fine grid")
            
            #TODO plot oxygen on a second axis
            #img4_0 = axarr4[count].plot(oxygen_grid[transect_index,:],interp_pressure,"g", label = "fine grid")
            
        
            img2_0b = axarr2[count].plot(np.nanmean(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),BBL[transect_index],"Dr")
            img2_0c = axarr2[count].plot(np.nanmean(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),bathymetrie[transect_index],"Dg")
            img2_0d = axarr2[count].plot(np.nanmean(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),BBL_range[transect_index],"ok")
        
            axarr2[count].set_xlabel("density")
        
        
            img2_0 = axarr2[count+1].plot(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:]),mid_point_pressure[int(list_of_BBL_range_indices[transect_index])-5:])
        

            img2_1b = axarr2[count+1].plot(np.nanmean(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:])),bathymetrie[transect_index],"Dg")
            img2_1c = axarr2[count+1].plot(np.nanmean(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:])),BBL_range[transect_index],"ok")
            img2_1d = axarr2[count+1].plot(np.nanmean(np.diff(density_grid[transect_index,int(list_of_BBL_range_indices[transect_index])-5:])),BBL[transect_index],"Dr")
        
            count+=2
        
        
        """
        transect_index = -5
        img4_0 = axarr4[0].plot(oxygen_grid[transect_index,:],interp_pressure,"g", label = "fine grid")
        img4_1 = axarr4[1].plot(salinity_grid[transect_index,:],interp_pressure,"g", label = "fine grid")

        img4_2 = axarr4[2].plot(consv_temperature_grid[transect_index,:],interp_pressure,"g", label = "fine grid")
        img4_2b = axarr4[2].plot(coarse_consv_temperature_grid[transect_index,:],interp_coarse_pressure, label = "coarse grid")

        img4_3 = axarr4[3].plot(BN_freq_squared_grid_gsw[transect_index,:],mid_point_pressure, "r", label = "midpoint grid")
        img4_3b = axarr4[3].plot(coarse_N_squared_grid[transect_index,:],interp_coarse_pressure, label = "coarse grid")

        img4_4 = axarr4[4].plot(coarse_viscosity_T_grid[transect_index,:]*10**6,interp_coarse_pressure,label = "Ilker coarse grid")
        img4_4b = axarr4[4].plot(coarse_viscosity_TS_grid[transect_index,:]*10**6,interp_coarse_pressure,"--",label = "Wikipedia coarse grid")
        img4_5 = axarr4[5].plot(np.log10(eps_grid[transect_index,:]),interp_coarse_pressure, label = "coarse grid")
        img4_6 = axarr4[6].plot(Reynolds_bouyancy_grid[transect_index,:],interp_coarse_pressure,label = "Ilker coarse grid")
        img4_6b = axarr4[6].plot(Reynolds_bouyancy_TS_grid[transect_index,:],interp_coarse_pressure,"--",label = "Wikipedia coarse grid")

        axarr4[0].set_xlabel("oxygen")
        axarr4[0].set_ylabel("pressure [dbar]")
        axarr4[1].set_xlabel("SA")
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
        """
        """
        #interpolate the temperature and  to the eps pressure grid
        coarse_consv_temperature_grid = np.zeros(np.shape(eps_grid)) #TODO DO I need the in-situ temperature here?
        coarse_N_squared_grid = np.zeros(np.shape(eps_grid))

        for i in range(number_of_transects):
            #print(np.shape(mid_point_pressure),np.shape(interp_pressure),np.shape(interp_coarse_pressure))
             
            #midpoint grid to coarse grid 
            coarse_N_squared_grid[i] = np.interp(interp_coarse_pressure,mid_point_pressure,BN_freq_squared_grid_gsw[i].flatten(), left = np.nan, right = np.nan)
             
            #fine grid to coarse grid
            coarse_consv_temperature_grid[i] = np.interp(interp_coarse_pressure,interp_pressure,consv_temperature_grid[i].flatten(), left = np.nan, right = np.nan)
            

        coarse_viscosity_grid = get_viscosity(coarse_consv_temperature_grid)

        #all used grids should be defined on the coarse pressure grid
        Reynolds_bouyancy_grid = eps_grid/(coarse_viscosity_grid*coarse_N_squared_grid) #from Maffioli 2014

        print("Test:",np.nanmedian(eps_grid),np.nanmedian(coarse_viscosity_grid),np.nanmedian(coarse_N_squared_grid))
        print("RB grid:",np.nanmax(Reynolds_bouyancy_grid),np.nanmin(Reynolds_bouyancy_grid),np.nanmedian(Reynolds_bouyancy_grid))

        print(np.max(interp_pressure),np.min(interp_pressure))

        
        """
         
         
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

        
        axarr1[0].plot(distance[first_profile:last_profile+1],-1*np.ones(5),"rD")
        axarr1[1].plot(distance[first_profile:last_profile+1],-1*np.ones(5),"rD")
        
            
            
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


