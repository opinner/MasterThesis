#-----------------------------------------------------------#
#Plots a Temperature and Dissipation Transect for every cruise
#with added isopycnals and the postion of the chain- and eddy-measurment
#
#Plots an oxygen transect plus a zoom of emb217
#with isopycnals and mooring positions
#-----------------------------------------------------------#


import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import geopy.distance as geo
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gsw 


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)

from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from mpl_toolkits.axes_grid1.inset_locator import (
    BboxPatch, BboxConnector, BboxConnectorPatch)


def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = {
            **prop_lines,
            "alpha": prop_lines.get("alpha", 1) * 0.2,
        }

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    # loc1a=3, loc2a=2, loc1b=4, loc2b=1
    p = BboxConnectorPatch(bbox1, bbox2,loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,**prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    Connect *ax1* and *ax2*. The *xmin*-to-*xmax* range in both axes will
    be marked.

    Parameters
    ----------
    ax1
        The main axes.
    ax2
        The zoomed axes.
    xmin, xmax
        The limits of the colored area in both plot axes.
    **kwargs
        Arguments passed to the patch constructor.
    """

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches = {**kwargs, "ec": "none", "alpha": 0.17}
    invisible_prop_patches = {**kwargs, "ec": "none", "alpha": 0.0}

    c1, c2, bbox_patch1, bbox_patch2, p = connect_bbox(mybbox1, mybbox2,loc1a=3, loc2a=2, loc1b=4, loc2b=1,prop_lines=kwargs, prop_patches=prop_patches)

    invisible_c1, invisible_c2, invisible_bbox_patch1, invisible_bbox_patch2, invisible_p = connect_bbox(mybbox1, mybbox2,loc1a=3, loc2a=2, loc1b=4, loc2b=1,prop_lines=kwargs, prop_patches=invisible_prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(invisible_bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p
    

#-------------------------------------------------------------------------------------------------------------------------------#
    
FILENAMES = ["/home/ole/share-windows/emb217_mss_data/TR1-4.mat","/home/ole/share-windows/emb177_mss_data/TS1_8.mat","/home/ole/share-windows/emb169_mss_data/MSS055/matlab/TS11.mat"]


#define the pictures
f1, axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = True, sharey = True)
f2, axarr2 = plt.subplots(nrows = 2)#, sharex = True, sharey = True)

cmap_hot = plt.get_cmap('hot_r')
cmap_hot.set_bad(color = 'lightgrey')

cmap_RdBu = plt.get_cmap('RdBu_r')
cmap_RdBu.set_bad(color = 'lightgrey')

for FILENAME in FILENAMES:

    splitted_filename = FILENAME.split("/")
    cruisename = splitted_filename[4][0:6]
    print("cruisename",cruisename)  
    if cruisename == "emb169":
        y_position = 0
        axarr1[y_position,0].set_title("Temperature from EMB169 (October 2017)",fontweight="bold")
        axarr1[y_position,1].set_title("Dissipation from EMB169 (October 2017)" ,fontweight="bold")

        axarr1[y_position,0].set_facecolor('lightgrey')
        axarr1[y_position,1].set_facecolor('lightgrey')

    elif cruisename == "emb177":
        y_position = 1
        axarr1[y_position,0].set_title("Temperature from EMB177 (March 2018)",fontweight="bold")
        axarr1[y_position,1].set_title("Dissipation from EMB177 (March 2018)",fontweight="bold")
        
        axarr1[y_position,0].set_facecolor('lightgrey')
        axarr1[y_position,1].set_facecolor('lightgrey')
                    
    elif cruisename == "emb217":
        y_position = 2
        axarr1[y_position,0].set_title("Temperature from EMB217 (July 2019)",fontweight="bold")
        axarr1[y_position,1].set_title("Dissipation from EMB217 (July 2019)",fontweight="bold")        

        axarr1[y_position,0].set_facecolor('lightgrey')
        axarr1[y_position,1].set_facecolor('lightgrey')
        
        axarr2[0].set_facecolor('lightgrey')
        axarr2[1].set_facecolor('lightgrey')
                
        
    print("Filename:",sio.whosmat(FILENAME))

    data = sio.loadmat(FILENAME)


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
    absolute_salinity = CTD_substructure["SA"][0] #is this unit sufficient
    consv_temperature = CTD_substructure["CT"][0] #TODO better use conservative temperature?
    oxygen = CTD_substructure["O2"][0]
    
    eps = MIX_substructure["eps"][0]
    eps_pressure = MIX_substructure["P"][0]

    number_of_profiles = np.shape(pressure)[-1]

    latitude = []
    longitude = []

    distance = np.zeros(number_of_profiles)
    origin = (float(lat[0][0][0]),float(lon[0][0][0])) #lots of brackets to get a number, not an array (workaround)
    for i in range(number_of_profiles):
        current_point = (float(lat[i][0][0]),float(lon[i][0][0]))
        latitude.append(float(lat[i][0][0]))
        longitude.append(float(lon[i][0][0]))
        distance[i] = geo.geodesic(origin,current_point).km #Distance in km, change to nautical miles?

    lat = np.asarray(latitude)
    lon = np.asarray(longitude)


    print("lat",max(lat),min(lat))
    print("lon",max(lon),min(lon))

    #remove data from a file, with overlapping positional points
    if (FILENAME == "/home/ole/share-windows/emb217_mss_data/TR1-8.mat"):
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
        
        number_of_profiles = np.shape(pressure)[-1]

    #remove excess data that already belongs to TS119
    if cruisename == "emb169" and FILENAME[:-4] == "TS118":
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
        number_of_profiles = np.shape(pressure)[-1]

    #removes the last data point, that dont seem to belong to the transect
    if cruisename == "emb169" and FILENAME[:-4] == "TRR109":
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
        number_of_profiles = np.shape(pressure)[-1]
                
                            
    if cruisename == "emb169" and FILENAME[:-4] == "TRR109":
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
        number_of_profiles = np.shape(pressure)[-1]

                
    #test if distance is monotonically increasing
    assert(np.all(np.diff(distance)>0))
    #plt.plot(lon)
    #plt.show()
    #print(np.diff(lon))
    assert(all(item >= 0 for item in lon) or all(item < 0 for item in lon))

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
    for i in range(number_of_profiles):

        
        assert(np.all(eps_pressure[i] == eps_pressure[0]))
      

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

    print("High resolution",min_pressure,max_pressure,min_size,max_size)
    print("Coarse resolution",min_eps_pressure,max_eps_pressure,min_eps_size,max_eps_size)

    #check if that worked correctly
    assert(max_size>= min_size)
    assert(max_eps_size>= min_eps_size)

    test_pressure = eps_pressure[0] + np.diff(eps_pressure[0])/2
        
    interp_pressure = np.linspace(min_pressure,max_pressure,min_size)

    #creates pressure grid, where every column is equal to interp_pressure
    pressure_grid = np.reshape(interp_pressure,(1,-1))*np.ones((np.shape(pressure)[-1],min_size))

    interp_coarse_pressure = np.linspace(min_eps_pressure,max_eps_pressure,min_eps_size)


    #create grids with distance on x and depth on y-axis
    salinity_grid = np.zeros((np.shape(pressure)[-1],min_size))
    consv_temperature_grid = np.copy(salinity_grid)
    oxygen_grid = np.copy(salinity_grid)

    #averaged of approx 5 depth bins (???)
    eps_grid = np.zeros((np.shape(pressure)[-1],min_eps_size))

    test_pressure = eps_pressure[0].flatten() + np.mean(np.diff(eps_pressure[0].flatten()))/2
    test_pressure = np.append(eps_pressure[0].flatten()[0]-np.mean(np.diff(eps_pressure[0].flatten()))/2, test_pressure)


    test_pressure_grid = np.reshape(test_pressure,(1,-1))*np.ones((np.shape(eps_pressure)[-1],test_pressure.size))
    test_salinity_grid = np.ones((np.shape(test_pressure_grid)))
    test_consv_temperature_grid = np.ones((np.shape(test_pressure_grid)))
    test_eps_grid = np.copy(eps_grid)

    eps_salinity_grid = np.ones((np.shape(eps_grid)))
    eps_consv_temperature_grid = np.ones(np.shape(eps_grid))

    #interpolation from measurement data to different grids
    for i in range(number_of_profiles): 
        
        #interpolation to a common fine grid
        salinity_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),absolute_salinity[i].flatten(), left = np.nan, right = np.nan)
        consv_temperature_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),consv_temperature[i].flatten(), left = np.nan, right = np.nan)
        oxygen_grid[i] = np.interp(interp_pressure,pressure[i].flatten(),oxygen[i].flatten(), left = np.nan, right = np.nan)  
         
        eps_grid[i] = eps[i].flatten() 
            
    density_grid = gsw.rho(salinity_grid,consv_temperature_grid,pressure_grid)

    #---------------------------------
    #number_of_isopycnals = 10
    #density_steps = np.linspace(np.nanmin(density_grid),np.nanmax(density_grid),number_of_isopycnals, endpoint = False)
    density_steps = np.arange(np.nanmin(density_grid),np.nanmax(density_grid),0.5)

    for density_step in density_steps:
        isopycnal = np.zeros(number_of_profiles)


        for i in range(number_of_profiles):
            iso_index = np.argmax(density_grid[i,:] >= density_step)
            
            isopycnal[i] = interp_pressure[iso_index]
            
            if iso_index == (density_grid[i,:].size -1):
                isopycnal[i] = np.nan
            if iso_index == 0:            
                isopycnal[i] = np.nan


        axarr1[y_position,0].plot(lon,isopycnal,"k")
        axarr1[y_position,1].plot(lon,isopycnal,"k")
        
        if cruisename == "emb217":
            axarr2[0].plot(lon,isopycnal,"k")  
            
    #-----------------------------------


    #------------------search for bathymetrie values starts from below:-------------------------------
    bathymetrie = np.zeros(number_of_profiles)-99 #fill value (or error value) of -99
    list_of_bathymetrie_indices = np.zeros(number_of_profiles)


    for i in range(number_of_profiles):

        #returns the pressure of the last nan value in a continuous row starting from high pressure (TODO:is it better to use the last index with data?)
        nan_index =  -np.argmax(np.flip(~np.isnan(salinity_grid[i,:]))) #at the moment the index is negative
        nan_index = density_grid[i,:].size + nan_index #now defined as positive index
        
        if nan_index == density_grid[i,:].size:
            if not np.isnan(salinity_grid[i,-1]): #if there are no NAN values towards the bottom
                nan_index = len(interp_pressure)-1
                
        list_of_bathymetrie_indices[i] = nan_index 
        bathymetrie[i] = interp_pressure[nan_index]
 
    
    
    #append the last distance plus the last difference (for plotting all the n profiles we need a distance array of size n+1 
    plotmesh_longitude = np.append(lon,2*lon[-1]-lon[-2])

      
    #print("TEST:",np.any(eps_grid==0))  
    #eps_grid[eps_grid==0] = 1e-10
       
    #Plot the data   
    img1_0 = axarr1[y_position,0].pcolormesh(plotmesh_longitude,interp_pressure,consv_temperature_grid.T, cmap = cmap_RdBu, vmin = 5, vmax = 15)
    img1_1 = axarr1[y_position,1].pcolormesh(plotmesh_longitude,interp_coarse_pressure,np.log10(eps_grid.T), vmax = -6.5, vmin = -9, cmap = cmap_hot)
    
    depth_at_tc_flach = bathymetrie[np.argmin(np.absolute(lon-20.6150))]
    depth_at_tc_tief = bathymetrie[np.argmin(np.absolute(lon-20.6000))]
    tc_flach = axarr1[y_position,0].vlines(20.6150,1,depth_at_tc_flach,color = "r", label = "TC-chains")
    tc_flach = axarr1[y_position,1].vlines(20.6150,1,depth_at_tc_flach,color = "r", label = "TC-chains")
    tc_tief = axarr1[y_position,0].vlines(20.6000,1,depth_at_tc_tief,color = "r")
    tc_tief = axarr1[y_position,1].vlines(20.6000,1,depth_at_tc_tief,color = "r")
    
    depth_at_eddy1 = bathymetrie[np.argmin(np.absolute(lon-20.6018))]
    depth_at_eddy2 = bathymetrie[np.argmin(np.absolute(lon-20.6151))]
    eddy1 = axarr1[y_position,0].vlines(20.6018,1,depth_at_eddy1, label = "Eddy")
    eddy1 = axarr1[y_position,1].vlines(20.6018,1,depth_at_eddy1, label = "Eddy", color = "g")
    eddy2 = axarr1[y_position,0].vlines(20.6155,1,depth_at_eddy2)
    eddy2 = axarr1[y_position,1].vlines(20.6155,1,depth_at_eddy2, color = "g")
    
    #axarr1[y_position,0].plot(lon,bathymetrie)
    #axarr1[y_position,1].plot(lon,bathymetrie)
     
    axarr1[y_position,0].legend(loc = "lower right")
    axarr1[y_position,1].legend(loc = "lower right")
        
    #axarr1[y_position,0].annotate("Temperature", xy=(-12, -12), xycoords='axes points',size=14, ha='left', va='bottom',bbox=dict(boxstyle='round', fc='w'))

    #axarr1[y_position,1].annotate("Dissipation rate", xy=(-12, -12), xycoords='axes points',size=14, ha='left', va='bottom',bbox=dict(boxstyle='round', fc='w'))


    colorbar(img1_0).set_label("conservative temperature [C]")
    colorbar(img1_1).set_label(r"log10($\epsilon$) [W/kg]") 
    
    if cruisename == "emb217":
        img2_0 = axarr2[0].pcolormesh(plotmesh_longitude,interp_pressure,oxygen_grid.T, cmap = cmap_RdBu)#, vmin = 5, vmax = 15)
        img2_1 = axarr2[1].pcolormesh(plotmesh_longitude,interp_pressure,oxygen_grid.T, cmap = cmap_RdBu)#, vmin = 5, vmax = 15)

        tc_flach = axarr2[0].vlines(20.6150,1,depth_at_tc_flach,color = "r", label = "TC-chains")
        tc_flach = axarr2[1].vlines(20.6150,1,depth_at_tc_flach,color = "r", label = "TC-chains")
        tc_tief = axarr2[0].vlines(20.6000,1,depth_at_tc_tief,color = "r")
        tc_tief = axarr2[1].vlines(20.6000,1,depth_at_tc_tief,color = "r")
        
        eddy1 = axarr2[0].vlines(20.6018,1,depth_at_eddy1, label = "Eddy", color = "k")
        eddy1 = axarr2[1].vlines(20.6018,1,depth_at_eddy1, label = "Eddy", color = "k")
        eddy2 = axarr2[0].vlines(20.6155,1,depth_at_eddy2, color = "k")
        eddy2 = axarr2[1].vlines(20.6155,1,depth_at_eddy2, color = "k")
         
        axarr2[0].legend(loc = "lower right")
        axarr2[1].legend(loc = "lower right")
    
        colorbar(img2_0).set_label("Oxygen [%]")
        colorbar(img2_1).set_label("Oxygen [%]")

        zoom_effect01(axarr2[0],axarr2[1],20.59,20.63, linestyle = "--")

        f2.suptitle("Oxygen saturation from EMB217 (July 2019)",fontweight="bold")
        #axarr2[1].set_title("Dissipation from EMB169",fontweight="bold")


axarr1[2,0].set_xlabel("Longitude")
axarr1[2,1].set_xlabel("Longitude")

axarr1[0,0].set_ylabel("pressure [dbar]")
axarr1[1,0].set_ylabel("pressure [dbar]")
axarr1[2,0].set_ylabel("pressure [dbar]")
    
axarr1[0,0].set_ylim((1,150))
axarr1[1,0].set_ylim((1,150))
axarr1[2,0].set_ylim((1,150))
    
f1.set_size_inches(1.618*7.2,7.2)
f1.set_size_inches(18,10.5)

#inverts all axis form figure 1, because of shared xy
axarr1[2,1].invert_yaxis()

f1.tight_layout() 

axarr2[0].set_xlabel("Longitude")
axarr2[1].set_xlabel("Longitude")

axarr2[0].set_ylabel("pressure [dbar]")
axarr2[1].set_ylabel("pressure [dbar]")


axarr2[1].set_xlim((20.59,20.63))
axarr2[1].set_ylim((1,120))

axarr2[0].invert_yaxis()
axarr2[1].invert_yaxis()

f2.subplots_adjust(top=0.95,bottom=0.056,left=0.039,right=0.961,hspace=0.2,wspace=0.2)

f2.set_size_inches(1.618*7.2,7.2)
f2.set_size_inches(18,10.5)

f2.tight_layout()
 
plt.show()


