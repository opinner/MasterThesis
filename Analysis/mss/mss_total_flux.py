import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt
import geopy.distance as geo
import mss_functions as thesis
#from scipy import integrate as sc
import scipy.stats as ss 
import warnings
warnings.filterwarnings('ignore')
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

def flux_look_up_table(flux, transect_bathymetry, depth_array, basin_bathymetry, interval_edge):

    #define 3 





    #because depth_array contains the center point of the bins, the bin edges has to be calculated
    bin_size = abs(depth_array[0]-depth_array[1])
    expanded_array = np.append(depth_array[0]-bin_size,depth_array)
    bin_edges = expanded_array + bin_size/2

    flux = flux * 1000000 #convert flux per m^2 to flux per km^2
    flux_wo_nans = flux[~np.isnan(flux)]
    transect_bathymetry_wo_nans = transect_bathymetry[~np.isnan(flux)]

    total_flux_mean_bin = 0
    total_flux_median_bin = 0
    binned_flux_mean, _bin_edges, _bin_number = ss.binned_statistic(transect_bathymetry_wo_nans,flux_wo_nans,"mean", bin_edges) 
    binned_flux_median, _bin_edges, _bin_number = ss.binned_statistic(transect_bathymetry_wo_nans,flux_wo_nans,"median", bin_edges) 
       
    #############################################################################   
    #TODO: Check if tiles below the lowest transect depth have a contribution. (As they should have)   
    #############################################################################   
    
    binned_flux_mean
        
    for column in range(np.shape(basin_bath)[0]):
        for row in range(np.shape(basin_bath)[1]): 
        
            #check if bathymetry tile has a defined depth
            if not np.isnan(basin_bath[column,row]):
                depth_index = np.nanargmin(np.abs(depth_array - abs(basin_bath[column,row]))) 
                
                #check if the bin is not empty
                if not np.isnan(binned_flux_mean[depth_index]): 
                    total_flux_mean_bin += binned_flux_mean[depth_index]
                    total_flux_median_bin += binned_flux_median[depth_index]
     
    #multiply with the tile size in km² to get rid of the area dependency / to get Gmol/y 
    total_flux_mean_bin *= pixel_area  
    total_flux_median_bin *= pixel_area
                    
    print("Total flux from depth binned mean fluxes",1e-12*365*total_flux_mean_bin,"Gmol/y")
    print("Total flux from depth binned median fluxes",1e-12*365*total_flux_median_bin,"Gmol/y")

def flux_ratio(flux, distance, background_flux, number_of_edge_pixels, number_of_interior_pixels, delta_X,pixel_area,side_length):
    
    flux = flux * 1000000 #convert flux per m^2 to flux per km^2
    background_flux = background_flux *1000000 #convert flux per m^2 to flux per km^2
    
    #set all nans to 0
    flux_without_nans = np.copy(flux)
    flux_without_nans[np.isnan(flux_without_nans)] = 0
    
    """
    unique_distance = []
    unique_flux = []
    averaged = False
    for i in range(len(distance)-1):
    
        if averaged == True:
            averaged = False
            continue
    
        if distance[i] == distance[i+1]:
            #print("test",flux_without_nans[i],flux_without_nans[i+1])
            unique_flux.append(np.nanmean(flux_without_nans[i]+flux_without_nans[i+1]))
            unique_distance.append(distance[i])
            averaged = True
        else:
            unique_flux.append(flux_without_nans[i])
            unique_distance.append(distance[i]) 
            
    assert np.all(np.diff(unique_distance) > 0)
    
    #for index in range(len(flux)):
    #    print(longitude[index],flux[index],np.round(delta_X[index],5),flux[index]*delta_X[index])
    """
        
    interior_flux = background_flux * pixel_area * number_of_interior_pixels
    print("flux through interior",1e-12*365*interior_flux,"Gmol/y")
    
    #edge_flux_own = np.nansum(flux * delta_X) *side_length * number_of_edge_pixels #in units of mmol/d
    edge_flux_trapz = np.trapz(flux_without_nans,distance) *side_length * number_of_edge_pixels #in units of mmol/d
    #edge_flux_simps = sc.simps(unique_flux,unique_distance) * side_length * number_of_edge_pixels #in units of mmol/d
    #print("Riemann:",np.nansum(flux * delta_X)*1e-12*365,"Gmol/y/km")
    
    print("Integral:",np.trapz(flux_without_nans,distance)*1e-12*365,"Gmol/y/km")
    #print("trapz2:",np.trapz(unique_flux,unique_distance)*1e-12*365,"Gmol/y/km")
    #print("cumtrapz:",sc.cumtrapz(unique_flux,unique_distance)[-1]*1e-12*365,"Gmol/y/km")
    #print("simps:",sc.simps(unique_flux,unique_distance)*1e-12*365,"Gmol/y/km")
    print("flux through edge",1e-12*365*edge_flux_trapz,"Gmol/y") #,1e-12*365*edge_flux_own,1e-12*365*edge_flux_simps)

    
    print("proportion of edge fluxes",np.round(100*edge_flux_trapz/(edge_flux_trapz + interior_flux),2), "%")
    
    print("Total flux from integral",1e-12*365*(edge_flux_trapz + interior_flux),"Gmol/y")
    
    
    """
    f,a = plt.subplots(1)
    a.plot(distance, flux, label = "original")
    a.plot(distance, flux_without_nans, ".")
    a.plot(unique_distance,unique_flux, ".")
    a.plot(distance, flux_without_nans, label = "no nans")
    a.plot(unique_distance,unique_flux, "k--", label = "unique")
    a.legend()
    f.set_size_inches(18.5,10)
    f.tight_layout()
    plt.show()
    """
    
    return 
        
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################        
###############################################################################################################################

out,axis = plt.subplots(ncols = 3, sharex = True, sharey = True)
out2,axis2 = plt.subplots(ncols = 3, sharex = True, sharey = True)
out3,axis3 = plt.subplots(1)

#TODO aktualisiere die Haloklinenwerte

#loop over the three cruises and their mean halocline depth
for cruise_index,cruise_name,set_depth,interval_edge in zip([0,1,2],["emb169","emb177","emb217"],[-73.23,-58.16,-72.15],[-79.6,-60.39,-77.58]):


    #load the bathymetry data        
    dset = nc.Dataset("/home/ole/windows/all_data/bathymetry_data/iowtopo2_rev03.nc")
    #print(dset)

    z_topo = dset.variables["Z_WATER"]
    XT_I = dset.variables["XT_I"]
    YT_J = dset.variables["YT_J"]

    #convert to a grid
    lon, lat = np.meshgrid(XT_I,YT_J)
    bathymetry = np.asarray(z_topo)

    #replace the fill_values with NaN
    bathymetry[bathymetry == -1.000e+34] = np.nan

    #set the desired depth of the <<halocline>> 
    

    #basin_shape = np.copy(bathymetry)
    basin_shape_list = bathymetry[bathymetry < set_depth]

    basin_lat_list = lat[bathymetry < set_depth]
    basin_lon_list = lon[bathymetry < set_depth]
    #print(np.shape(basin_lat_list),np.shape(basin_lon_list),np.shape(basin_shape_list))


    #only inside the rectangle defined by this values
    left = np.argmin(np.abs(np.mean(lon,axis=0)-18.5))
    right = np.argmin(np.abs(np.mean(lon,axis=0)-21))
    bottom = np.argmin(np.abs(np.mean(lat,axis=1)-56.5))
    top = np.argmin(np.abs(np.mean(lat,axis=1)-57.75))

    #left = np.argmin(np.abs(np.mean(lon,axis=0)-16.8))
    #right = np.argmin(np.abs(np.mean(lon,axis=0)-25))
    #bottom = np.argmin(np.abs(np.mean(lat,axis=1)-54.0))
    #top = np.argmin(np.abs(np.mean(lat,axis=1)-59.8))
    
    #print(lon[0,left],lon[0,right],lat[bottom,0],lat[top,0])

    #print("lat",np.min(np.mean(lat,axis=1)),np.max(np.mean(lat,axis=1)),np.shape(np.mean(lat,axis=1)))
    #print("lon",np.min(np.mean(lon,axis=0)),np.max(np.mean(lon,axis=0)),np.shape(np.mean(lon,axis=0)))

    #print(left,right,bottom,top)

    #create a subset, a smaller matrix only with the area of interest
    basin_bath = bathymetry[bottom:top,left:right]
    basin_lat = lat[bottom:top,left:right]
    basin_lon = lon[bottom:top,left:right]
    basin_mask = (basin_bath < set_depth) #True for everything deeper that the set_depth

    #calculate average pixel size:
    lat_array = np.mean(basin_lat,axis=1)
    mean_lat = np.mean(lat_array)
    mean_lat_diff = np.mean(np.diff(lat_array))
    lon_array = np.mean(basin_lon,axis=0)
    mean_lon = np.mean(lon_array)
    mean_lon_diff = np.mean(np.diff(lon_array))

    x_distance = geo.distance((mean_lat,mean_lon),(mean_lat,mean_lon+mean_lon_diff)).km
    y_distance = geo.distance((mean_lat,mean_lon),(mean_lat+mean_lat_diff,mean_lon)).km

    pixel_area = x_distance*y_distance

    #print("lat",mean_lat,mean_lat_diff)
    #print("lon",mean_lon,mean_lon_diff)
    #print("result",x_distance,y_distance,pixel_area)

    #calculating the edge
    #--------------------------------------
    horizontal_edges = 0 #Number of horizontal edges
    vertical_edges = 0 #Number of vertical edges
                    
    edge_array = np.zeros(np.shape(basin_bath))
    for column in range(np.shape(basin_bath)[0]):
        for row in range(np.shape(basin_bath)[1]):  
        
            #check if the pixel falls inside the basin boundary
            #if not, it can be skipped
            if basin_mask[column,row] != True:
                continue
        
            #try,except construct to account for the edges of the array
            try:
            
                #number how many sides of the pixel are on the edge
                count = 0

                #check if neighboring pixels are not in the basin with a defined height
                if basin_mask[column+1,row] != True:
                    count+=1
                    vertical_edges += 1
                if basin_mask[column-1,row] != True: 
                    count+=1
                    vertical_edges += 1
                if basin_mask[column,row+1] != True:
                    count+=1
                    horizontal_edges += 1
                if basin_mask[column,row-1] != True:    
                    count+=1
                    horizontal_edges += 1
                    
                #check also diagonal neigbors?
                #if basin_mask[column+1,row+1] != True: count+=1
                #if basin_mask[column-1,row-1] != True: count+=1
                #if basin_mask[column+1,row-1] != True: count+=1
                #if basin_mask[column-1,row+1] != True: count+=1
                
              
                edge_array[column,row] = count
                
            except IndexError:
                continue
            
    edge_mask = (edge_array != 0) #True if the point is on the edge

    #print(np.shape(lat),np.shape(lon),np.shape(bathymetry))
    #print(np.shape(basin_lat),np.shape(basin_lon),np.shape(basin_bath))
    #print("basin lat",np.min(np.mean(basin_lat,axis=1)),np.max(np.mean(basin_lat,axis=1)),np.shape(np.mean(basin_lat,axis=1)))
    #print("basin lon",np.min(np.mean(basin_lon,axis=0)),np.max(np.mean(basin_lon,axis=0)),np.shape(np.mean(basin_lon,axis=0)))

    assert(np.all(np.shape(basin_mask) == np.shape(basin_lon)))

    #left = np.argmin(np.abs(np.mean(lon,axis=0)-16.8))
    #right = np.argmin(np.abs(np.mean(lon,axis=0)-25))
    #bottom = np.argmin(np.abs(np.mean(lat,axis=1)-54.0))
    #top = np.argmin(np.abs(np.mean(lat,axis=1)-59.8))
    
    """
    #set evertything outside this box to 0
    basin_mask[basin_lat < 56.5] = 0
    edge_mask[basin_lat < 56.5] = 0
    basin_mask[basin_lat > 57.75] = 0
    edge_mask[basin_lat >57.75] = 0
    basin_mask[basin_lon < 18.5] = 0
    edge_mask[basin_lon < 18.5] = 0
    basin_mask[basin_lon > 21] = 0
    edge_mask[basin_lon >21] = 0
    """
    
    pixels_edge =  len(basin_lat[edge_mask])
    total_pixels_basin = len(basin_lat[basin_mask])
    pixels_interior = total_pixels_basin-pixels_edge

    weighted_mean_pixel_side_length = (horizontal_edges * x_distance + vertical_edges * y_distance)/(horizontal_edges + vertical_edges)
    #print("weighted_mean_pixel_side_length",weighted_mean_pixel_side_length)        



    ratio = pixels_edge/(total_pixels_basin-pixels_edge)
    #print("\nbasin:",total_pixels_basin,"\tedge",pixels_edge,"\tratio",ratio,"\tarea",total_pixels_basin*pixel_area,"\n")


    #longitude,distance,transect_bathymetry,raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih = np.loadtxt("./data/"+cruise_name+"_bin_flux_results.txt", unpack=True)
    longitude,distance,transect_bathymetry,interval_ground_distance,raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih = np.loadtxt("./data/"+cruise_name+"_iso_flux_results.txt", unpack=True)



    ##########################################################################################
    #the colorscheme ['#d95f02','#7570b3','#1b9e77'] stems from colorbrewer (colorbrewer2.org) to be friendly to color blindness and colored printing
    for color,label_name,cruise in zip(['#d95f02','#7570b3','#1b9e77'],["summer cruise emb217","winter cruise emb177","autumn cruise emb169"],["emb217","emb177","emb169"]):
        if cruise_name == cruise:
            break
            
            
    sorted_transect_bathymetry = sorted(transect_bathymetry)
    sorted_raw_Shih = [x for _,x in sorted(zip(transect_bathymetry,raw_Shih))]
    sorted_mean_Shih = [x for _,x in sorted(zip(transect_bathymetry,rolling_mean_Shih))]
      
    axis[cruise_index].plot(sorted_raw_Shih,sorted_transect_bathymetry,"-", label = r"$\langle$Shih flux$\rangle_{\mathrm{z}}$ ")
    axis[cruise_index].plot(sorted_mean_Shih,sorted_transect_bathymetry,"-", label = r"$\langle$Shih flux$\rangle_{\mathrm{z,lon}}$ ")
    axis[cruise_index].set_title(cruise_name)
    axis[cruise_index].legend(loc = "lower center")
    axis[cruise_index].set_xlabel("Shih Oxygen flux")
    
    #remove nans
    mask = np.isnan(interval_ground_distance)
    plot_interval_ground_distance = interval_ground_distance[~mask]
    plot_raw_Shih = raw_Shih[~mask]
    plot_rolling_mean_Shih = rolling_mean_Shih[~mask]
    
    sorted_distance_from_ground = sorted(plot_interval_ground_distance)
    assert np.all(np.diff(sorted_distance_from_ground)>= 0)
    sorted_raw_Shih = [x for _,x in sorted(zip(plot_interval_ground_distance,plot_raw_Shih))]
    sorted_mean_Shih = [x for _,x in sorted(zip(plot_interval_ground_distance,plot_rolling_mean_Shih))]
    
    axis2[cruise_index].plot(sorted_raw_Shih,sorted_distance_from_ground,"-", label = r"$\langle$Shih flux$\rangle_{\mathrm{z}}$ ")
    axis2[cruise_index].plot(sorted_mean_Shih,sorted_distance_from_ground,"-", label = r"$\langle$Shih flux$\rangle_{\mathrm{z,lon}}$ ")
    axis2[cruise_index].set_title(cruise_name)
    axis2[cruise_index].legend(loc = "lower center")
    axis2[cruise_index].set_xlabel("Shih Oxygen flux")    

    axis3.plot(longitude,raw_Shih,":", color = color, alpha = 0.7, label = cruise_name +r" $\langle$Shih flux$\rangle_{\mathrm{z}}$ ")
    axis3.plot(longitude,rolling_mean_Shih,"-", color = color, label = cruise_name+r" $\langle$Shih flux$\rangle_{\mathrm{z,lon}}$ ")
    axis3.set_title(cruise_name)
    axis3.legend(loc = "lower center")
    axis3.set_xlabel("Shih Oxygen flux")    
        
    ##########################################################################################

    
    depth_array = np.arange(0,400,5) #sorted(list(set(transect_bathymetry)))


    diff_distance = thesis.central_differences(distance)
    #diff_distance = np.diff(distance)
    #assert np.all(diff_distance > 0)
    #print(np.round(diff_distance,3))
    #print(diff_distance >= 0)
    
    
    #integral_axis[0,cruise_index].bar(distance,raw_Osborn,alpha = 0.5, width = diff_distance, align = "edge")
    #integral_axis[0,cruise_index].plot(distance,raw_Osborn,"k-")
    #integral_axis[1,cruise_index].bar(distance,rolling_mean_Shih,alpha = 0.5, width = diff_distance, align = "edge")
    #integral_axis[1,cruise_index].plot(distance,rolling_mean_Shih,"k-")


    print("-------------------------------------------------------------------------------")
    print(cruise_name)
    print("For a basin with a halocline at depth",set_depth,"m with a total area of circa",int(total_pixels_basin*pixel_area),"km^2 consisting of",total_pixels_basin,"data points\n")
    for name,flux in zip(["raw_Osborn","rolling_mean_Osborn","raw_Shih","rolling_mean_Shih"],[raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih]):
        print(name)
        print("mean flux:", np.nanmean(flux),"mmol/m^2/d")
        print("max flux:",np.nanmin(flux),"mmol/m^2/d")
        flux_ratio(flux, distance, -1,pixels_edge, pixels_interior, diff_distance,pixel_area,weighted_mean_pixel_side_length)
        
        flux_look_up_table(flux, transect_bathymetry, depth_array, basin_bath, interval_edge)
        print("\n")



axis[0].set_ylabel("sea floor depth")
axis[0].invert_yaxis()
out.set_size_inches(6.2012,6.2012/1.618)
out.subplots_adjust(top=0.904,bottom=0.161,left=0.122,right=0.962,hspace=0.2,wspace=0.182)

axis2[0].set_ylabel("distance to ground")
axis2[0].invert_yaxis()
out2.set_size_inches(6.2012,6.2012/1.618)
out2.subplots_adjust(top=0.904,bottom=0.161,left=0.122,right=0.962,hspace=0.2,wspace=0.182)


#out.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/flux_bathymetry_relation", dpi = 600)
plt.show()



