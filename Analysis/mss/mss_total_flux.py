import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt
import geopy.distance as geo
import mss_functions as thesis
from scipy import integrate as sc
import warnings
warnings.filterwarnings('ignore')

#integral,integral_axis = plt.subplots(ncols = 3, nrows = 2,sharex = True, sharey = True)

for cruise_index,cruisename,set_depth in zip([0,1,2],["emb169","emb177","emb217"],[-71.07,-58.52,-70.35]):


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

    pixels_edge =  len(basin_lat[edge_mask])
    total_pixels_basin = len(basin_lat[basin_mask])
    pixels_interior = total_pixels_basin-pixels_edge

    weighted_mean_pixel_side_length = (horizontal_edges * x_distance + vertical_edges * y_distance)/(horizontal_edges + vertical_edges)
    #print("weighted_mean_pixel_side_length",weighted_mean_pixel_side_length)        



    ratio = pixels_edge/(total_pixels_basin-pixels_edge)
    #print("\nbasin:",total_pixels_basin,"\tedge",pixels_edge,"\tratio",ratio,"\tarea",total_pixels_basin*pixel_area,"\n")


    longitude,distance,raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih = np.loadtxt("./"+cruisename+"_flux_results.txt", unpack=True)

    diff_distance = thesis.central_differences(distance)
    #diff_distance = np.diff(distance)
    #assert np.all(diff_distance > 0)
    #print(np.round(diff_distance,3))
    #print(diff_distance >= 0)
    
    
    hist_diff = np.copy(diff_distance)
    hist_diff[0] = 0
    hist_diff[-1] = hist_diff[-2]
    #print(hist_diff > 0)    
    
    #print(distance+hist_diff)
    #print(np.diff(distance) > 0)
    
    #integral_axis[0,cruise_index].bar(distance,raw_Osborn,alpha = 0.5, width = diff_distance, align = "edge")
    #integral_axis[0,cruise_index].plot(distance,raw_Osborn,"k-")
    #integral_axis[1,cruise_index].bar(distance,rolling_mean_Shih,alpha = 0.5, width = diff_distance, align = "edge")
    #integral_axis[1,cruise_index].plot(distance,rolling_mean_Shih,"k-")

    def flux_ratio(flux, distance, background_flux, number_of_edge_pixels, number_of_interior_pixels, delta_X,pixel_area,side_length):
        
        flux = flux * 1000000 #convert flux per m^2 to flux per km^2
        background_flux = background_flux *1000000 #convert flux per m^2 to flux per km^2
        
        flux_without_nans = np.copy(flux)
        flux_without_nans[np.isnan(flux_without_nans)] = 0
        
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
        

        
        edge_flux_own = np.nansum(flux * delta_X) *side_length * number_of_edge_pixels #in units of mmol/d
        edge_flux_trapz = np.trapz(flux_without_nans,distance) *side_length * number_of_edge_pixels #in units of mmol/d
        edge_flux_simps = sc.simps(unique_flux,unique_distance) * side_length * number_of_edge_pixels #in units of mmol/d
        print("integral:",np.nansum(flux * delta_X)*1e-12*365,"Gmol/y/km")
        print("trapz:",np.trapz(flux_without_nans,distance)*1e-12*365,"Gmol/y/km")
        print("trapz2:",np.trapz(unique_flux,unique_distance)*1e-12*365,"Gmol/y/km")
        print("cumtrapz:",sc.cumtrapz(unique_flux,unique_distance)[-1]*1e-12*365,"Gmol/y/km")
        print("simps:",sc.simps(unique_flux,unique_distance)*1e-12*365,"Gmol/y/km")
        print("flux through edge",1e-12*365*edge_flux_trapz,"Gmol/y",1e-12*365*edge_flux_own)#,1e-12*365*edge_flux_simps)
        interior_flux = background_flux * pixel_area * number_of_interior_pixels
        print("flux through interior",1e-12*365*interior_flux,"Gmol/y")
        print("proportion of edge fluxes",np.round(100*edge_flux_trapz/(edge_flux_trapz + interior_flux),2), "%")
        
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
        
        return 

    print("-------------------------------------------------------------------------------")
    print(cruisename)
    print("For a basin with a halocline at depth",set_depth,"m with a total area of circa",int(total_pixels_basin*pixel_area),"km^2 consisting of",total_pixels_basin,"data points\n")
    for name,flux in zip(["raw_Osborn","rolling_mean_Osborn","raw_Shih","rolling_mean_Shih"],[raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih]):
        print(name)
        print("mean flux:", np.nanmean(flux),"mmol/m^2/d")
        print("max flux:",np.nanmin(flux),"mmol/m^2/d")
        flux_ratio(flux, distance, -1,pixels_edge, pixels_interior, diff_distance,pixel_area,weighted_mean_pixel_side_length)
        print("\n")

#plt.show()
