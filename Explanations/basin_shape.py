import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import geopy.distance as geo

transect_path = "/home/ole/windows/processed_mss/emb217/TR1-4.npz"
transect_data = np.load(transect_path)
transect_lat = transect_data["lat"] #Latitude of the profiles
transect_lon = transect_data["lon"] #Longitude of the profiles
        
dset = nc.Dataset("/home/ole/windows/all_data/bathymetry_data/iowtopo2_rev03.nc")
#print(dset)

z_topo = dset.variables["Z_WATER"]
XT_I = dset.variables["XT_I"]
YT_J = dset.variables["YT_J"]

lon_0 = 20.1
lat_0 = 57.1
output_picture, axis = plt.subplots(1)
map_ax = Basemap(ax = axis, width=550000,height=620000,resolution='c',projection='stere',lat_ts=50,lat_0=lat_0,lon_0=lon_0)

#Plot with contourf
#---------------------------------------------------
lon, lat = np.meshgrid(XT_I,YT_J)
bathymetry = np.asarray(z_topo)

#replace the fill_values with NaN
bathymetry[bathymetry == -1.000e+34] = np.nan

for index,set_depth in enumerate([-71.07,-58.52,-70.35]):

    #set the desired depth of the <<halocline>> 
    #set_depth = -70

    #basin_shape = np.copy(bathymetry)
    basin_shape_list = bathymetry[bathymetry < set_depth]

    basin_lat_list = lat[bathymetry < set_depth]
    basin_lon_list = lon[bathymetry < set_depth]
    print(np.shape(basin_lat_list),np.shape(basin_lon_list),np.shape(basin_shape_list))

    #calculating the edge
    #--------------------------------------
    """
    left = np.argmin(np.abs(np.mean(lon,axis=0)-18.0))
    right = np.argmin(np.abs(np.mean(lon,axis=0)-21.0))
    bottom = np.argmin(np.abs(np.mean(lat,axis=1)-56.0))
    top = np.argmin(np.abs(np.mean(lat,axis=1)-59.0))
    """

    left = np.argmin(np.abs(np.mean(lon,axis=0)-16.8))
    right = np.argmin(np.abs(np.mean(lon,axis=0)-25))
    bottom = np.argmin(np.abs(np.mean(lat,axis=1)-54.0))
    top = np.argmin(np.abs(np.mean(lat,axis=1)-59.8))

    #only inside the rectangle defined by this values
    #left = np.argmin(np.abs(np.mean(lon,axis=0)-18.5))
    #right = np.argmin(np.abs(np.mean(lon,axis=0)-21))
    #bottom = np.argmin(np.abs(np.mean(lat,axis=1)-56.5))
    #top = np.argmin(np.abs(np.mean(lat,axis=1)-57.75))
        
    print(lon[0,left],lon[0,right],lat[bottom,0],lat[top,0])

    print("lat",np.min(np.mean(lat,axis=1)),np.max(np.mean(lat,axis=1)),np.shape(np.mean(lat,axis=1)))
    print("lon",np.min(np.mean(lon,axis=0)),np.max(np.mean(lon,axis=0)),np.shape(np.mean(lon,axis=0)))

    print(left,right,bottom,top)

    #create smaller matrix only with the area of interest
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

    print("lat",mean_lat,mean_lat_diff)
    print("lon",mean_lon,mean_lon_diff)
    print("result",x_distance,y_distance,pixel_area)


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

    print(np.shape(lat),np.shape(lon),np.shape(bathymetry))
    print(np.shape(basin_lat),np.shape(basin_lon),np.shape(basin_bath))
    print("basin lat",np.min(np.mean(basin_lat,axis=1)),np.max(np.mean(basin_lat,axis=1)),np.shape(np.mean(basin_lat,axis=1)))
    print("basin lon",np.min(np.mean(basin_lon,axis=0)),np.max(np.mean(basin_lon,axis=0)),np.shape(np.mean(basin_lon,axis=0)))

    assert(np.all(np.shape(basin_mask) == np.shape(basin_lon)))

    pixels_edge =  len(basin_lat[edge_mask])
    pixels_basin = len(basin_lat[basin_mask])


    weighted_mean_pixel_side_length = (horizontal_edges * x_distance + vertical_edges * y_distance)/(horizontal_edges + vertical_edges)
            



    ratio = pixels_edge/(pixels_basin-pixels_edge)
    print("\nbasin:",pixels_basin,"\tedge",pixels_edge,"\tratio",ratio,"\tarea",pixels_basin*pixel_area,"\n")








    #plotting
    #-------------------------------
    levels = np.arange(-220,20,5)

    cs2 = map_ax.plot(basin_lon, np.ma.masked_where(~basin_mask, basin_lat) , "r.", latlon = True)
    #cs3 = map_ax.plot(basin_lon, np.ma.masked_where(~edge_mask, basin_lat) , "kx", latlon = True)
    cs = map_ax.contourf(lon, lat, bathymetry,levels, extend = "min", cmap="Blues_r", latlon = True)           


    map_ax.plot(basin_lon[edge_mask], basin_lat[edge_mask], "k.", latlon = True)
    """
    for pointx,pointy,neighbours in zip(basin_lon[edge_mask], basin_lat[edge_mask],edge_array[edge_mask]):
        #map_ax.plot(pointx, pointy,"kx", latlon = True)
        map_ax.plot(pointx, pointy , marker = '$'+str(int(neighbours))+'$', color = "k", latlon = True)
    """    


    map_ax.plot(transect_lon,transect_lat,c = "tab:green", latlon = True,  linewidth = 3, label = "Transect")

    # Add Grid Lines
    map_ax.drawparallels(np.arange(0, 70, 1.), labels=[1,0,0,0], fontsize=10)
    map_ax.drawmeridians(np.arange(0, 100, 2), labels=[0,0,0,1], fontsize=10)

    # Add Coastlines, States, and Country Boundaries
    map_ax.drawcoastlines()
    #m.drawstates()
    map_ax.drawcountries()
    map_ax.fillcontinents(color="lightgrey")

    map_ax.drawmapscale(23.5,54.8,lon_0,lat_0,50,barstyle='fancy')

    # Add Colorbar
    cbar = map_ax.colorbar(cs, location='right', pad="2%")
    cbar.set_label("depth [m]",size = 14)

    # Add Title
    axis.set_title('Eastern Gotland Basin with hypoxic area below '+str(set_depth)+" m", size = 16, weight = "bold")
    x,y = map_ax(18.45,57.4)
    axis.text(x,y,'Gotland',ha = "center", va = "center",rotation=40, size = 14 )
    x,y = map_ax(21.7,57.2)
    axis.text(x,y,'Latvia',ha = "center", va = "center",rotation=0, size = 14 )

    axis.legend()

    output_picture.set_size_inches(9,10.5)

    output_picture.tight_layout()
    output_picture.savefig("./basin_shape"+str(index),dpi=300)





    #plt.show()
