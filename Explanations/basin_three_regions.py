import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import geopy.distance as geo
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

def z_from_p(pressure_values):
    import gsw.conversions 
    
    center_gotland_basin_lat = 57.0
    return gsw.conversions.z_from_p(pressure_values,center_gotland_basin_lat * np.ones(np.shape(pressure_values)))
    
    
transect_path = "/home/ole/windows/processed_mss/emb217/TR1-4.npz"
transect_data = np.load(transect_path)
transect_lat = transect_data["lat"] #Latitude of the profiles
transect_lon = transect_data["lon"] #Longitude of the profiles
        
dset = nc.Dataset("/home/ole/windows/all_data/bathymetry_data/iowtopo2_rev03.nc")
#print(dset)

z_topo = dset.variables["Z_WATER"]
XT_I = dset.variables["XT_I"]
YT_J = dset.variables["YT_J"]



#Plot with contourf
#---------------------------------------------------
lon, lat = np.meshgrid(XT_I,YT_J)
bathymetry = np.asarray(z_topo)

#replace the fill_values with NaN
bathymetry[bathymetry == -1.000e+34] = np.nan

#halocline_depths = z_from_p([73.23,58.16,72.15])
#mean_interval_edges = z_from_p([79.6,60.39,77.58])

halocline_depths = [-72.24,-57.52,-71.38]
mean_interval_edges = [-74.84,-57.85,-73.123]
regions_list = [0,10],[10,40],[40,999]]


#loop over the three cruises and their mean halocline depth
for cruise_index,cruise_name,halocline_depth,mean_interval_edge in zip([0,1,2],["emb169","emb177","emb217"],halocline_depths,mean_interval_edges):

    lon_0 = 19.8
    lat_0 = 57.1
    output_picture, axis = plt.subplots(1)
    map_ax = Basemap(ax = axis, width=230000,height=280000,resolution='f',projection='stere',lat_ts=50,lat_0=lat_0,lon_0=lon_0)

    #basin_shape = np.copy(bathymetry)
    basin_shape_list = bathymetry[bathymetry < halocline_depth]

    basin_lat_list = lat[bathymetry < halocline_depth]
    basin_lon_list = lon[bathymetry < halocline_depth]
    print(np.shape(basin_lat_list),np.shape(basin_lon_list),np.shape(basin_shape_list))

    #calculating the edge
    #--------------------------------------
    """
    left = np.argmin(np.abs(np.mean(lon,axis=0)-18.0))
    right = np.argmin(np.abs(np.mean(lon,axis=0)-21.0))
    bottom = np.argmin(np.abs(np.mean(lat,axis=1)-56.0))
    top = np.argmin(np.abs(np.mean(lat,axis=1)-59.0))
    """

    #left = np.argmin(np.abs(np.mean(lon,axis=0)-16.8))
    #right = np.argmin(np.abs(np.mean(lon,axis=0)-25))
    #bottom = np.argmin(np.abs(np.mean(lat,axis=1)-54.0))
    #top = np.argmin(np.abs(np.mean(lat,axis=1)-59.8))

    #only inside the rectangle defined by this values
    left = np.argmin(np.abs(np.mean(lon,axis=0)-18.5))
    right = np.argmin(np.abs(np.mean(lon,axis=0)-21))
    bottom = np.argmin(np.abs(np.mean(lat,axis=1)-56.5))
    top = np.argmin(np.abs(np.mean(lat,axis=1)-57.75))
        
    print(lon[0,left],lon[0,right],lat[bottom,0],lat[top,0])

    print("lat",np.min(np.mean(lat,axis=1)),np.max(np.mean(lat,axis=1)),np.shape(np.mean(lat,axis=1)))
    print("lon",np.min(np.mean(lon,axis=0)),np.max(np.mean(lon,axis=0)),np.shape(np.mean(lon,axis=0)))

    print(left,right,bottom,top)

    #create smaller matrix only with the area of interest
    basin_bath = bathymetry[bottom:top,left:right]
    basin_lat = lat[bottom:top,left:right]
    basin_lon = lon[bottom:top,left:right]
    basin_mask = (basin_bath < halocline_depth) #True for everything deeper that the halocline_depth


    #get the values below the halocline
    below_halocline_bathymetry = basin_bath[basin_mask]
    below_halocline_lat = basin_lat[basin_mask]
    below_halocline_lon = basin_lon[basin_mask]
    

    below_halocline_distance_to_ground = np.maximum(0,abs(below_halocline_bathymetry) - abs(mean_interval_edge))

    print("-"*50)
    print(np.shape(below_halocline_distance_to_ground),min(below_halocline_distance_to_ground),max(below_halocline_distance_to_ground))



    #plotting
    #-------------------------------
    levels = np.arange(-220,20,5)

    #cs2 = map_ax.plot(basin_lon, np.ma.masked_where(~basin_mask, basin_lat) , "r.", latlon = True)
    #cs2 = map_ax.plot(below_halocline_lon, below_halocline_lat , "r.", latlon = True)
    #cs3 = map_ax.plot(basin_lon, np.ma.masked_where(~edge_mask, basin_lat) , "kx", latlon = True)
    #print(np.shape(lon),np.shape(lat),np.shape(bathymetry))
    cs = map_ax.contourf(lon, lat, bathymetry,levels, extend = "min", cmap="Blues_r", latlon = True)           

    #"#ffc8c0","#f78a7e","#e54742"
      
    for region,color,label_name in zip([regions_list,['#fed98e','#fe9929','#cc4c02'],["edge region","intermediate region","interior region"]):
        region_mask = np.logical_and(below_halocline_distance_to_ground >= region[0],below_halocline_distance_to_ground < region[1])
        #print(np.shape(below_halocline_lat[region_mask]))
        #print(below_halocline_lat[region_mask])
        map_ax.plot(below_halocline_lon[region_mask], below_halocline_lat[region_mask], ".", c = color, marker = "s", markersize = 9, latlon = True, alpha = 0.7, label = label_name)
        
        """
        shape = (len(below_halocline_lon[region_mask]),len(below_halocline_lon[region_mask]))
        
        #mask = region_mask[region_mask.mask == False]
        x = below_halocline_lon #[below_halocline_lon.mask == False]
        y = below_halocline_lat #[below_halocline_lat.mask == False]
        

        #x = np.array(x.flatten().tolist())
        #y = np.array(y.flatten().tolist())
        #region_mask = np.array(region_mask.flatten().tolist())
        
        
        x =  x[region_mask]
        y =  y[region_mask]
        
        xx, yy = np.meshgrid(x, y)
        z = np.ones(np.shape(xx))
        #print(xx,yy,z)
        map_ax.contourf(yy, xx, z, [0,1,2], latlon = True) #, colors = color)
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

    #map_ax.drawmapscale(23.5,54.8,lon_0,lat_0,50,barstyle='fancy')

    # Add Colorbar
    cbar = map_ax.colorbar(cs, location='right', pad="2%")
    cbar.set_label("depth [m]",size = 14)

    # Add Title
    axis.set_title('Eastern Gotland Basin with hypoxic area below '+str(-halocline_depth)+" m", size = 16, weight = "bold")
    x,y = map_ax(18.45,57.4)
    axis.text(x,y,'Gotland',ha = "center", va = "center",rotation=40, size = 14 )
    x,y = map_ax(21.4,56.8)
    axis.text(x,y,'Latvia',ha = "center", va = "center",rotation=0, size = 14 )

    axis.set_title(cruise_name+": Halocline at "+str(abs(halocline_depth))+" m")
    axis.legend(loc = "lower left")

    #beamer figure sizes
    width = 6.2012 #f_iso.set_size_inches(6.2012,6.2012/1.618)
    height = 6.2012 #* 4/3 #1.618

    output_picture.set_size_inches(width,height)

    output_picture.tight_layout()
    output_picture.savefig("./"+cruise_name+"three_regions.pdf",dpi=600)
    output_picture.savefig("./"+cruise_name+"three_regions.png",dpi=600)




plt.show()
