#TODO Python plot coordinates axis fraction/relation? non distorted images? (Mercator für kleine Region?)

import scipy.io as sio
import datetime as dt
import pathlib
import numpy as np
import netCDF4 as CDF
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import geopy.distance as geo #https://pypi.org/project/geopy/

#contains the MSS Data
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/emb217_mss_data"]#,"/home/ole/share-windows/emb177_mss_data/","/home/ole/share-windows/emb169_mss_data/MSS055/matlab/","/home/ole/share-windows/emb169_mss_data/MSS038/matlab/"]

#TODO
#define the figure
f1, axis = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = True)

def coord2float(degrees,minutes):
    return degrees.astype("float") + (minutes.astype("float"))/60

def coord2float2(degrees,minutes):
    return degrees + (minutes)/60
    
dset = CDF.Dataset("/home/ole/windows/all_data/bathymetry_data/iowtopo2_rev03.nc")
#print(dset)

z_topo = dset.variables["Z_WATER"]
XT_I = dset.variables["XT_I"]
YT_J = dset.variables["YT_J"]
#print(np.shape(z_topo))
#print(np.shape(XT_I))

read_positions = np.genfromtxt("/home/ole/Thesis/Preprocessing_TC_stations/station_positions_formatted", dtype = str)

#print(read_positions[:,1],read_positions[:,2])

longitude_list = coord2float(read_positions[:,1],read_positions[:,2])

latitude_list = coord2float(read_positions[:,3],read_positions[:,4])

emb169_lars = [coord2float2(57,07.984),coord2float2(20,31.971)]
print("emb169_lars",emb169_lars)

emb169_flach = [coord2float2(57,19.213),coord2float2(20,37.348)]
print("emb169_flach",emb169_flach)

emb169_tief = [coord2float2(57,19.182),coord2float2(20,35.934)]
print("emb169_tief",emb169_tief)

emb177_flach = [coord2float2(57,19.20),coord2float2(20,37.27)]
print("emb177_flach",emb177_flach)

emb177_tief = [coord2float2(57,19.232),coord2float2(20,36.011)]
print("emb177_tief",emb177_tief)

emb217_flach = [57.3200,20.6150]
print("emb217_flach",emb217_flach)

emb217_tief = [57.3200,20.600]
print("emb217_tief",emb217_tief)

#aspect_ratio = XT_I.size/YT_J.size

distance = geo.distance((emb169_tief[0],emb169_tief[1]),(emb169_flach[0],emb169_flach[1])).km
print("distance emb169 tief/flach in km = ",distance)

distance = geo.distance((emb177_tief[0],emb177_tief[1]),(emb177_flach[0],emb177_flach[1])).km
print("distance emb177 tief/flach in km = ",distance)

distance = geo.distance((emb169_tief[0],emb169_tief[1]),(emb177_tief[0],emb177_tief[1])).km
print("distance tief/tief in km = ",distance)

"""
#Plot with pcolormesh
#------------------------------------------------------
plt.pcolormesh(XT_I,YT_J,z_topo, vmin = -300)
plt.colorbar(label = "depth [m]")
plt.plot(emb169_lars[1],emb169_lars[0],"m.")
plt.plot(emb169_flach[1],emb169_flach[0],"r.")
plt.plot(emb169_tief[1],emb169_tief[0],"r.")
plt.plot(emb177_flach[1],emb177_flach[0],"k.")
plt.plot(emb177_tief[1],emb177_tief[0],"k.")

plt.ylim([56,58.2])
plt.xlim([17.8,21.8])

plt.xlabel("longitude [degrees]")
plt.ylabel("latitude [degrees]")
plt.title("adcp positions")
#plt.savefig("pictures/adcp_positions.png")
plt.show()
plt.close()
"""
#TODO unverzerrt plotten
#TODO limit to colorbar

#Plot with contourf
#---------------------------------------------------
long, lat = np.meshgrid(XT_I,YT_J)
bathymetrie = np.asarray(z_topo)

#replace the fill_values with NaN
bathymetrie[bathymetrie == -1.000e+34] = np.nan

#TODO replace this dirty trick
#bathymetrie[bathymetrie < -350] = -350

#itemindex = numpy.where(array==item)

axis.contourf(long, lat, bathymetrie, 15, extend = "min", vmin = -300, vmax = 0)

#axis.plot(emb169_lars[1],emb169_lars[0],"m.")
#axis.plot(emb169_flach[1],emb169_flach[0],"gD")
#axis.plot(emb169_tief[1],emb169_tief[0],"gD")
#axis.plot(emb177_flach[1],emb177_flach[0],"k.")
#axis.plot(emb177_tief[1],emb177_tief[0],"k.")
axis.plot(emb217_flach[1],emb217_flach[0],"k.")
axis.plot(emb217_tief[1],emb217_tief[0],"k.")

#print(longitude_list)
#print(latitude_list)


#axis.colorbar(label = "depth [m]")
#axis.ylim([56,58.2])
#axis.xlim([17.8,21.8])
#axis.cm.set_clim(-300,0)

axis.set_xlabel("longitude [degrees]")
axis.set_ylabel("latitude [degrees]")
axis.set_title("adcp positions")
#axis.savefig("adcp_positions.png")


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
    
        if DATAFILENAME == "TS11_TODL_merged.mat":
            continue
            
        print(cruisename,DATAFILENAME[:-4])
    
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
                               
        #print(sio.whosmat(FILENAME))
        data = sio.loadmat(datafile_path) 

        STA_substructure = data["STA"]
        CTD_substructure = data["CTD"]

        lat = STA_substructure["LAT"][0]
        lon = STA_substructure["LON"][0]

        pressure = CTD_substructure["P"][0]
        
        date = STA_substructure["date"][0]
 
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


        #print("lat",max(lat),min(lat))
        #print("lon",max(lon),min(lon))

        #remove data from a file, with overlapping positional points
        if (datafile_path == "/home/ole/share-windows/emb217_mss_data/TR1-8.mat"):
            lat = np.delete(lat,np.s_[33:47])
            lon = np.delete(lon,np.s_[33:47])
            distance = np.delete(distance,np.s_[33:47])
            
            pressure = np.delete(pressure,np.s_[33:47],axis=0)
            number_of_transects = np.shape(pressure)[-1]
            
        #remove excess data that already belongs to TS119       
        if cruisename == "emb169" and  DATAFILENAME[:-4] == "TS118":
            lat = lat[:21]
            lon = lon[:21]
            distance = distance[:21]
            pressure = pressure[:21]
            number_of_transects = np.shape(pressure)[-1]
        
        #removes the last data point, that dont seem to belong to the transect    
        if cruisename == "emb169" and  DATAFILENAME[:-4] == "TRR109":
            lat = lat[:-1]
            lon = lon[:-1]
            distance = distance[:-1]
            pressure = pressure[:-1]
            number_of_transects = np.shape(pressure)[-1]   

        #f1.set_size_inches(14,14)    
        #f1.set_size_inches(1.618*7.2,7.2)
            
        if cruisename == "emb169":
            axis.plot(lon,lat,"g")
            axis.plot(lon,lat,"gx")

        if cruisename == "emb177":
            axis.plot(lon,lat,"r")
            axis.plot(lon,lat,"rx")            
            
        if cruisename == "emb217":
            axis.plot(lon,lat,"b")
            axis.plot(lon,lat,"bx")   
                
        if DATAFILENAME[:-4] == "S106-1":
            axis.plot(lon,lat,"k")
            axis.plot(lon,lat,"kx")       

#axis.plot(longitude_list,latitude_list,"kD")           
plt.show()          
            
            
            
            
            
            
            
