import gsw.conversions #https://teos-10.github.io/GSW-Python/
import gsw 
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates

def SA_from_C(C, t, p, lon, lat):
    SP = gsw.conversions.SP_from_C(C, t, p)
    return gsw.SA_from_SP_Baltic(SP, lon, lat)
    
#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/Seabird/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/Seabird/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/Seabird/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/Seabird/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/Seabird/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/seabird/data"]

LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/Seabird/data"]


#depths of all the sensors, used for labeling in the plots
sensor_positions_path = "/home/ole/thesis/Preprocessing_TC_stations/Seabird/seabird_properties"
sensor_positions = np.genfromtxt(sensor_positions_path,skip_header= 1, usecols = (0,1,2,3,4,5,6), dtype = "str")


for FOLDERNAME in LIST_OF_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    cruisename = splitted_foldername[5]
    flach_or_tief = splitted_foldername[8]
    
    #go through all files of specified folder and select only files ending with .txt
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])
        if ((all_files_name[-4:] == ".txt") or (all_files_name[-4:] == ".TXT")):
            DATAFILENAMES.append(str(p.parts[-1]))


    #sort the files alphanumerically
    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    
    #create ouput pictures, which will be filled later in the code
    
    #figure 1 for the measurements
    f1, axarr1 = plt.subplots(3,sharex = True)

    #figure 2 for the check:
    f2, axarr2 = plt.subplots(2, sharex = True)
    
    #for the save file of the cleaned data
    data_from_all_sensors = []
    label_list = []
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        sensor_ID = DATAFILENAME.split("_")[-1][2:6]
        
        #exception for the files with depth data
        if DATAFILENAME == DATAFILENAMES[0]:
            print("########################################\nfile:")
        else:
            print("--------------------------------------\nfile:")
        print(datafile_path[25:])
        
        temporary_data = np.genfromtxt(datafile_path, comments = "%", usecols = (0,1,2,3,4,5,6,7,8,9), encoding="iso8859_15")
        
        
        #search for measurement properties in the file
        for i in np.arange(np.shape(sensor_positions)[0]):
            #print(str(sensor_positions[i,0][0:6]),str(cruisename),str(sensor_positions[i,1]),str(sensor_ID))
            
            #checks cruisename and sensor ID from the filename with the description in the properties file
            if (str(sensor_positions[i,0][0:6]) == str(cruisename) and str(sensor_positions[i,1]) == str(sensor_ID)):
                
                sensor_depth = sensor_positions[i,2]
                start_time = dt.datetime.strptime(sensor_positions[i,3], "%Y-%m-%d_%X")
                stop_time = dt.datetime.strptime(sensor_positions[i,4], "%Y-%m-%d_%X")
                start_index = int(sensor_positions[i,5])
                stop_index = int(sensor_positions[i,6])
                
                print(str(cruisename),str(sensor_ID),sensor_depth,start_time,stop_time,start_index,stop_index)
                break
                
            
        #assign the data to more readerfriendly arrays     
        full_temperature = temporary_data[:,0]
        full_conductivity = temporary_data[:,1]
        full_salinity = temporary_data[:,2]
        full_depth = temporary_data[:,3]
        
        #trim the data to the desired time span known from the properties file
        temperature = full_temperature[start_index:stop_index]
        conductivity = full_conductivity[start_index:stop_index]
        salinity = full_salinity[start_index:stop_index]
        depth = full_depth[start_index:stop_index]
        
        #temperature can't be reasonable negative
        temperature[temperature<0] = np.nan
        
        #convert the numbers for year, month, day, hour, minute, second to datetime objects
        full_utc = np.asarray([dt.datetime(int(temporary_data[i,4]),int(temporary_data[i,5]),int(temporary_data[i,6]),int(temporary_data[i,7]),int(temporary_data[i,8]),int(temporary_data[i,9])) for i in np.arange(np.shape(temporary_data)[0])])      

        #trim to the desired length
        utc = full_utc[start_index:stop_index]

        
        print(full_utc[0],full_utc[-1])
        print(start_time,stop_time)




        #unit conversion because of different standards for the conductivity
        if np.nanmean(conductivity) < 3:
            conductivity*= 10
        

        constant_pressure = gsw.p_from_z(-float(sensor_depth),54.32) #latitude of the transect
        #print("pressure",constant_pressure)
        computed_salinity = SA_from_C(conductivity,temperature,constant_pressure,20.6,54.32)
        #print("comparison:",np.nanmean(salinity),np.nanmean(computed_salinity))
        #print("check:",np.nanmean(conductivity))
        

              
        #TODO Start- und Stopzeit so waehlen, dass hier was rauskommt
        print("TEST:",np.arange(0,full_utc.size)[full_utc==start_time])
        print("TEST:",np.arange(0,full_utc.size)[full_utc==stop_time])
        
        temperature = full_temperature[start_index:stop_index]
        utc = full_utc[start_index:stop_index]
        
        print("start",utc[0],start_time, utc[0] ==  start_time)
        print("stop",utc[-1],stop_time, utc[-1] ==  stop_time)
        print(utc.size)
        
        label = sensor_ID+" bottom + "+str(sensor_depth)+"m"
        label_list.append([sensor_ID,sensor_depth])
        
         
            
        #fill the plots with data
        axarr1[0].plot(utc,temperature, label = label)
        axarr1[1].plot(utc,salinity, label = label)
        axarr1[1].plot(utc,computed_salinity, "--",label = "calc_"+label)
        axarr1[2].plot(utc,depth, label = label)
    
        axarr2[0].plot(full_utc,full_temperature, label = label)
        axarr2[0].axvline(utc[0], c = "k")
        axarr2[0].axvline(utc[-1], c= "k")
        
        axarr2[1].plot(full_temperature)
        axarr2[1].axvline(start_index, c = "k")
        axarr2[1].axvline(stop_index, c= "k")
    
    
    #data_from_all_sensors = np.asarray(data_from_all_sensors)
    
    
    #TODO save the data  
    #np.savez("/home/ole/thesis/Preprocessing_TC_stations/Seabird/data/seabird_"+cruisename+"_"+flach_or_tief+"_temperature", temperature = data_from_all_sensors, label_list = label_list, utc = utc)       
        
    #theoretically only the year from the last datafile, but all of them belong to the same measurements
    axarr1[2].set_xlabel(utc[0].strftime("%Y")) 
    
    #label    
    axarr1[0].set_ylabel("temperature [C]")
    axarr1[1].set_ylabel("conductivity [mS/cm]")
    axarr1[1].set_ylabel("salinity [SA]")
    axarr1[2].set_ylabel("depth [m]")

    title_fig1 = "Seabird "+cruisename+" "+flach_or_tief+" temperature"
    #title_fig1_2 = "Seabird "+cruisename+" "+flach_or_tief+" depth"
    
    axarr1[0].set_title(title_fig1)

    axarr1[0].legend()
    
    #format xticks
    hfmt = mdates.DateFormatter('%d %b')
    axarr1[0].xaxis.set_major_locator(mdates.DayLocator())
    axarr1[0].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
    axarr1[0].xaxis.set_major_formatter(hfmt)
    
    f1.set_size_inches(12,7)

    plot1_name = "/home/ole/thesis/Preprocessing_TC_stations/Seabird/pictures/"+"seabird_"+cruisename+"_"+flach_or_tief+" temperature"
    f1.savefig(plot1_name)
    
    #figure 2
    axarr2[0].set_ylabel("temperature")
    axarr2[0].set_xlabel(utc[0].strftime("%Y")) #label the axis with the corresponding year
      
    title_fig2 = "PME "+cruisename+" "+flach_or_tief+" untrimmed data comparison"
    
    axarr2[0].set_title(title_fig2)

    axarr2[0].legend()
    f2.set_size_inches(12,7)

    plot2_name = "/home/ole/thesis/Preprocessing_TC_stations/Seabird/pictures/"+"Seabird_"+cruisename+"_"+flach_or_tief+" untrimmed data comparison"
    f2.savefig(plot2_name)

plt.show()          
