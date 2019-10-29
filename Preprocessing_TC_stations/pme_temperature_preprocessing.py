import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import pylab as pl
import datetime as dt
import matplotlib.dates as mdates

LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/PME/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/PME/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/PME/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/PME/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/PME/data"]

#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/PME/data"]

#File with SensorID, Depth and measurement period information
sensor_positions_path = "/home/ole/thesis/Preprocessing_TC_stations/PME/pme_properties"
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


    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    
    #create ouput pictures, which will be filled later in the code
    
    #figure 1 for the measurements
    f1, axarr1 = plt.subplots(1)

    #figure 2 for the untrimmed data
    f2, axarr2 = plt.subplots(2)
    
    temperature_from_all_sensors = []
    label_list = []
    
    #loop over the files in a specified folder
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        sensor_ID = DATAFILENAME.split("_")[-1][3:-4]
        
        if DATAFILENAME == DATAFILENAMES[0]:
            print("####################################################\nfile:")
            print(datafile_path[25:])
        else:
            print("----------------------------------------------------\nfile:")
            print(datafile_path[25:])
        
        #exeption because of a different data structure
        if (datafile_path == "/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/PME/data/EMB169_TC2_PME041954.TXT"):
            temporary_data = np.genfromtxt(datafile_path,skip_header= 9, usecols = (0,5,6), delimiter = ",", encoding="iso8859_15")    
        
        elif (datafile_path == "/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/PME/data/EMB177_tc-tief_PME299629.txt"):
            print("EMB177_tc-tief_PME299629.txt is skipped, because its length doesn't fit to the other data")
            continue  
                
        #data from emb217 has a different coloumn structure
        elif FOLDERNAME in ["/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/PME/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/PME/data"]:
            temporary_data = np.genfromtxt(datafile_path,skip_header= 9, usecols = (0,4), delimiter = ",", encoding="iso8859_15")
            
        else:
            temporary_data = np.genfromtxt(datafile_path,skip_header= 9, usecols = (0,3), delimiter = ",", encoding="iso8859_15")
        

        
        fulltime = temporary_data[:,0]


        #search for measurement properties in the properties file
        for i in np.arange(np.shape(sensor_positions)[0]):

            #checks cruisename and sensor ID from the filename with the description in the properties file
            if (str(sensor_positions[i,0][0:6]) == str(cruisename) and str(sensor_positions[i,1]) == str(sensor_ID)):
                
                sensor_depth = sensor_positions[i,2]
                start_time = float(sensor_positions[i,3])
                stop_time = float(sensor_positions[i,4])
                #start_time = dt.datetime.strptime(sensor_positions[i,3], "%Y-%m-%d_%X")
                #stop_time = dt.datetime.strptime(sensor_positions[i,4], "%Y-%m-%d_%X")
                start_index = int(sensor_positions[i,5])
                stop_index = int(sensor_positions[i,6])
                
                print(str(cruisename),str(sensor_ID),sensor_depth,start_time,stop_time,start_index,stop_index)
                break
                 
        print("TEST:",np.arange(0,fulltime.size)[fulltime==start_time])
        print("TEST:",np.arange(0,fulltime.size)[fulltime==stop_time])
        
        trimmed_data =  temporary_data[start_index:stop_index,:]
        temperature = trimmed_data[:,1]
        
        #temperaure can't be reasonable negative
        temperature[temperature<0] = np.nan
        
        unix_time = fulltime[start_index:stop_index]
        #time axis of the measurement in coordinated universal time 
        utc = [dt.datetime.utcfromtimestamp(ts) for ts in unix_time]        
        utc = np.asarray(utc)    
        
            
        print("start",utc[0],start_time, unix_time[0] ==  start_time)
        print("stop",utc[-1],stop_time, unix_time[-1] ==  stop_time)
        
     
                 
        label = sensor_ID+" bottom + "+str(sensor_depth)+"m"

        temperature_from_all_sensors.append(temperature)
        label_list.append([sensor_ID,sensor_depth])

        #fill plots with data
        axarr1.plot(utc,temperature, label = label) #temperature plot
        
        axarr2[0].plot(fulltime,temporary_data[:,1], label = label)
        axarr2[0].axvline(unix_time[0], c = "k")
        axarr2[0].axvline(unix_time[-1], c= "k")
        
        axarr2[1].plot(temporary_data[:,1], label = label)
        axarr2[1].axvline(start_index, c = "k")
        axarr2[1].axvline(stop_index, c= "k")
        
    

    #save the trimmed data for easier analysis afterwards
    np.savez("/home/ole/thesis/Preprocessing_TC_stations/PME/data/PME_"+cruisename+"_"+flach_or_tief+"_temperature", temperature = np.asarray(temperature_from_all_sensors),  label_list = label_list, utc = utc)
    
    
    #set the xticks (dateformat and tick location)
    hfmt = mdates.DateFormatter('%d %b')
    axarr1.xaxis.set_major_locator(mdates.DayLocator())
    axarr1.xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
    axarr1.xaxis.set_major_formatter(hfmt)
    
    axarr1.set_ylabel("temperature [deg C]")
    axarr1.set_xlabel(utc[0].strftime("%Y"))
        
    title_fig1 = "PME "+cruisename+" "+flach_or_tief+" temperature"
    
    
    axarr1.set_title(title_fig1)
    axarr1.legend()
    
    f1.set_size_inches(12,7)

    plot1_name = "/home/ole/thesis/Preprocessing_TC_stations/PME/pictures/"+"PME_"+cruisename+"_"+flach_or_tief+" temperature"
    f1.savefig(plot1_name)
    
    #figure 2
    axarr2[0].set_ylabel("temperature [deg C]")
    axarr2[1].set_ylabel("temperature [deg C]")
    axarr2[1].set_xlabel(utc[0].strftime("%Y")) #label the axis with the corresponding year
      
    title_fig2 = "PME "+cruisename+" "+flach_or_tief+" untrimmed data comparison"
    
    axarr2[0].set_title(title_fig2)

    axarr2[0].legend()
    f2.set_size_inches(12,7)

    plot2_name = "/home/ole/thesis/Preprocessing_TC_stations/PME/pictures/"+"PME_"+cruisename+"_"+flach_or_tief+" untrimmed temperature comparison"
    f2.savefig(plot2_name)
      
#plt.show()    
