#################################
#TODO Calibration               #
#TODO Temperature measurement   #
#TODO save data                 #
#################################
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import pylab as pl
import datetime as dt
import matplotlib.dates as mdates


#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/PME/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/PME/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/PME/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/PME/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/PME/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/PME/data"]
#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/PME/data"]

#virtual machine
LIST_OF_FOLDERS = ["/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_flach/PME/data","/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_tief/PME/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-flach/PME/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-tief/PME/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/PME/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/PME/data"]
#LIST_OF_FOLDERS = ["/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/PME/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/PME/data"]

#File with SensorID, Depth and measurement period information
sensor_positions_path = "/home/ole/windows/Preprocessing_TC_stations/PME/pme_properties"
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
    f1, axarr1 = plt.subplots(2, sharex=True)

    #figure 2 for the untrimmed data
    f2, axarr2 = plt.subplots(2)
    
    concentration_from_all_sensors = []
    saturation_from_all_sensors = []
    label_list = []
    
    #loop over the files in a specified folder
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        sensor_ID = DATAFILENAME.split("_")[-1][3:-4]
        
        if DATAFILENAME == DATAFILENAMES[0]:
            print("####################################################\nfile:")
            print(datafile_path[27:])
        else:
            print("----------------------------------------------------\nfile:")
            print(datafile_path[27:])
        
        #exeption because of a different data structure
        if (datafile_path == "/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_flach/PME/data/EMB169_TC2_PME041954.TXT"):
            temporary_data = np.genfromtxt(datafile_path,skip_header= 9, usecols = (0,5,6), delimiter = ",", encoding="iso8859_15")    
        
        elif (datafile_path == "/home/ole/windows/all_data/emb177/deployments/moorings/TC-tief/PME/data/EMB177_tc-tief_PME299629.txt"):
            print("EMB177_tc-tief_PME299629.txt is skipped, because its length doesn't fit to the other data")
            continue  
                
        #data from emb217 has a different coloumn structure
        elif FOLDERNAME in ["/home/ole/windows/all_data/emb177/deployments/moorings/TC-tief/PME/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/PME/data"]:
            temporary_data = np.genfromtxt(datafile_path,skip_header= 9, usecols = (0,5,6), delimiter = ",", encoding="iso8859_15")
            
        else:
            temporary_data = np.genfromtxt(datafile_path,skip_header= 9, usecols = (0,4,5), delimiter = ",", encoding="iso8859_15")
        
        #oxygen levels and saturation can't be reasonable negative
        temporary_data[temporary_data<0] = np.nan
        
        fulltime = temporary_data[:,0]


        #search for measurement properties in the file
        for i in np.arange(np.shape(sensor_positions)[0]):

            print(str(sensor_positions[i,0][0:6]),str(cruisename),str(sensor_positions[i,1]),str(sensor_ID))

            #checks cruisename and sensor ID from the filename with the description in the properties file
            if (str(sensor_positions[i,0][0:6]) == str(cruisename) and str(sensor_positions[i,1]) == str(sensor_ID)):
                
                sensor_depth = sensor_positions[i,2]

                try:
                    start_time = dt.datetime.strptime(sensor_positions[i,3], "%Y-%m-%d_%X")
                    stop_time = dt.datetime.strptime(sensor_positions[i,4], "%Y-%m-%d_%X")
                except ValueError:
                    start_time = float(sensor_positions[i,3])
                    stop_time = float(sensor_positions[i,4])
                    
                start_index = int(sensor_positions[i,5])
                stop_index = int(sensor_positions[i,6])
                
                print(str(cruisename),str(sensor_ID),sensor_depth,start_time,stop_time,start_index,stop_index)
                break
                 
        print("TEST:",np.arange(0,fulltime.size)[fulltime==start_time])
        print("TEST:",np.arange(0,fulltime.size)[fulltime==stop_time])
        
        trimmed_data =  temporary_data[start_index:stop_index,:]
        O_concentration = trimmed_data[:,1]
        O_saturation = trimmed_data[:,2]
        
        unix_time = fulltime[start_index:stop_index]
        #time axis of the measurement in coordinated universal time 
        utc = [dt.datetime.utcfromtimestamp(ts) for ts in unix_time]        
        utc = np.asarray(utc)    
        
        
        """
        #calibration through linear stretching of the data
        if FOLDERNAME in ["/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/PME/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/PME/data"]:
            #Assumption that during the calibration measurement the lowest 
            #measurement corresponds to 0 saturationand the highest to 100  
            max_saturation = np.max(temporary_data[:start_index,2])
            min_saturation = np.min(temporary_data[:start_index,2])
            
            m = 100/(max_saturation-min_saturation) 
            b = -100*min_saturation/(max_saturation-min_saturation)
            
            print("Max/Min",np.max(temporary_data[:start_index,2]),np.min(temporary_data[:start_index,2]))
            trimmed_data[:,1:2] = trimmed_data[:,1:2] * m + b
            
            print("Max/Min",np.max(temporary_data[:start_index,2]),np.min(temporary_data[:start_index,2]))
        """
            
        print("start",utc[0],start_time, utc[0] ==  start_time)
        print("stop",utc[-1],stop_time, utc[-1] ==  stop_time)
        
     
                 
        label = sensor_ID+" bottom + "+str(sensor_depth)+"m"

        concentration_from_all_sensors.append(O_concentration)
        saturation_from_all_sensors.append(O_saturation)
        label_list.append([sensor_ID,sensor_depth])

        #fill plots with data
        axarr1[0].plot(utc,O_concentration, label = label) #oxygen concentration
        axarr1[1].plot(utc,O_saturation, label = label) #oxygen saturation
        
        
        axarr2[0].plot(fulltime,temporary_data[:,2], label = label)
        axarr2[0].axvline(unix_time[0], c = "k")
        axarr2[0].axvline(unix_time[-1], c= "k")
        
        axarr2[1].plot(temporary_data[:,2], label = label)
        axarr2[1].axvline(start_index, c = "k")
        axarr2[1].axvline(stop_index, c= "k")
        
    

    #save the trimmed data for easier analysis afterwards
    np.savez("/home/ole/windows/Preprocessing_TC_stations/PME/data/PME_"+cruisename+"_"+flach_or_tief+"_oxygen", concentration = np.asarray(concentration_from_all_sensors), saturation = np.asarray(saturation_from_all_sensors),  label_list = label_list, utc = utc)
    
    
    #set the xticks (dateformat and tick location)
    hfmt = mdates.DateFormatter('%d %b')
    axarr1[1].xaxis.set_major_locator(mdates.DayLocator())
    axarr1[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
    axarr1[1].xaxis.set_major_formatter(hfmt)
    
    axarr1[0].set_ylabel("dissolved oxygen [mg/l]")
    axarr1[1].set_ylabel("dissolved oxygen saturation [%]")
    axarr1[1].set_xlabel(utc[0].strftime("%Y"))
        
    title_fig1_1 = "PME "+cruisename+" "+flach_or_tief+" dissolved oxygen"
    title_fig1_2 = "PME "+cruisename+" "+flach_or_tief+" dissolved oxygen saturation"
    
    axarr1[0].set_title(title_fig1_1)
    axarr1[1].set_title(title_fig1_2)

    axarr1[0].legend()
    axarr1[1].legend()
    f1.set_size_inches(12,7)

    plot1_save_windows_name = "/home/ole/windows/Preprocessing_TC_stations/PME/pictures/"+"PME_"+cruisename+"_"+flach_or_tief+" oxygen"
    plot1_save_ubuntu_name =  "/home/ole/Thesis/Preprocessing_TC_stations/PME/pictures/"+"PME_"+cruisename+"_"+flach_or_tief+" oxygen"
        
    f1.savefig(plot1_save_windows_name)
    f1.savefig(plot1_save_ubuntu_name)
        
    #figure 2
    axarr2[0].set_ylabel("dissolved oxygen [mg/l]")
    axarr2[0].set_xlabel(utc[0].strftime("%Y")) #label the axis with the corresponding year
      
    title_fig2 = "PME "+cruisename+" "+flach_or_tief+" untrimmed data comparison"
    
    axarr2[0].set_title(title_fig2)

    axarr2[0].legend()
    f2.set_size_inches(12,7)

    plot2_save_windows_name = "/home/ole/windows/Preprocessing_TC_stations/PME/pictures/"+"PME_"+cruisename+"_"+flach_or_tief+" untrimmed data comparison"
    plot2_save_ubuntu_name =  "/home/ole/Thesis/Preprocessing_TC_stations/PME/pictures/"+"PME_"+cruisename+"_"+flach_or_tief+" untrimmed data comparison"
    
    f2.savefig(plot2_save_windows_name)
    f2.savefig(plot2_save_ubuntu_name)
      
#plt.show()    
