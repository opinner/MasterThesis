import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
#from methods_plots import *
import pylab as pl
import datetime as dt
import matplotlib.dates as mdates
from numpy.core.defchararray import add
import locale                                                           

print(locale.getlocale())

#ubuntu laptop
#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/RBR/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/RBR/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/RBR/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/RBR/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/RBR/data","/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Tief/RBR/data"]

#virtual machine
#LIST_OF_FOLDERS = ["/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_flach/RBR/data","/home/ole/windows/all_data/emb169/deployments/moorings/Peter_TC_tief/RBR/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-flach/RBR/data","/home/ole/windows/all_data/emb177/deployments/moorings/TC-tief/RBR/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/RBR/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/RBR/data"]

LIST_OF_FOLDERS = ["/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/RBR/data","/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/RBR/data"]

#depths of all the sensors, used for labeling in the plots
sensor_positions_path = "/home/ole/windows/Preprocessing_TC_stations/RBR/RBR_properties"
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
    f1, axarr1 = plt.subplots(1)

    #figure 2 for the check:
    f2, axarr2 = plt.subplots(2)
    
    #for the save file of the cleaned data
    data_from_all_sensors = []
    label_list = []
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        sensor_ID = DATAFILENAME.split("_")[-2]
        
        #exception for the files with depth data
        if DATAFILENAME == DATAFILENAMES[0]:

            print("########################################\nfile:")
            print(datafile_path[25:])
        
            skip_header= 25
            if FOLDERNAME == "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/RBR/data" or FOLDERNAME == "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/RBR/data":
                skip_header= 28
            

            temporary_data = np.genfromtxt(datafile_path,skip_header= skip_header, usecols = (0,1,2), encoding="iso8859_15", dtype= "str")

        else:
            print("--------------------------------------\nfile:")
            print(datafile_path[25:])
            skip_header = 22 #default value
            if DATAFILENAME == "EMB177_tc-tief_rbr_016172_eng.txt":
                skip_header = 25
                             
            if FOLDERNAME == "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Flach/RBR/data" or FOLDERNAME == "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/RBR/data":
                skip_header = 28
               
            temporary_data = np.genfromtxt(datafile_path,skip_header= skip_header, usecols = (0,1,2), encoding="iso8859_15", dtype= "str")
            
        locale.setlocale(locale.LC_TIME, "C")
        print(locale.getlocale())
        #convert the strings in the data to datetime objects
        days = temporary_data[:,0]
        hours = temporary_data[:,1]
        string_time = np.asarray(add(days,add("-",hours)))
        full_utc = np.asarray([dt.datetime.strptime(string_time[i], "%d-%b-%Y-%X.%f") for i in np.arange(string_time.size)])
        
        #Error in the time log of one sensor. Wrong about one hour
        #if datafile_path[25:] == "/emb217/deployments/moorings/TC_Flach/RBR/data/EMB217_TC-Chain-flach_016172_eng.txt":
        #    full_utc = full_utc - dt.timedelta(hours=1) 
       
        full_temperature = temporary_data[:,2].astype("float")
        #temperature can't be reasonable negative
        full_temperature[full_temperature<0] = np.nan

        
        
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
              
  
        print("TEST:",np.arange(0,full_utc.size)[full_utc==start_time])
        print("TEST:",np.arange(0,full_utc.size)[full_utc==stop_time])
        
        temperature = full_temperature[start_index:stop_index]
        utc = full_utc[start_index:stop_index]
        
        print("start",utc[0],start_time, utc[0] ==  start_time)
        print("stop",utc[-1],stop_time, utc[-1] ==  stop_time)
        print(utc.size)
        
        label = sensor_ID+" bottom + "+str(sensor_depth)+"m"
        label_list.append([sensor_ID,sensor_depth])
        
        data_from_all_sensors.append(temperature)
        print(np.shape(data_from_all_sensors))
        

        #data is too big to load at once, so I split it here and save the first part (the first 4 data files)
        if FOLDERNAME == "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/RBR/data" and np.shape(data_from_all_sensors)[0] >= 4:
            np.savez("/home/ole/windows/Preprocessing_TC_stations/RBR/data/RBR_"+cruisename+"_"+flach_or_tief+"_Part1_temperature", temperature = data_from_all_sensors, label_list = label_list, utc = utc)
            data_from_all_sensors = []

            
        #fill the plots with data
        axarr1.plot(utc,temperature, label = label)
    
        axarr2[0].plot(full_utc,full_temperature, label = label)
        axarr2[0].axvline(utc[0], c = "k")
        axarr2[0].axvline(utc[-1], c= "k")
        
        axarr2[1].plot(full_temperature)
        axarr2[1].axvline(start_index, c = "k")
        axarr2[1].axvline(stop_index, c= "k")
    
    
    #data_from_all_sensors = np.asarray(data_from_all_sensors)
    #save the trimmed data for easier analysis afterwards
    #TODO Better compressing or not?
    
    #np.savez("/home/ole/windows/Preprocessing_TC_stations/RBR/data/RBR_"+cruisename+"_"+flach_or_tief+"_temperature", temperature = data_from_all_sensors, label_list = label_list, utc = utc) 
 
    if FOLDERNAME == "/home/ole/windows/all_data/emb217/deployments/moorings/TC_Tief/RBR/data":
        np.savez("/home/ole/windows/Preprocessing_TC_stations/RBR/data/RBR_"+cruisename+"_"+flach_or_tief+"_Part2_temperature", temperature = data_from_all_sensors, label_list = label_list, utc = utc)
            
            
    else:
        np.savez("/home/ole/windows/Preprocessing_TC_stations/RBR/data/RBR_"+cruisename+"_"+flach_or_tief+"_temperature", temperature = data_from_all_sensors, label_list = label_list, utc = utc)       

     

        
    #axarr1[1].invert_yaxis()
    #axarr1[0].set_xlim(deploy_time,recovery_time)

    title_fig1 = "RBR "+cruisename+" "+flach_or_tief+" temperature"
    #title_fig1_2 = "RBR "+cruisename+" "+flach_or_tief+" depth"
    
    axarr1.set_ylabel("temperature [C]")
    axarr1.set_xlabel(utc[0].strftime("%Y"))
    axarr1.set_title(title_fig1)

    axarr1.legend()
    
    #format xticks
    hfmt = mdates.DateFormatter('%d %b')
    axarr1.xaxis.set_major_locator(mdates.DayLocator())
    axarr1.xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
    axarr1.xaxis.set_major_formatter(hfmt)
    
    f1.set_size_inches(12,7)

    plot1_save_windows_name = "/home/ole/windows/Preprocessing_TC_stations/RBR/pictures/"+"RBR_"+cruisename+"_"+flach_or_tief+"_temperature"
    plot1_save_ubuntu_name =  "/home/ole/Thesis/Preprocessing_TC_stations/RBR/pictures/"+"RBR_"+cruisename+"_"+flach_or_tief+"_temperature"
        
    f1.savefig(plot1_save_windows_name)
    f1.savefig(plot1_save_ubuntu_name)
    
    #figure 2
    axarr2[0].set_ylabel("temperature")
    axarr2[0].set_xlabel(utc[0].strftime("%Y")) #label the axis with the corresponding year
      
    title_fig2 = "RBR "+cruisename+" "+flach_or_tief+" untrimmed data comparison"
    
    axarr2[0].set_title(title_fig2)

    axarr2[0].legend()
    f2.set_size_inches(12,7)

    plot2_save_windows_name = "/home/ole/windows/Preprocessing_TC_stations/RBR/pictures/"+"RBR_"+cruisename+"_"+flach_or_tief+"_untrimmed_data_comparison"
    plot2_save_ubuntu_name =  "/home/ole/Thesis/Preprocessing_TC_stations/RBR/pictures/"+"RBR_"+cruisename+"_"+flach_or_tief+"_untrimmed_data_comparison"
        
    f2.savefig(plot2_save_windows_name)
    f2.savefig(plot2_save_ubuntu_name)
    
#plt.show()          
