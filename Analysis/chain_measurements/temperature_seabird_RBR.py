import pathlib
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates
from operator import itemgetter

FOLDERNAME = "/home/ole/thesis/Preprocessing_TC_stations/Seabird/data"

#File with SensorID, Depth and measurement period information
sensor_positions_path = "/home/ole/thesis/Preprocessing_TC_stations/Seabird/seabird_properties"
sensor_positions = np.genfromtxt(sensor_positions_path,skip_header= 1, usecols = (0,1,2,3,4,5,6), dtype = "str")

path = pathlib.Path(FOLDERNAME)


#splitted_foldername = FOLDERNAME.split("/")

#cruisename = splitted_foldername[5]
#flach_or_tief = splitted_foldername[8]

DATAFILENAMES = []
#go through all files of specified folder and select only files ending with "temperature.npz"
for p in path.iterdir():
    FILENAME = str(p.parts[-1])
    splitted_filename = FILENAME.split("_")
    ending = splitted_filename[-1]
    
    if (ending == "temperature.npz"): 
        DATAFILENAMES.append(FILENAME)


#DATAFILENAMES = sorted(DATAFILENAMES) 
print(DATAFILENAMES)

#create ouput pictures, which will be filled later in the code

#figure 1 for TC_Flach
f1, axarr1 = plt.subplots(nrows = 3, ncols = 2)#, sharey= "col")

#figure 2 for TC_Tief
f2, axarr2 = plt.subplots(3)

max_time_delta = dt.timedelta(minutes=1)

#loop over the the seabird files in a specified folder
for DATAFILENAME in DATAFILENAMES:
    datafile_path = FOLDERNAME+"/"+DATAFILENAME
    
    print(datafile_path)
    
    splitted_filename = DATAFILENAME.split("_")

    device = splitted_filename[0]
    cruisename = splitted_filename[1]
    flach_or_tief = splitted_filename[-2]
    ending = splitted_filename[-1]
    print(device,cruisename,flach_or_tief,ending)

    seabird_data = np.load(datafile_path)
    print(seabird_data.files)
    
    seabird_temperature = seabird_data["temperature"]
    seabird_label_list = seabird_data["label_list"]
    seabird_depth = seabird_label_list[:,1].astype("float")
    seabird_utc = seabird_data["utc"]
    
    print(np.shape(seabird_temperature))
    print(seabird_utc[0],seabird_utc[-1])
    
    print("seabird_label_list")
    print(seabird_label_list)
    
    #sort the depth array in the case the depth are not increasing
    sorted_depth = np.asarray(sorted(seabird_depth))
    number_of_depth_bins = seabird_depth.size
    
    #Prepend a 0 to the beginning (index 0)
    #This is needed to plot all the depths later in the pcolormesh plot
    plot_edges_depth = np.insert(sorted_depth,0,0) 
    
    #print(seabird_label_list)
    #print(np.reshape(seabird_label_list[:,1],(number_of_depth_bins,1)))
    
    #stack the index to the label_list
    sorting_array = np.hstack((np.reshape(seabird_label_list[:,1],(number_of_depth_bins,1)),np.reshape(np.arange(number_of_depth_bins),(number_of_depth_bins,1)))).astype("float")
    
    print("shape seabird")
    print(np.shape(seabird_temperature))#,np.shape(seabird_label_list),np.shape(seabird_utc))
    
    """
    print("original:")
    for depth_index in range(np.shape(seabird_temperature)[0]):
        print(depth_index,seabird_depth[depth_index],np.nanmean(seabird_temperature[depth_index,:]))
    """
    

    
    #print(sorting_array)
    #print(sorting_key)
    
    
    #sort after the first column, eg the depths (sorts automatically the indices)
    #from https://stackoverflow.com/questions/2173797/how-to-sort-2d-array-by-row-in-python
    sorting_key = np.asarray(sorted(sorting_array, key= itemgetter(0)))
    
    #sort the data with the help of the sorting key
    sorted_temperature = []
    for index in sorting_key[:,1]:
        sorted_temperature.append(seabird_temperature[int(index),:])
    sorted_temperature = np.asarray(sorted_temperature)
    
    """
    print("sorted")
    for depth_index in range(np.shape(sorted_temperature)[0]):
        print(depth_index,sorted_depth[depth_index],np.nanmean(sorted_temperature[depth_index,:]))
 
        
    print("------------------------------------------------------------")    
    """    
    #---------------------------------------------------------------------------------------------------------------------------------
    #load the corresponding RBR data
    RBR_FOLDERNAME = "/home/ole/thesis/Preprocessing_TC_stations/RBR/data"

    #File with SensorID, Depth and measurement period information
    sensor_positions_path = "/home/ole/thesis/Preprocessing_TC_stations/RBR/RBR_properties"
    sensor_positions = np.genfromtxt(sensor_positions_path,skip_header= 1, usecols = (0,1,2,3,4,5,6), dtype = "str")

    path = pathlib.Path(FOLDERNAME)


    if cruisename == "emb169":
        RBR_FILENAME = "/home/ole/thesis/Preprocessing_TC_stations/RBR/data/RBR_"+cruisename+"_Peter_TC_"+flach_or_tief+"_temperature.npz"
    
    elif cruisename == "emb177":
        RBR_FILENAME = "/home/ole/thesis/Preprocessing_TC_stations/RBR/data/RBR_"+cruisename+"_"+flach_or_tief+"_temperature.npz"
    
    elif cruisename == "emb217":
        RBR_FILENAME = "/home/ole/thesis/Preprocessing_TC_stations/RBR/data/RBR_"+cruisename+"_TC_"+flach_or_tief+"_temperature.npz"
            
    else:    
        print("cruisename not found")
        
    print(RBR_FILENAME)
    RBR_data = np.load(RBR_FILENAME)
    print(RBR_data.files)
    
    try:
        RBR_temperature = RBR_data["temperature"]
    except KeyError:
        RBR_temperature = RBR_data["data"]
    RBR_label_list = RBR_data["label_list"]
    RBR_depth = RBR_label_list[:,1].astype("float")
    RBR_utc = RBR_data["utc"]
    
    print("RBR_label_list")
    print(RBR_label_list)
    
    # RBR should be the denser data, therefore RBR gets interpolated
    assert(seabird_utc.size < RBR_utc.size)
    
    print("RBR shape")
    print(np.shape(RBR_temperature))      
    
    seabird_interpolation_time = np.asarray([dt.datetime.timestamp(x) for x in seabird_utc])
    RBR_interpolation_time = np.asarray([dt.datetime.timestamp(x) for x in RBR_utc])
    print("Interpolation from ",np.shape(RBR_interpolation_time)," points to ",np.shape(seabird_interpolation_time)," points")
    
    #interpolate the RBR Data to the seabird time axis 
    # loop through the depths
    for i in range(RBR_depth.size):
        current_depth = RBR_depth[i]
        interpolated_RBR_temperature = np.interp(seabird_interpolation_time,RBR_interpolation_time,RBR_temperature[i,:])
        #print(np.shape(RBR_interpolation_time),np.shape(RBR_temperature[i,:]),,np.shape(interpolated_RBR_temperature))
    
    
        #Plots to compare the original RBR data to the interpolated data
        if i == 0:
            if cruisename == "emb169":
                plot_number = 0
            elif cruisename == "emb177":
                plot_number = 1
            elif cruisename == "emb217":
                plot_number = 2
            else:
                print("ERROR: ID der Ausfahrt nicht bestimmbar")
        
            #plot the orignal data with the points in black, line in blue
            axarr2[plot_number].plot(RBR_interpolation_time,RBR_temperature[i,:],"k.")
            axarr2[plot_number].plot(RBR_interpolation_time,RBR_temperature[i,:],"b")
            
            #plot the interpolated data with points in red, line in yellow
            axarr2[plot_number].plot(seabird_interpolation_time,interpolated_RBR_temperature,"r.")
            axarr2[plot_number].plot(seabird_interpolation_time,interpolated_RBR_temperature,"y")

      
        #print(current_depth)
        #print("####")
        #print(sorted_depth)
        
        
        #insert the interpolated data at the right index   
        for j in range(sorted_depth.size):
            
            if current_depth == sorted_depth[j]:
                print("###################################")
                print("two datasets for ",current_depth," m")
                print("###################################")
                break #basically skips the insertion
            
        
        
            if current_depth < sorted_depth[j]:
            
                sorted_depth = np.insert(sorted_depth,j,current_depth)
                
                print(np.shape(sorted_temperature),np.shape(interpolated_RBR_temperature))
                
                sorted_temperature = np.insert(sorted_temperature,j,interpolated_RBR_temperature,axis = 0)
                break #point of insertion found, loop can be broken
                
            if  sorted_depth[j] ==  sorted_depth[-1]:
                sorted_depth = np.append(sorted_depth,current_depth)
                #print(np.shape(sorted_temperature),np.shape(interpolated_RBR_temperature))
                sorted_temperature = np.append(sorted_temperature,np.asarray([interpolated_RBR_temperature]),axis = 0)

        
        #print("####")
        #print(sorted_depth)    
    
    print("RBR----------------CHECK-----------------------")
    for depth_index in range(np.shape(sorted_temperature)[0]):
        
        print(depth_index,sorted_depth[depth_index],np.mean(sorted_temperature[depth_index,:]))
            
        
    #check if sorted depth is still sorted, eg all insertions were correct    
    assert(np.all(np.diff(sorted_depth) > 0))

    
    #-----------------------------------------------------------------------------------
    #Code for plotting
    
    
    #Prepend a 0 to the beginning (index 0)
    plot_edges_depth = np.insert(sorted_depth,0,0) 
    utc = np.copy(seabird_utc)
    
    #sort the data to the right plot position
    if flach_or_tief.lower() == "flach" or flach_or_tief.lower() == "tc-flach":
        x_position = 0
        if cruisename == "emb169":
            y_position = 0
            start_flach_2017 = utc[0]
        elif cruisename == "emb177":
            y_position = 1
            start_flach_2018 = utc[0]
        elif cruisename == "emb217":
            y_position = 2
            start_flach_2019 = utc[0]
            
            #because emb217 had the longest measurement period
            max_time_delta = max(max_time_delta,utc[-1]-utc[-0])
        else:
            print("ERROR: ID der Ausfahrt nicht bestimmbar")
        
    elif flach_or_tief.lower() == "tief" or flach_or_tief.lower() == "tc-tief":
        x_position = 1
        if cruisename == "emb169":
            y_position = 0
            start_tief_2017 = utc[0]
        elif cruisename == "emb177":
            y_position = 1
            start_tief_2018 = utc[0]
        elif cruisename == "emb217":
            y_position = 2
            start_tief_2019 = utc[0]
            
            #because emb217 had the longest measurement period
            max_time_delta = max(max_time_delta,utc[-1]-utc[-0])
        else:
            print("ERROR: ID der Ausfahrt nicht bestimmbar")
        
    else:
        print(flach_or_tief.lower())
        print("ERROR: Ort (flach_or_tief) nicht feststellbar")    
            

            
            
    #fill figure 1 with data
    img1_1 = axarr1[y_position,x_position].pcolormesh(utc,plot_edges_depth,sorted_temperature, cmap = plt.cm.RdYlBu_r)

    #add a colorbar to every plot
    f1.colorbar(img1_1, ax=axarr1[y_position,x_position])



#Preferences of the plots

#figure 1
#set the xticks (dateformat and tick location)
hfmt = mdates.DateFormatter('%d %b')
for row in range(3):
    for column in range(2):
        axarr1[row,column].xaxis.set_major_locator(mdates.DayLocator())
        axarr1[row,column].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
        axarr1[row,column].xaxis.set_major_formatter(hfmt)

        #axarr1[row,column].set_ylim(bottom = 0, top = 80)
        #print("xlim:",row,column,axarr1[row,column].get_xlim())





"""
for row in range(3):
    for column in range(2):
        print("xlim:",row,column,axarr1[row,column].get_xlim())
"""

axarr1[0,0].set_xlabel("2017")
axarr1[0,1].set_xlabel("2017")
axarr1[1,0].set_xlabel("2018")
axarr1[1,1].set_xlabel("2018")
axarr1[2,0].set_xlabel("2019")        
axarr1[2,1].set_xlabel("2019")        


title_fig1 = "flach"
title_fig2 = "tief"


#f1.subplots_adjust(right=0.9)
#cbar_ax = f1.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
#f1.colorbar(img1_2, cax=cbar_ax).set_label("velocity [m/s]") 

#f2.subplots_adjust(right=0.9)
#cbar_ax = f2.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
#f2.colorbar(img2_2, cax=cbar_ax).set_label("velocity [m/s]") 

        
axarr1[0,0].set_title(title_fig1)
axarr1[0,1].set_title(title_fig2)    
     
f1.suptitle("temperature from RBR/seabird comparison")
f2.suptitle("comparison original data to interpolated data")

f1.set_size_inches(18,10.5)
f1.tight_layout(rect=[0, 0.03, 1, 0.95])

f2.set_size_inches(18,10.5)
f2.tight_layout(rect=[0, 0.03, 1, 0.99])
    
#Save the plot as png
plot1_name = "./pictures/"+"temperature_RBR_seabird_comparison2" 
f1.savefig(plot1_name)

 
#preferences of the x-limit
axarr1[0,0].set_xlim(left = mdates.date2num(start_flach_2017), right = mdates.date2num(start_flach_2017 + max_time_delta))
axarr1[0,1].set_xlim(left = mdates.date2num(start_tief_2017), right = mdates.date2num(start_tief_2017 + max_time_delta))
axarr1[1,0].set_xlim(left = mdates.date2num(start_flach_2018), right = mdates.date2num(start_flach_2018 + max_time_delta))
axarr1[1,1].set_xlim(left = mdates.date2num(start_tief_2018), right = mdates.date2num(start_tief_2018 + max_time_delta))
axarr1[2,0].set_xlim(left = mdates.date2num(start_flach_2019), right = mdates.date2num(start_flach_2019 + max_time_delta))
axarr1[2,1].set_xlim(left = mdates.date2num(start_tief_2019), right = mdates.date2num(start_tief_2019 + max_time_delta))


#Save the changend plot again as png
plot1_name = "./pictures/"+"temperature_RBR_seabird_comparison_same_scale2" 
f1.savefig(plot1_name)
  
plt.show()
  
       
    
    
    
    
