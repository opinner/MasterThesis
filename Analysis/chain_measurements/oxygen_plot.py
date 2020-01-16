######################################################
#TODO 
#Read in the type of measurement device from the filename
#determine the cruise from the filename
#plot as pcolor 
# + sorting?
######################################################




import pathlib
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates
from operator import itemgetter

FOLDERNAME = "/home/ole/thesis/Preprocessing_TC_stations/PME/data"

#File with SensorID, Depth and measurement period information
sensor_positions_path = "/home/ole/thesis/Preprocessing_TC_stations/PME/pme_properties"
sensor_positions = np.genfromtxt(sensor_positions_path,skip_header= 1, usecols = (0,1,2,3,4,5,6), dtype = "str")

path = pathlib.Path(FOLDERNAME)


#splitted_foldername = FOLDERNAME.split("/")

#cruisename = splitted_foldername[5]
#flach_or_tief = splitted_foldername[8]

DATAFILENAMES = []
#go through all files of specified folder and select only files ending with "oxygen.npz"
for p in path.iterdir():
    FILENAME = str(p.parts[-1])
    splitted_filename = FILENAME.split("_")
    ending = splitted_filename[-1]
    
    if (ending == "oxygen.npz"): 
        DATAFILENAMES.append(FILENAME)


#DATAFILENAMES = sorted(DATAFILENAMES) 
print(DATAFILENAMES)

#create ouput pictures, which will be filled later in the code

#figure 1 for TC_Flach
f1, axarr1 = plt.subplots(nrows = 3, ncols = 2, sharey= "col")#, sharex = "row", )

#figure 2 for TC_Tief
#f2, axarr2 = plt.subplots(nrows = 3, ncols = 2, sharex = "row", sharey=True)

max_time_delta = dt.timedelta(minutes=1)

#loop over the files in a specified folder
for DATAFILENAME in DATAFILENAMES:
    datafile_path = FOLDERNAME+"/"+DATAFILENAME
    
    print(datafile_path)
    
    splitted_filename = DATAFILENAME.split("_")

    device = splitted_filename[0]
    cruisename = splitted_filename[1]
    flach_or_tief = splitted_filename[-2]
    print(device,cruisename,flach_or_tief,ending)

    data = np.load(datafile_path)
    print(data.files)
    
    #concentration = data["concentration"]
    saturation = data["saturation"]
    label_list = data["label_list"]
    depth = label_list[:,1].astype("float")
    utc = data["utc"]
   
    sorted_depth = sorted(depth)
    #print(depth)
    #print(sorted_depth)
    number_of_depth_bins = depth.size
    
    plot_edges_depth = np.insert(sorted_depth,0,0) #Prepend a 0 to he beginning (index 0)
    
    print(label_list)
    print(np.reshape(label_list[:,1],(number_of_depth_bins,1)))
    
    #stack the index to the label_list
    sorting_array = np.hstack((np.reshape(label_list[:,1],(number_of_depth_bins,1)),np.reshape(np.arange(number_of_depth_bins),(number_of_depth_bins,1)))).astype("float")
    
    
    print(np.shape(saturation),np.shape(label_list),np.shape(utc))
   
    """
    print("original:")
    for depth_index in range(np.shape(saturation)[0]):
        print(depth_index,depth[depth_index],np.nanmean(saturation[depth_index,:]))
    """
    
    #sort after the first column, eg the depths (sorts automatically the indices)
    #from https://stackoverflow.com/questions/2173797/how-to-sort-2d-array-by-row-in-python
    sorting_key = np.asarray(sorted(sorting_array, key= itemgetter(0)))
    
    print(sorting_array)
    print(sorting_key)
    
    sorted_saturation = []
    for index in sorting_key[:,1]:
        sorted_saturation.append(saturation[int(index),:])
    sorted_saturation = np.asarray(sorted_saturation)
    
    """
    print("sorted")
    for depth_index in range(np.shape(sorted_saturation)[0]):
        print(depth_index,sorted_depth[depth_index],np.nanmean(sorted_saturation[depth_index,:]))
    """
    
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
    img1_1 = axarr1[y_position,x_position].pcolormesh(utc,plot_edges_depth,sorted_saturation, cmap = plt.cm.RdYlBu_r, vmin = 0, vmax=50)

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

"""
f1.subplots_adjust(right=0.9)
cbar_ax = f1.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
f1.colorbar(img1_2, cax=cbar_ax).set_label("velocity [m/s]") 

f2.subplots_adjust(right=0.9)
cbar_ax = f2.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
#f2.colorbar(img2_2, cax=cbar_ax).set_label("velocity [m/s]") 
"""
        
axarr1[0,0].set_title(title_fig1)
axarr1[0,1].set_title(title_fig2)    
     
f1.suptitle("Oxygen saturation comparison")


f1.set_size_inches(18,10.5)
f1.tight_layout()
    
#Save the plot as png
plot1_name = "./pictures/"+"oxygen_comparison" 
f1.savefig(plot1_name)

 
#preferences of the x-limit
axarr1[0,0].set_xlim(left = mdates.date2num(start_flach_2017), right = mdates.date2num(start_flach_2017 + max_time_delta))
axarr1[0,1].set_xlim(left = mdates.date2num(start_tief_2017), right = mdates.date2num(start_tief_2017 + max_time_delta))
axarr1[1,0].set_xlim(left = mdates.date2num(start_flach_2018), right = mdates.date2num(start_flach_2018 + max_time_delta))
axarr1[1,1].set_xlim(left = mdates.date2num(start_tief_2018), right = mdates.date2num(start_tief_2018 + max_time_delta))
axarr1[2,0].set_xlim(left = mdates.date2num(start_flach_2019), right = mdates.date2num(start_flach_2019 + max_time_delta))
axarr1[2,1].set_xlim(left = mdates.date2num(start_tief_2019), right = mdates.date2num(start_tief_2019 + max_time_delta))


#Save the changend plot again as png
plot1_name = "./pictures/"+"oxygen_comparison_same_scale" 
f1.savefig(plot1_name)
  
plt.show()
    
       
    
    
    
    
