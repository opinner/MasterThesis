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
f1, axarr1 = plt.subplots(nrows = 3, ncols = 2, sharex = "row", sharey= "col")

#figure 2 for TC_Tief
f2, axarr2 = plt.subplots(nrows = 3, ncols = 2, sharex = "row", sharey=True)

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
    print(depth)
    print(sorted_depth)
    number_of_depth_bins = depth.size
    
    #stack the index to the label_list
    sorting_array = np.hstack((label_list,np.reshape(np.arange(number_of_depth_bins),(number_of_depth_bins,1)).astype("str")))
    
    
    print(np.shape(saturation),np.shape(label_list),np.shape(utc))
   

    
    #sort after the second column, eg the depths
    #from https://stackoverflow.com/questions/2173797/how-to-sort-2d-array-by-row-in-python
    sorting_key = np.asarray(sorted(sorting_array, key=itemgetter(1)))
    
    
    sorted_saturation = []
    for index in sorting_key[:,2]:
        sorted_saturation.append(saturation[int(index),:])
    sorted_saturation = np.asarray(sorted_saturation)
    
    if flach_or_tief.lower() == "flach" or flach_or_tief.lower() == "tc-flach":
        x_position = 0
    elif flach_or_tief.lower() == "tief" or flach_or_tief.lower() == "tc-tief":
        x_position = 1
    else:
        print(flach_or_tief.lower())
        print("ERROR: Ort (flach_or_tief) nicht feststellbar")    
            
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
            
            
    #fill figure 1 with data
    img1_1 = axarr1[y_position,x_position].pcolormesh(utc,sorted_depth,sorted_saturation, cmap = plt.cm.RdYlBu_r)
    img1_2 = axarr1[y_position,x_position].pcolormesh(utc,sorted_depth,sorted_saturation, cmap = plt.cm.RdYlBu_r)
        
        #colorbar(img1_1).set_label('velocity [m/s]')
        #colorbar(img1_2).set_label('velocity [m/s]')
        
        
    #plt.pcolormesh(utc,sorted_depth,sorted_saturation)


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
        print("xlim:",row,column,axarr1[row,column].get_xlim())
        #axarr1[0,column].invert_yaxis()


        axarr2[row,column].xaxis.set_major_locator(mdates.DayLocator())
        axarr2[row,column].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
        axarr2[row,column].xaxis.set_major_formatter(hfmt)

        
        #axarr2[row,column].set_ylim(bottom = 0, top = 90)
        axarr2[row,column].invert_yaxis()





for row in range(3):
    for column in range(2):
        print("xlim:",row,column,axarr1[row,column].get_xlim())

axarr1[0,0].set_xlabel("2017")
axarr1[0,1].set_xlabel("2017")
axarr1[1,0].set_xlabel("2018")
axarr1[1,1].set_xlabel("2018")
axarr1[2,0].set_xlabel("2019")        
axarr1[2,1].set_xlabel("2019")        

axarr2[0,0].set_xlabel("2017")
axarr2[0,1].set_xlabel("2017")
axarr2[1,0].set_xlabel("2018")
axarr2[1,1].set_xlabel("2018")
axarr2[2,0].set_xlabel("2019")        
axarr2[2,1].set_xlabel("2019") 

title_fig1 = "west component u"
title_fig2 = "north component v"

f1.subplots_adjust(right=0.9)
cbar_ax = f1.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
f1.colorbar(img1_2, cax=cbar_ax).set_label("velocity [m/s]") 

f2.subplots_adjust(right=0.9)
cbar_ax = f2.add_axes([0.92, 0.15, 0.02, 0.7]) #[left, bottom, width, height] as fraction of figure size
#f2.colorbar(img2_2, cax=cbar_ax).set_label("velocity [m/s]") 
        
axarr1[0,0].set_title(title_fig1)
axarr1[0,1].set_title(title_fig2)

axarr2[0,0].set_title(title_fig1)
axarr2[0,1].set_title(title_fig2)
     
f1.suptitle("ADCP Comparison Flach")
f2.suptitle("ADCP Comparison Tief")


f1.set_size_inches(18,10.5)
f2.set_size_inches(18,10.5)
    
#Save the plot as png
plot1_name = "./pictures/"+"ADCP_comparison_flach" 
#f1.savefig(plot1_name)

plot2_name = "./pictures/"+"ADCP_comparison_tief" 
#f2.savefig(plot2_name)

"""    
#preferences of the x-limit
axarr1[0,0].set_xlim(left = mdates.date2num(start_flach_2017), right = mdates.date2num(start_flach_2017 + max_time_delta))
axarr1[0,1].set_xlim(left = mdates.date2num(start_flach_2017), right = mdates.date2num(start_flach_2017 + max_time_delta))
axarr1[1,0].set_xlim(left = mdates.date2num(start_flach_2018), right = mdates.date2num(start_flach_2018 + max_time_delta))
axarr1[1,1].set_xlim(left = mdates.date2num(start_flach_2018), right = mdates.date2num(start_flach_2018 + max_time_delta))
axarr1[2,0].set_xlim(left = mdates.date2num(start_flach_2019), right = mdates.date2num(start_flach_2019 + max_time_delta))
axarr1[2,1].set_xlim(left = mdates.date2num(start_flach_2019), right = mdates.date2num(start_flach_2019 + max_time_delta))

axarr2[0,0].set_xlim(left = mdates.date2num(start_tief_2017), right = mdates.date2num(start_tief_2017 + max_time_delta))
axarr2[0,1].set_xlim(left = mdates.date2num(start_tief_2017), right = mdates.date2num(start_tief_2017 + max_time_delta))
axarr2[1,0].set_xlim(left = mdates.date2num(start_tief_2018), right = mdates.date2num(start_tief_2018 + max_time_delta))
axarr2[1,1].set_xlim(left = mdates.date2num(start_tief_2018), right = mdates.date2num(start_tief_2018 + max_time_delta))
axarr2[2,0].set_xlim(left = mdates.date2num(start_tief_2019), right = mdates.date2num(start_tief_2019 + max_time_delta))
axarr2[2,1].set_xlim(left = mdates.date2num(start_tief_2019), right = mdates.date2num(start_tief_2019 + max_time_delta))

#Save the changend plot again as png
plot1_name = "./pictures/"+"ADCP_comparison_flach_same_scale" 
f1.savefig(plot1_name)

plot2_name = "./pictures/"+"ADCP_comparison_tief_same_scale" 
f2.savefig(plot2_name)
"""    
plt.show()
    
    
        
plt.show()    
    
    
    
    
