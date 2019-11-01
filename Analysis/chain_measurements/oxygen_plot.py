######################################################
#TODO 
#Read in the type of measurement device from the filename
#determine the cruise from the filename
#plot as pcolor 
# + sorting?
######################################################





import numpy as np
import matplotlib.pyplot as plt

LIST_OF_FOLDERS = ["/home/ole/thesis/Preprocessing_TC_stations/PME/data"]

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
        all_files_name_ending = str(p.parts[-1])
        if ((all_files_name_ending == "oxygen.npz") 
            DATAFILENAMES.append(str(p.parts[-1]))


    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    
    #create ouput pictures, which will be filled later in the code
    
    #figure 1 for the measurements
    f1, axarr1 = plt.subplots(3, sharex=True)

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
            print(datafile_path[25:])
        else:
            print("----------------------------------------------------\nfile:")
            print(datafile_path[25:])
