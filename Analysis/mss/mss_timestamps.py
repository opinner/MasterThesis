import scipy.io as sio
import datetime as dt
import pathlib
import numpy as np

#contains the MSS Data
LIST_OF_MSS_FOLDERS = ["/home/ole/share-windows/emb217_mss_data","/home/ole/share-windows/emb177_mss_data/","/home/ole/share-windows/emb169_mss_data/MSS055/matlab/","/home/ole/share-windows/emb169_mss_data/MSS038/matlab/"]

emb169_mss_start = []
emb169_mss_stop = []
emb169_transsect_list = []

emb177_mss_start = []
emb177_mss_stop = []
emb177_transsect_list = []

emb217_mss_start = []
emb217_mss_stop = [] 
emb217_transsect_list = []


 
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
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        if DATAFILENAME == "TS11_TODL_merged.mat":
            continue

        
        #print(sio.whosmat(FILENAME))
        data = sio.loadmat(datafile_path) 

        STA_substructure = data["STA"]
        #print("\nSTA\n",STA_substructure.dtype)
        
        date = STA_substructure["date"][0]
 
        #print(date[0][0])
        #print(date[-1][0])
        
        start_time = dt.datetime.strptime(date[0][0],"%d-%b-%Y %X")
        stop_time = dt.datetime.strptime(date[-1][0],"%d-%b-%Y %X")
        
        print(cruisename)
        print(DATAFILENAME)
        print(start_time)
        print(stop_time)
        
        
        if cruisename == "emb169":
            emb169_mss_start.append(start_time)
            emb169_mss_stop.append(stop_time)
            emb169_transsect_list.append(DATAFILENAME)
            
        if cruisename == "emb177":
            emb177_mss_start.append(start_time)
            emb177_mss_stop.append(stop_time)
            emb177_transsect_list.append(DATAFILENAME)
            
        if cruisename == "emb217":
            emb217_mss_start.append(start_time)
            emb217_mss_stop.append(stop_time)
            emb217_transsect_list.append(DATAFILENAME)
            
#TODO: Test if names lists are unique

print(emb169_mss_start)
print(emb177_mss_start)
print(emb217_mss_start)
            
np.savez("./emb169_mss_timestamps",emb169_mss_start = emb169_mss_start,  emb169_mss_stop = emb169_mss_stop, emb169_transsect_list = emb169_transsect_list)
np.savez("./emb177_mss_timestamps",emb177_mss_start = emb177_mss_start,  emb177_mss_stop = emb177_mss_stop, emb177_transsect_list = emb177_transsect_list)            
np.savez("./emb217_mss_timestamps",emb217_mss_start = emb217_mss_start,  emb217_mss_stop = emb217_mss_stop, emb217_transsect_list = emb217_transsect_list)
            
            
            
            
            
            
            
            
            
            
            
            
