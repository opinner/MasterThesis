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
 
        if (datafile_path == "/home/ole/share-windows/emb217_mss_data/TR1-8.mat"):
            date = np.delete(date,np.s_[33:47])
    
        if cruisename == "emb169" and  DATAFILENAME[:-4] == "TS118":
            date = date[:21]
            
        if cruisename == "emb169" and  DATAFILENAME[:-4] == "TRR109":
            date = date[:-1]
                
        start_time = dt.datetime.strptime(date[0][0],"%d-%b-%Y %X")
        stop_time = dt.datetime.strptime(date[-1][0],"%d-%b-%Y %X")
        
        print(cruisename,DATAFILENAME[:-4], start_time, stop_time)
        
        #if DATAFILENAME[:-4] in ["TS118","TS119"]:
        #    for i in range(date.size):
        #        print(DATAFILENAME[:-4],dt.datetime.strptime(date[i][0],"%d-%b-%Y %X"))
        
        if cruisename == "emb169":
            emb169_mss_start.append(start_time)
            emb169_mss_stop.append(stop_time)
            emb169_transsect_list.append(DATAFILENAME[:-4])
            
        if cruisename == "emb177":
            emb177_mss_start.append(start_time)
            emb177_mss_stop.append(stop_time)
            emb177_transsect_list.append(DATAFILENAME[:-4])
            
        if cruisename == "emb217":
            emb217_mss_start.append(start_time)
            emb217_mss_stop.append(stop_time)
            emb217_transsect_list.append(DATAFILENAME[:-4])
            
#TODO: Test if names lists are unique
assert(len(emb169_transsect_list) == len(set(emb169_transsect_list)))
assert(len(emb177_transsect_list) == len(set(emb177_transsect_list)))
assert(len(emb217_transsect_list) == len(set(emb217_transsect_list)))

            
np.savez("./emb169_mss_timestamps",emb169_mss_start = emb169_mss_start,  emb169_mss_stop = emb169_mss_stop, emb169_transsect_list = emb169_transsect_list)
np.savez("./emb177_mss_timestamps",emb177_mss_start = emb177_mss_start,  emb177_mss_stop = emb177_mss_stop, emb177_transsect_list = emb177_transsect_list)            
np.savez("./emb217_mss_timestamps",emb217_mss_start = emb217_mss_start,  emb217_mss_stop = emb217_mss_stop, emb217_transsect_list = emb217_transsect_list)
            
            
            
            
            
            
            
            
            
            
            
            
