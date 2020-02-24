import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import datetime as dt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit 
import numpy.ma as ma


#TODO change append to np.append

FOLDERNAME = "/home/ole/share-windows/adcp_slices"
path = pathlib.Path(FOLDERNAME)
DATAFILENAMES = []

#define the pictures
picture_emb169_flach_1, emb169_flach_axarr1 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb169_flach_2, emb169_flach_axarr2 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb169_flach_3, emb169_flach_axarr3 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")

picture_emb169_tief_1, emb169_tief_axarr1 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb169_tief_2, emb169_tief_axarr2 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb169_tief_3, emb169_tief_axarr3 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")

picture_emb177_flach_1, emb177_flach_axarr1 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb177_flach_2, emb177_flach_axarr2 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb177_flach_3, emb177_flach_axarr3 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")

picture_emb177_tief_1, emb177_tief_axarr1 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb177_tief_2, emb177_tief_axarr2 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb177_tief_3, emb177_tief_axarr3 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")

picture_emb217_flach_1, emb217_flach_axarr1 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb217_flach_2, emb217_flach_axarr2 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb217_flach_3, emb217_flach_axarr3 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")

picture_emb217_tief_1, emb217_tief_axarr1 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb217_tief_2, emb217_tief_axarr2 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
picture_emb217_tief_3, emb217_tief_axarr3 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")  


#go through all files of specified folder and select only files ending with .mat
for p in path.iterdir():
    all_files_name = str(p.parts[-1])

    if all_files_name[-4:] == ".mat":
        DATAFILENAMES.append(str(p.parts[-1]))

DATAFILENAMES = sorted(DATAFILENAMES) 

for DATAFILENAME in DATAFILENAMES:
    datafile_path = FOLDERNAME+"/"+DATAFILENAME
    
    #print("--------------------------------------\nfile:")
    #print(datafile_path[25:])
    #print(cruisename,flach_or_tief[3:])
    

    data = sio.loadmat(datafile_path)    
    
    print(DATAFILENAME)
    print(sio.whosmat(datafile_path))
    print(data.keys())
    
    print("\n\n\n\n")
    
    if DATAFILENAME[:4] == "ADCP":
        matlab_time = data["matlab_time"]
        utc = np.asarray(pl.num2date(matlab_time-366))
        patched_u = data["patched_u"]
        patched_v = data["patched_v"]
        depths_of_cuts = ["depths_of_cuts"]
        patched_mask = ["nan_mask"]
        
        
    elif DATAFILENAME[:3] == "PSD":
        F = data["F_1"]
        Pu = data["Pu_1"]
        Pv = data["Pv_1"]
        Ctop = data["Ctop_1"]
        Cbot = data["Cbot_1"]
        PW = data["Pw_1"]
        clockwise = data["clockwise_1"]
        anticlockwise = data["anticlockwise_1"]
        
        F = np.append(F,data["F_2"])
        Pu.append(data["Pu_2"])
        Pv.append(data["Pv_2"])
        Ctop.append(data["Ctop_2"])
        Cbot.append(data["Cbot_2"])
        Pw.append(data["Pw_2"])
        clockwise.append(data["clockwise_2"])
        anticlockwise.append(data["anticlockwise_2"])
        
        F.append(data["F_3"])
        Pu.append(data["Pu_3"])
        Pv.append(data["Pv_3"])
        Ctop.append(data["Ctop_3"])
        Cbot.append(data["Cbot_3"])
        Pw.append(data["Pw_3"])
        clockwise.append(data["clockwise_3"])
        anticlockwise.append(data["anticlockwise_3"])
        
    else:
        print("File couldn't be classifed")
