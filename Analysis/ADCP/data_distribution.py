import pathlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl
import datetime as dt
import matplotlib.dates as mdates
from scipy.optimize import curve_fit 

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)



#Constants
#---------------------------------------------------------------------------


#Set which fraction of missing data per depth is acceptable
acceptance_limit = 0.2

#LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_flach/adcp/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Peter_TC_tief/adcp/data","/home/ole/thesis/all_data/emb169/deployments/moorings/Lars_TC/ADCP/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/adcp/data","/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/adcp/data"]

LIST_OF_FOLDERS = ["/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data"]

for FOLDERNAME in LIST_OF_FOLDERS:
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    cruisename = splitted_foldername[5]
    flach_or_tief = splitted_foldername[8]
    
    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".mat":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    
    #create ouput pictures, which will be filled later in the code
    
    #figure 1 for the measurements
    f1, axarr1 = plt.subplots(1)

    #figure 2 for the test
    #f2, axarr2 = plt.subplots(2, sharex=True)#, sharey = True)
    
    
    for DATAFILENAME in DATAFILENAMES:
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        print("--------------------------------------\nfile:")
        print(datafile_path[25:])
        
        #print(sio.whosmat(FILENAME))
        data = sio.loadmat(datafile_path)

        
        #print(data.keys())
        data = data["adcpavg"]
        substructure = data.dtype
        #print(substructure)
        
        
        rtc = data["rtc"][0][0].flatten()
        curr = data["curr"][0][0]
        vertical_v = data["vu"][0][0].T
        
        print(np.shape(vertical_v))

        
        
        #exceptional trim for the recovery process
        if FOLDERNAME == "/home/ole/thesis/all_data/emb177/deployments/moorings/TC-tief/adcp/data":
            curr = curr[0:5488]
            rtc = rtc[0:5488]
            vertical_v = vertical_v[0:5488]
            
            
        if FOLDERNAME == "/home/ole/thesis/all_data/emb177/deployments/moorings/TC-flach/adcp/data":
            curr = curr[0:12775]
            rtc = rtc[0:12775]
            vertical_v = vertical_v[0:12775]
            
        west_east = np.real(curr).T
        north_south = np.imag(curr).T

        #convert matlab time to utc
        utc = np.asarray(pl.num2date(rtc-366))
        
        
        axarr1.hist(vertical_v.flatten(), bins = 150)
        print(np.nanstd(vertical_v.flatten()))
        plt.show()

