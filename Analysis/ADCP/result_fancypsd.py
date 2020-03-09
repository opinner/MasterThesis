#----------------------------------------------------------------
# plots the 3 horizontal adcp measured velocities timeseries (per mooring)
#plus their power spectral density and rotary spectra

#----------------------------------------------------------------
#TODO What is wrong with emb217_flach?

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

FOLDERNAME = "/home/ole/share-windows/adcp_slices"
path = pathlib.Path(FOLDERNAME)
DATAFILENAMES = []


#go through all files of specified folder and select only files ending with .mat
for p in path.iterdir():
    all_files_name = str(p.parts[-1])

    if all_files_name[-4:] == ".mat":
        DATAFILENAMES.append(str(p.parts[-1]))

DATAFILENAMES = sorted(DATAFILENAMES) 


mooring_names = ["emb169_flach","emb169_tief","emb177_flach","emb177_tief","emb217_tief"]#,"emb217_flach"] #["emb217_flach"] #

for mooring in  mooring_names:

    f1, axarr1 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
    f2, axarr2 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
    f3, axarr3 = plt.subplots(nrows = 2, ncols = 2, sharex = "row", sharey = "row")
    for DATAFILENAME in DATAFILENAMES:


        datafile_path = FOLDERNAME+"/"+DATAFILENAME        
        if DATAFILENAME[:4] == "ADCP":
            current_mooring = DATAFILENAME[22:-4]
            
        elif DATAFILENAME[:3] == "PSD":
            current_mooring = DATAFILENAME[4:-4]
            
        #print(current_mooring,mooring)  
        if current_mooring != mooring:
            continue
        
        data = sio.loadmat(datafile_path)    
        
        print(DATAFILENAME)
        #print(sio.whosmat(datafile_path))
        #print(data.keys())
        

        if DATAFILENAME[:4] == "ADCP":
            matlab_time = data["matlab_time"]
            
            if np.shape(matlab_time)[0] == 1:
                matlab_time = matlab_time[0]
            
            """
            print(matlab_time)
            print(np.shape(matlab_time))
            print(np.shape(matlab_time[0]))
            print(type(matlab_time))
            print(type(matlab_time[0][0]))
            """
        
    
            days = []
            for i in range(3):
                days.append(matlab_time[i].flatten()-matlab_time[i][0])
                
            utc = []
            for i in range(3):
                utc.append(np.asarray(pl.num2date(matlab_time[i]-366)).flatten())

            
            patched_u = data["patched_u"]
            patched_v = data["patched_v"]
            depths_of_cuts = data["depths_of_cuts"].flatten()
            patched_mask = data["nan_mask"]
            
            if np.shape(patched_u)[0] == 1:
                patched_u = patched_u[0]
                patched_v = patched_v[0]
                patched_mask = patched_mask[0]            
            #print(np.shape(utc))
            #print(np.shape(patched_u[0]))
            
            #top=0.954,bottom=0.081,left=0.05,right=0.985,hspace=0.231,wspace=0.02
            print("test",np.shape(days[0]),np.shape(patched_u[0]))
        
            linewidth = 0.5
            axarr1[0,0].plot(days[0],patched_u[0].flatten(),"b",label = "u",linewidth=linewidth)                                          
            axarr1[0,1].plot(days[0],patched_v[0].flatten(),"r",label = "v",linewidth=linewidth)
            nan_u  = ma.masked_where(patched_mask[0],patched_u[0]).flatten()
            nan_v = ma.masked_where(patched_mask[0],patched_v[0]).flatten()
            axarr1[0,0].plot(days[0],nan_u,"k",linewidth=linewidth)
            axarr1[0,1].plot(days[0],nan_v,"k",linewidth=linewidth)
            
            axarr1[0,0].set_title(current_mooring+" cut at "+str(np.round(depths_of_cuts[0],2))+" dbar")
            #f1.suptitle(current_mooring+" cut at "+str(depths_of_cuts[0])+" dbar")
            
            axarr2[0,0].plot(days[1],patched_u[1].flatten(),"b",label = "u",linewidth=linewidth)                                           
            axarr2[0,1].plot(days[1],patched_v[1].flatten(),"r",label = "v",linewidth=linewidth)
            nan_u  = ma.masked_where(patched_mask[1],patched_u[1]).flatten()
            nan_v = ma.masked_where(patched_mask[1],patched_v[1]).flatten()
            axarr2[0,0].plot(days[1],nan_u,"k",linewidth=linewidth)
            axarr2[0,1].plot(days[1],nan_v,"k",linewidth=linewidth)
            axarr2[0,0].set_title(current_mooring+" cut at "+str(np.round(depths_of_cuts[1],2))+" dbar")
            #f2.suptitle(current_mooring+" cut at "+str(depths_of_cuts[1])+" dbar")
                        
            axarr3[0,0].plot(days[2],patched_u[2].flatten(),"b",label = "u",linewidth=linewidth)                                               
            axarr3[0,1].plot(days[2],patched_v[2].flatten(),"r",label = "v",linewidth=linewidth)
            nan_u  = ma.masked_where(patched_mask[2],patched_u[2]).flatten()
            nan_v = ma.masked_where(patched_mask[2],patched_v[2]).flatten()
            axarr3[0,0].plot(days[2],nan_u,"k",linewidth=linewidth)
            axarr3[0,1].plot(days[2],nan_v,"k",linewidth=linewidth)       
            axarr3[0,0].set_title(current_mooring+" cut at "+str(np.round(depths_of_cuts[2],2))+" dbar")  
            #f3.suptitle(current_mooring+" cut at "+str(depths_of_cuts[2])+" dbar")
                        
        elif DATAFILENAME[:3] == "PSD":
            print("Found PSD")
            F = [data["F_1"].flatten()]
            Pu = [data["Pu_1"].flatten()]
            Pv = [data["Pv_1"].flatten()]
            Ctop = [data["Ctop_1"].flatten()]
            Cbot = [data["Cbot_1"].flatten()]
            Pw = [data["Pw_1"].flatten()]
            clockwise = [data["clockwise_1"].flatten()]
            anticlockwise = [data["anticlockwise_1"].flatten()]
                        
            F.append(data["F_2"].flatten())
            Pu.append(data["Pu_2"].flatten())
            Pv.append(data["Pv_2"].flatten())
            Ctop.append(data["Ctop_2"].flatten())
            Cbot.append(data["Cbot_2"].flatten())
            Pw.append(data["Pw_2"].flatten())
            clockwise.append(data["clockwise_2"].flatten())
            anticlockwise.append(data["anticlockwise_2"].flatten())
                        
            F.append(data["F_3"].flatten())
            Pu.append(data["Pu_3"].flatten())
            Pv.append(data["Pv_3"].flatten())
            Ctop.append(data["Ctop_3"].flatten())
            Cbot.append(data["Cbot_3"].flatten())
            Pw.append(data["Pw_3"].flatten())
            clockwise.append(data["clockwise_3"].flatten())
            anticlockwise.append(data["anticlockwise_3"].flatten())
            
            axarr1[1,0].loglog(F[0],Pu[0],"b", label = "u")                                               
            axarr1[1,0].loglog(F[0],Pv[0],"r", label = "v")
            axarr1[1,1].loglog(F[0],clockwise[0],label = "clockwise")                                               
            axarr1[1,1].loglog(F[0],anticlockwise[0], label = "anticlockwise")
            axarr1[1,0].loglog(F[0],Ctop[0],"k", alpha = 0.5)            
            axarr1[1,0].loglog(F[0],Cbot[0],"k", alpha = 0.5)      
            axarr1[1,1].loglog(F[0],Ctop[0],"k", alpha = 0.5)            
            axarr1[1,1].loglog(F[0],Cbot[0],"k", alpha = 0.5)                   
             
            print(np.shape(F[0]),np.shape(clockwise[0])) 
            print(np.shape(F[1]),np.shape(clockwise[1]))
            print(np.shape(F[2]),np.shape(clockwise[2]))
                                   
            axarr2[1,0].loglog(F[1],Pu[1],"b",label = "u")                                               
            axarr2[1,0].loglog(F[1],Pv[1],"r",label = "v")
            axarr2[1,1].loglog(F[1],clockwise[1],label = "clockwise")                                               
            axarr2[1,1].loglog(F[1],anticlockwise[1], label = "anticlockwise")
            axarr2[1,0].loglog(F[1],Ctop[1],"k", alpha = 0.5)            
            axarr2[1,0].loglog(F[1],Cbot[1],"k", alpha = 0.5)    
            axarr2[1,1].loglog(F[1],Ctop[1],"k", alpha = 0.5)            
            axarr2[1,1].loglog(F[1],Cbot[1],"k", alpha = 0.5)    
                                   
            axarr3[1,0].loglog(F[2],Pu[2],"b",label = "u")                                               
            axarr3[1,0].loglog(F[2],Pv[2],"r",label = "v")
            axarr3[1,1].loglog(F[2],clockwise[2],label = "clockwise")                                               
            axarr3[1,1].loglog(F[2],anticlockwise[2], label = "anticlockwise")
            axarr3[1,0].loglog(F[2],Ctop[2],"k", alpha = 0.5)            
            axarr3[1,0].loglog(F[2],Cbot[2],"k", alpha = 0.5)    
            axarr3[1,1].loglog(F[2],Ctop[2],"k", alpha = 0.5)            
            axarr3[1,1].loglog(F[2],Cbot[2],"k", alpha = 0.5)    
                           
        else:
            print("File couldn't be classifed")

    for row in range(2):
        for column in range(2):
            axarr1[row,column].legend()
            axarr2[row,column].legend()
            axarr3[row,column].legend()   
    
    axarr1[0,0].set_xlabel("days")
    axarr1[0,1].set_xlabel("days")
    axarr2[0,0].set_xlabel("days")
    axarr2[0,1].set_xlabel("days")
    axarr3[0,0].set_xlabel("days")
    axarr3[0,1].set_xlabel("days")                   

    axarr1[0,1].set_ylim(-0.2,0.2)
    axarr1[0,0].set_ylim(-0.2,0.2)
    axarr2[0,1].set_ylim(-0.2,0.2)
    axarr2[0,0].set_ylim(-0.2,0.2)    
    axarr3[0,1].set_ylim(-0.2,0.2)
    axarr3[0,0].set_ylim(-0.2,0.2)

    axarr1[1,0].set_xlabel("cycles per hour")
    axarr1[1,1].set_xlabel("cycles per hour")
    axarr2[1,0].set_xlabel("cycles per hour")
    axarr2[1,1].set_xlabel("cycles per hour")
    axarr3[1,0].set_xlabel("cycles per hour")
    axarr3[1,1].set_xlabel("cycles per hour")   
    
    
    f1.set_size_inches(1.618*7.2,7.2)
    f2.set_size_inches(1.618*7.2,7.2)
    f3.set_size_inches(1.618*7.2,7.2)
            
    plot_name = "/home/ole/share-windows/adcp_slices/"+mooring
    f1.savefig(plot_name+str(1))
    f2.savefig(plot_name+str(2))    
    f3.savefig(plot_name+str(3)) 
    plt.close(fig = "all")
        
    print("\n")
    
plt.show()
