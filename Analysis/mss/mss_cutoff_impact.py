import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt
import geopy.distance as geo
import mss_functions as thesis
from scipy import integrate as sc
import warnings
warnings.filterwarnings('ignore')

output,axis = plt.subplots(1) #plt.subplots(ncols = 3, nrows = 2,sharex = True, sharey = True)

cruisenames = ["emb169","emb177","emb217"]
cruisename = cruisenames[0]

longitude,distance,raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih = np.loadtxt("./"+cruisename+"_flux_results_wo_cutoff.txt", unpack=True)
longitude,distance,cutoff_raw_Osborn,cutoff_rolling_mean_Osborn,cutoff_raw_Shih,cutoff_rolling_mean_Shih = np.loadtxt("./"+cruisename+"_flux_results_cutoff.txt", unpack=True) 


#axis.plot(longitude,raw_Osborn,"x")
axis.plot(longitude,rolling_mean_Osborn, "k--", label = "no cut off (Osborn coarse grid)")
#axis.plot(longitude,raw_Shih,"o")
axis.plot(longitude,rolling_mean_Shih, "k", label = "no cut off (Shih coarse grid)")

#axis.plot(longitude,cutoff_raw_Osborn,"x")
axis.plot(longitude,cutoff_rolling_mean_Osborn, "--", c = "tab:red", label = "cut off at -500 mmol/m²/d (Osborn coarse grid)")
#axis.plot(longitude,cutoff_raw_Shih,"o")
axis.plot(longitude,cutoff_rolling_mean_Shih, "-", c = "tab:red", label = "cut off at -500 mmol/m²/d (Shih coarse grid)")

axis.legend()
axis.set_title("Impact of a cutoff value for oxygen fluxes")
axis.set_xlabel("longitude")
axis.set_ylabel("pressure [dbar]")
plt.show()
    

