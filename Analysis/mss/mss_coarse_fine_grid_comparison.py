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

longitude,distance,raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih = np.loadtxt("./"+cruisename+"_coarse_flux_results.txt", unpack=True)
longitude,distance,fine_raw_Osborn,fine_rolling_mean_Osborn,fine_raw_Shih,fine_rolling_mean_Shih = np.loadtxt("./"+cruisename+"_fine_flux_results.txt", unpack=True)


#axis.plot(longitude,raw_Osborn,"x")
axis.plot(longitude,rolling_mean_Osborn, "k--", label = "Osborn coarse grid (eps interpolated to 02)")
#axis.plot(longitude,raw_Shih,"o")
axis.plot(longitude,rolling_mean_Shih, "k", label = "Shih coarse grid (eps interpolated to 02)")

#axis.plot(longitude,fine_raw_Osborn,"x")
axis.plot(longitude,fine_rolling_mean_Osborn, "--", c = "tab:red", label = "Osborn fine grid (02 interpolated to eps)")
#axis.plot(longitude,fine_raw_Shih,"o")
axis.plot(longitude,fine_rolling_mean_Shih, "-", c = "tab:red", label = "Shih fine grid (02 interpolated to eps)")

axis.legend()
axis.set_title("Impact of coarse and fine grid")
axis.set_xlabel("longitude")
axis.set_ylabel("pressure [dbar]")
plt.show()
    

