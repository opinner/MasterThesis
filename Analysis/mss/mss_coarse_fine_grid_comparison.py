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

longitude,distance,bathymetry,raw_Osborn,rolling_mean_Osborn,raw_Shih,rolling_mean_Shih = np.loadtxt("./data"+cruisename+"_coarse_flux_results.txt", unpack=True)
fine_longitude,fine_distance,fine_bathymetry, fine_raw_Osborn,fine_rolling_mean_Osborn,fine_raw_Shih,fine_rolling_mean_Shih = np.loadtxt("./data"+cruisename+"_fine_flux_results.txt", unpack=True)
bin_longitude,bin_distance,bin_bathymetry, bin_raw_Osborn,bin_rolling_mean_Osborn,bin_raw_Shih,bin_rolling_mean_Shih = np.loadtxt("./data"+cruisename+"_bin_flux_results.txt", unpack=True)

"""
print(len(longitude),len(distance))
print(len(fine_longitude),len(fine_distance))
print(len(bin_longitude),len(bin_distance))
for a,b in zip(longitude,bin_longitude):
    print(a,b,a==b)
"""

assert np.all(longitude == bin_longitude)

#axis.plot(longitude,raw_Osborn,"x")
#axis.plot(longitude,rolling_mean_Osborn, "k--", label = "Osborn coarse grid (eps interpolated to 02)")
axis.plot(longitude,raw_Shih,"ko", alpha = 0.2)
axis.plot(longitude,rolling_mean_Shih, "k", label = "Shih coarse grid (eps interpolated to 02)")

#axis.plot(longitude,fine_raw_Osborn,"x")
#axis.plot(fine_longitude,fine_rolling_mean_Osborn, "--", c = "tab:red", label = "Osborn fine grid (02 interpolated to eps)")
#axis.plot(longitude,fine_raw_Shih,"o")
#axis.plot(fine_longitude,fine_rolling_mean_Shih, "-", c = "tab:red", label = "Shih fine grid (02 interpolated to eps)")

#axis.plot(longitude,fine_raw_Osborn,"x")
#axis.plot(bin_longitude,bin_rolling_mean_Osborn, "--", c = "tab:blue", label = "Osborn bin grid (02 binned to eps)")
axis.plot(bin_longitude,bin_raw_Shih,"x", alpha = 0.2, c= "tab:blue")
axis.plot(bin_longitude,bin_rolling_mean_Shih, "-", c = "tab:blue", label = "Shih bin grid (02 binned to eps)")

axis.legend()
axis.set_title("Impact of coarse, fine and binned grid")
axis.set_xlabel("longitude")
axis.set_ylabel("pressure [dbar]")
output.suptitle(cruisename+": Impact of interpolation methods")

plt.show()
    

