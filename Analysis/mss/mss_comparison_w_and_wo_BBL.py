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
color_scheme = ["tab:green","tab:blue","tab:red"]

for index in range(3):
    cruisename = cruisenames[index]
    color = color_scheme[index]

    bin_longitude,bin_distance,bin_bathymetry,bin_raw_Osborn,bin_rolling_mean_Osborn,bin_raw_Shih,bin_rolling_mean_Shih = np.loadtxt("./data/"+cruisename+"_bin_flux_results.txt", unpack=True)
    bin_longitude_wob,bin_distance_wob,bin_bathymetry_wob, bin_raw_Osborn_wob,bin_rolling_mean_Osborn_wob,bin_raw_Shih_wob,bin_rolling_mean_Shih_wob = np.loadtxt("./data/"+cruisename+"_bin_flux_results_wo_BBL.txt", unpack=True)

    """
    print(len(longitude),len(distance))
    print(len(fine_longitude),len(fine_distance))
    print(len(bin_longitude),len(bin_distance))
    for a,b in zip(longitude,bin_longitude):
        print(a,b,a==b)
    """

    assert np.all(bin_longitude_wob == bin_longitude)

    #axis.plot(bin_longitude,bin_raw_Osborn,"x")
    #axis.plot(bin_longitude,bin_rolling_mean_Osborn, "k--", label = "Osborn coarse grid (eps interpolated to 02)")
    #axis.plot(bin_longitude,bin_raw_Shih,"o", alpha = 0.2)
    axis.plot(bin_longitude,bin_rolling_mean_Shih, c = color, label = cruisename+" Shih")

    #axis.plot(bin_longitude,fine_raw_Osborn,"x")
    #axis.plot(bin_longitude_wob,bin_rolling_mean_Osborn_wob, "--", c = "tab:blue", label = "Osborn bin grid (02 binned to eps)")
    #axis.plot(bin_longitude_wob,bin_raw_Shih_wob,"x", alpha = 0.2, c= color)
    axis.plot(bin_longitude_wob,bin_rolling_mean_Shih_wob, "--", c = color, label = cruisename+" Shih wob")

axis.legend()
axis.set_title("Comparison with and w/o the lowermost 1.5m above the sea floor")
axis.set_xlabel("longitude")
axis.set_ylabel("pressure [dbar]")
output.suptitle(cruisename+": Impact of near bottom fluxes")

plt.show()
    

