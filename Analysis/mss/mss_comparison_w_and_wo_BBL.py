import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt
import geopy.distance as geo
import mss_functions as thesis
from scipy import integrate as sc
import warnings
warnings.filterwarnings('ignore')
#matplotlib preferences:
MINI_SIZE = 9
SMALL_SIZE = 10.95
MEDIUM_SIZE = 12
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE, titleweight = "bold")     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MINI_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight = "bold")  # fontsize of the figure title

width = 6.2012
height = width / 1.618

output,axis = plt.subplots(1) #plt.subplots(ncols = 3, nrows = 2,sharex = True, sharey = True)
output.set_size_inches(width,height)

cruisenames = ["emb169","emb177","emb217"]
color_scheme = ['#1b9e77','#7570b3','#d95f02']

for index in range(3):
    cruisename = cruisenames[index]
    color = color_scheme[index]

    bin_longitude,bin_distance,bin_bathymetry,bin_raw_Osborn,bin_rolling_mean_Osborn,bin_raw_Shih,bin_rolling_mean_Shih = np.loadtxt("./data/"+cruisename+"_bin_flux_results.txt", unpack=True)
    bin_longitude_wob,bin_distance_wob,bin_bathymetry_wob, bin_raw_Osborn_wob,bin_rolling_mean_Osborn_wob,bin_raw_Shih_wob,bin_rolling_mean_Shih_wob = np.loadtxt("./data/"+cruisename+"_bin_flux_results_wo_BBL.txt", unpack=True)

    assert np.all(bin_longitude == bin_longitude_wob)


    #axis.plot(bin_longitude,bin_raw_Osborn,"x")
    axis.plot(bin_longitude,bin_rolling_mean_Osborn, "-", c = color, label = cruisename+" Shih with BBL")
    #axis.plot(bin_longitude,bin_raw_Shih,"o", alpha = 0.2)
    #axis.plot(bin_longitude,bin_rolling_mean_Shih, c = color, label = cruisename+" Shih with BBL")

    #axis.plot(bin_longitude,fine_raw_Osborn,"x")
    axis.plot(bin_longitude_wob,bin_rolling_mean_Osborn_wob, "--", c = color, label = cruisename+" Shih")
    #axis.plot(bin_longitude_wob,bin_raw_Shih_wob,"x", alpha = 0.2, c= color)
    #axis.plot(bin_longitude_wob,bin_rolling_mean_Shih_wob, "--", c = color, label = cruisename+" Shih")

axis.legend()
axis.set_title("Comparison of using a different flux model in the BBL")
axis.set_xlabel("longitude")
axis.set_ylabel("pressure [dbar]")
output.suptitle(cruisename+": Impact of near bottom fluxes")

plt.show()
    

