import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

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

#width = 6.2012
#height = width / 1.618

#beamer figure sizes
width = 1.5*4.252 #6.2012
height = 1.5*3.7341 #* 4/3 #1.618

output,axis = plt.subplots(1) #plt.subplots(ncols = 3, nrows = 2,sharex = True, sharey = True)
output.set_size_inches(width,height)

cruisenames = ["emb169","emb177","emb217"]
color_scheme = ['#1b9e77','#7570b3','#d95f02']

for index in range(3):
    cruisename = cruisenames[index]
    color = color_scheme[index]

    bin_longitude,bin_distance,bin_bathymetry,interval_ground_distance,bin_raw_Osborn,bin_rolling_mean_Osborn,bin_raw_Shih,bin_rolling_mean_Shih = np.loadtxt("./data/"+cruisename+"_iso_flux_results.txt", unpack=True)
    bin_longitude_wob,bin_distance_wob,bin_bathymetry_wob, interval_ground_distance, bin_raw_Osborn_wob,bin_rolling_mean_Osborn_wob,bin_raw_Shih_wob,bin_rolling_mean_Shih_wob = np.loadtxt("./data/"+cruisename+"_iso_flux_results_woBBL.txt", unpack=True)

    assert np.all(bin_longitude == bin_longitude_wob)


    #axis.plot(bin_longitude,bin_raw_Osborn,"x")
    axis.plot(bin_longitude,bin_rolling_mean_Osborn, ":", c = color, label = cruisename+" Osborm with BBL")
    #axis.plot(bin_longitude,bin_raw_Shih,"o", alpha = 0.2)
    #axis.plot(bin_longitude,bin_rolling_mean_Shih, c = color, label = cruisename+" Shih with BBL")

    #axis.plot(bin_longitude,fine_raw_Osborn,"x")
    axis.plot(bin_longitude_wob,bin_rolling_mean_Osborn_wob, "--", c = color, label = cruisename+" Osborn")
    #axis.plot(bin_longitude_wob,bin_raw_Shih_wob,"x", alpha = 0.2, c= color)
    #axis.plot(bin_longitude_wob,bin_rolling_mean_Shih_wob, ":", c = color, label = cruisename+" Shih")

    if cruisename == "emb217":
        bathy_axis = axis.twinx()
        bathy_axis.set_ylim((np.nanmin(bin_bathymetry)-5,np.nanmax(bin_bathymetry)))
        bathy_axis.invert_yaxis()
        bathy_axis.set_ylabel("pressure [dbar]")
        bathy_axis.fill_between(bin_longitude,bin_bathymetry, np.ones(len(bin_longitude))*max(bin_bathymetry),color = "lightgrey", zorder = -10, alpha = 0.8, label = "bathymetry")
        # Set ax's patch invisible
        axis.patch.set_visible(False)
        # Set axtwin's patch visible and colorize it in grey
        bathy_axis.patch.set_visible(True)

        # move ax in front
        axis.set_zorder(axis.get_zorder() + 1)

        
emb169_label = mpatches.Patch(color='#1b9e77', label='autumn cruise emb169')
emb177_label = mpatches.Patch(color='#7570b3', label='winter cruise emb177')
emb217_label = mpatches.Patch(color='#d95f02', label='summer cruise emb217')
BBL_label = mlines.Line2D([], [], ls = ":", lw = 2.5, c = "k", label = "accounting for BBLs")
wo_BBL_label = mlines.Line2D([], [], ls = "--", lw = 2.5, c = "k", label = "ignoring BBLs")
bathy_label = mpatches.Patch(color='lightgrey', alpha = 0.6, label='bathymetry')

axis.legend(handles=[emb169_label,emb177_label,emb217_label,BBL_label,wo_BBL_label,bathy_label],loc = "lower right")
axis.set_title("Impact of accounting for BBLs")
axis.set_xlabel("longitude")
axis.set_xlim((20.57,20.64))
axis.set_ylabel(r"$\langle$ Osborn Oxygen flux$\rangle$ [mmol/(m$^2$d)]")
#output.suptitle(cruisename+": Impact of fluxes in BBLs")

output.subplots_adjust(top=0.934,bottom=0.11,left=0.132,hspace=0.058,wspace=0.185)
output.savefig("/home/ole/Thesis/Analysis/mss/pictures/statistics/beamer_BBL_comparison", dpi = 600)

plt.show()
    

