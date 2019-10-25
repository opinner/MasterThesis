import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt
import matplotlib.dates as mdates
import scipy.io as sio
import pylab as pl

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #cax.set_label('Temperature / Celsius')
    return fig.colorbar(mappable, cax=cax)
    
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
        
data = sio.loadmat("/home/ole/thesis/all_data/emb217/deployments/moorings/TC_Flach/ADCP600/data/EMB217_TC-flach_adcp600_val.mat")

#print(data.keys())
data = data["adcpavg"]
substructure = data.dtype
#print(substructure)

rtc = data["rtc"][0][0].flatten()
curr = data["curr"][0][0]
vertical_v = data["vu"][0][0].T
#convert matlab time to utc
utc = np.asarray(pl.num2date(rtc-366))        
        
npzfile = np.load("dissipation_rate_estimation.npz")
print(npzfile.files)
depth = npzfile["depth"]
utc_chunks = npzfile["utc"]
dissipation_rate = npzfile["dissipation_rate"]


#figure 1 for the measurements
f1, axarr1 = plt.subplots(2, sharex=True, sharey = True)

#figure 2 for the test
f2, axarr2 = plt.subplots(2, sharex=True, sharey = True)

#figure 3 for the test
f3, axarr3 = plt.subplots(ncols = 2)

#figure 4 for the test
#f4, axarr4 = plt.subplots(2, sharex=True, sharey = True)

#set plot limit to either zero or the minimum depth
bottom_limit = max(0,min(depth))


#Einstellungen Plot:
vmin = -0.3
vmax = +0.3     

#fill figure 1 with data
img1_1 = axarr1[0].pcolormesh(utc,depth,vertical_v, vmin = -0.0075, vmax = 0.0075, cmap = plt.cm.RdYlBu_r)
img1_2 = axarr1[1].pcolormesh(utc_chunks,depth, np.log10(dissipation_rate), cmap = plt.cm.RdYlBu_r)

#img2_1 = axarr2[0].pcolormesh(utc_chunks,depth,dissipation_rate, cmap = plt.cm.RdYlBu_r)
#img2_2 = axarr2[1].pcolormesh(utc_chunks,depth,np.log10(dissipation_rate), cmap = plt.cm.RdYlBu_r)
img2_1 = axarr2[0].pcolormesh(utc_chunks,depth, np.log10(dissipation_rate), vmin = -23, vmax = -15, cmap = plt.cm.RdYlBu_r)
img2_2 = axarr2[1].pcolormesh(utc_chunks,depth, np.log10(dissipation_rate), vmin = -15, vmax = -6.5, cmap = plt.cm.RdYlBu_r)

img3_1 = axarr3[0].hist(dissipation_rate.flatten(), bins = 1500)
img3_2 = axarr3[1].hist(np.log10(dissipation_rate.flatten()), bins = 150)      


#img4_1 = axarr4[0].pcolormesh(utc_chunks,depth, np.log10(dissipation_rate), vmin = -23, vmax = -6.5, cmap = plt.cm.seismic)
#img4_2 = axarr4[1].pcolormesh(utc_chunks,depth, np.log10(dissipation_rate), vmin = -23, vmax = -6.5, cmap = shiftedColorMap(plt.cm.seismic, midpoint = 1 - 23 / (23 + abs(-6.5))))

1 - vmax / (vmax + abs(vmin))
plt.cm.seismic

print(np.nanmean(dissipation_rate.flatten()),np.nanstd(dissipation_rate.flatten()),np.nanmax(dissipation_rate.flatten()),np.nanmin(dissipation_rate.flatten()))

print(np.log10(np.nanmean(dissipation_rate.flatten())),np.log10(np.nanstd(dissipation_rate.flatten())),np.log10(np.nanmax(dissipation_rate.flatten())),np.log10(np.nanmin(dissipation_rate.flatten())))

print(np.nanmean(np.log10(dissipation_rate.flatten())),np.nanstd(np.log10(dissipation_rate.flatten())),np.nanmax(np.log10(dissipation_rate.flatten())),np.nanmin(np.log10(dissipation_rate.flatten())))
        
#Preferences of the plots

#figure 1
#set the xticks (dateformat and tick location)
hfmt = mdates.DateFormatter('%d %b')
axarr1[1].xaxis.set_major_locator(mdates.DayLocator())
axarr1[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
axarr1[1].xaxis.set_major_formatter(hfmt)

axarr1[1].set_ylim(bottom = bottom_limit)
axarr1[1].invert_yaxis()

axarr2[1].xaxis.set_major_locator(mdates.DayLocator())
axarr2[1].xaxis.set_minor_locator(mdates.HourLocator(byhour = [0,6,12,18],interval = 1))
axarr2[1].xaxis.set_major_formatter(hfmt)

axarr2[1].set_ylim(bottom = bottom_limit)
axarr2[1].invert_yaxis()

#axarr4[1].invert_yaxis()

#title_fig1_1 = "adcp "+cruisename+" "+flach_or_tief+" east component"
#title_fig1_2 = "adcp "+cruisename+" "+flach_or_tief+" north component"

#axarr1[0].set_title(title_fig1_1)
#axarr1[1].set_title(title_fig1_2)

colorbar(img1_1).set_label('velocity [m/s]')
colorbar(img1_2).set_label('log10(dissipation rate)')

colorbar(img2_1).set_label('dissipation rate')
colorbar(img2_2).set_label('log10(dissipation rate)')

colorbar(img4_1).set_label('not shifted')
colorbar(img4_2).set_label('shifted')

f1.set_size_inches(12,7)
f2.set_size_inches(12,7)
f3.set_size_inches(12,7)
f4.set_size_inches(12,7)

#Save the plot as png
#plot1_name = "./pictures/"+"adcp_"+cruisename+"_"+flach_or_tief 
#f1.savefig(plot1_name)


plt.show()
