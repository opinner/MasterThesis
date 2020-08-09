############################################################
#this program loads all profiles from a cruise, removes outliers
#and retrieves the maximum oxyen flux in the lowermost meters 
#of the water column in choosable longitude intervals

#TODO plot mean dissipation per transect
##############################################################
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
#matplotlib preferences:
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
import datetime as dt
  
#LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217"]#,"/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb217","/home/ole/windows/processed_mss/emb169"]


search_depth_range = [52,90]
 
f,axis = plt.subplots(1) 
axis.invert_yaxis()

y,year_axis = plt.subplots(1) 
year_axis.invert_yaxis()

months = np.arange(1,13)
carstensen_results = [73.39516129032258,72.69354838709677,72.97177419354838,73.37096774193549,72.76612903225806,72.60887096774194,72.29435483870968,72.37903225806451,71.52016129032258,70.79435483870967,70.86693548387096,71.24193548387096]

year_axis.plot(months,carstensen_results,"k--", label = "Carstensen2014")
year_axis.plot(months,carstensen_results,"k.")

end_results = []
         
for FOLDERNAME in LIST_OF_MSS_FOLDERS:

    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print("----------------------------------------------")
    print(cruisename)
      
    #the pickle flag is needed to be able to load objects
    timestamps = np.load("/home/ole/Thesis/Analysis/mss/data/"+cruisename+"_mss_timestamps.npz",allow_pickle= True) 

    mss_start = timestamps[cruisename+"_mss_start"]
    mss_stop = timestamps[cruisename+"_mss_stop"]
    transect_list = timestamps[cruisename+"_transsect_list"]        
        
    cruise_start = min(mss_start) + 0.5 *( min(mss_stop) - min(mss_start))
    cruise_halocline = []
    cruise_halocline_std = []
    cruise_time_delta = []
    cruise_time = []
    
    #go through all files of specified folder and select only files ending with .mat
    for p in path.iterdir():
        all_files_name = str(p.parts[-1])

        if all_files_name[-4:] == ".npz":
            DATAFILENAMES.append(str(p.parts[-1]))

    DATAFILENAMES = sorted(DATAFILENAMES) 
    
    #print(DATAFILENAMES)
    
    
    for DATAFILENAME in DATAFILENAMES:
    
        datafile_path = FOLDERNAME+"/"+DATAFILENAME
        
        transect_name = DATAFILENAME[:-4]
    
    
        #skip the short "S106" transect
        if transect_name[0:4] == "S106":
            print(transect_name,"skipped")
            continue
            
        #something is not correct with this measurement
        if cruisename == "emb169" and transect_name[0:4] == "TS13":
            print(transect_name,"skipped, measurement looks wrong")
            continue
                
        print("\n",transect_name)
        
        #test,taxis = plt.subplots(1) 
        #taxis.invert_yaxis()    
        
        #find point in time of the transect
        for index,element in enumerate(transect_list):
            if transect_name == element:
                timestamp_index = index
                break
    
            

        transect_point_in_time =  mss_start[timestamp_index] + 0.5 * (mss_stop[timestamp_index] - mss_start[timestamp_index])
        
        transect_point_since_cruise_start = transect_point_in_time - cruise_start
        
        transect_halocline = []
        
        data = np.load(datafile_path)
        
        try:
            number_of_profiles = data["number_of_profiles"] #
            lat = data["lat"] #Latitude of the profiles
            lon = data["lon"] #Longitude of the profiles
            distance = data["distance"] #distance from the starting profile (monotonically increasing)
            
            interp_pressure = data["interp_pressure"]
            oxygen_sat_grid = data["oxygen_sat_grid"]
            salinity_grid = data["salinity_grid"]
            consv_temperature_grid = data["consv_temperature_grid"]
            density_grid = data["density_grid"]
            
            eps_pressure = data["eps_pressure"]
            eps_grid = data["eps_grid"]
            corrected_eps_grid = data["corrected_eps_grid"]
            corrected_eps_wiki_grid = data["corrected_eps_wiki_grid"]
            eps_consv_temperature_grid = data["eps_consv_temperature_grid"]
            eps_oxygen_sat_grid = data["eps_oxygen_sat_grid"]
            eps_oxygen_grid = data["eps_oxygen_grid"] 
            
            eps_N_squared_grid = data["eps_N_squared_grid"]
            eps_density_grid = data["eps_density_grid"]
            eps_pot_density_grid = data["eps_pot_density_grid"]
            #eps_viscosity_grid = data["eps_viscosity_grid"]
            eps_Reynolds_bouyancy_grid = data["eps_Reynolds_bouyancy_grid"]
            corrected_eps_Reynolds_bouyancy_grid = data["corrected_eps_Reynolds_bouyancy_grid"]
            eps_wiki_Reynolds_bouyancy_grid = data["eps_wiki_Reynolds_bouyancy_grid"]
            corrected_eps_wiki_Reynolds_bouyancy_grid = data["corrected_eps_wiki_Reynolds_bouyancy_grid"]
            
            """
            number_of_profiles              number of profiles/casts in the transect
            lat                             latitude in degrees (as a float) of the casts
            lon                             longitude in degrees (as a float) of the casts
            distance                        distance in km from the starting point of the transect
            
            interp_pressure                 equidistant 1D pressure array between the highest and the lowest measured pressure value
            oxygen_grid                     oxygen concentration in in micromol per litre as a grid (number_of_profiles x len(interp_pressure))
            salinity_grid                   salinity in g/kg as a grid (number_of_profiles x len(interp_pressure)) 
            consv_temperature_grid          conservative temperature in degrees Celsius as a grid (number_of_profiles x len(interp_pressure))
            density_grid                    density in kg/m^3 as a grid (number_of_profiles x len(interp_pressure))
            
            eps_pressure                    pressure values to the dissipation rate values (the pressure distance between points is bigger than in interp_pressure) 
            eps_grid                        measured dissipation rate values (number_of_profiles x len(eps_pressure))
            eps_consv_temperature_grid      conservative temperature as a grid (number_of_profiles x len(eps_pressure))
            eps_oxygen_grid                 oxygen concentration in micromol per litre as a grid (number_of_profiles x len(eps_pressure))
            eps_N_squared_grid              N^2, the Brunt-Vaisala frequency in 1/s^2 as a grid (number_of_profiles x len(eps_pressure))
            eps_density_grid                density in kg/m^3 as a grid (number_of_profiles x len(eps_pressure))
            
            eps_viscosity_grid
            eps_Reynolds_bouyancy_grid
            corrected_eps_Reynolds_bouyancy_grid 
            eps_wiki_Reynolds_bouyancy_grid
            corrected_eps_wiki_Reynolds_bouyancy_grid 
            """
        
        except KeyError:
            print(transect_name," is skipped, Error during loading data")
            continue
        
            
        print("Number of profiles:",number_of_profiles)
        
        #print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
        
        #calculate the idices of the bottom and some meters above that
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_density_grid,eps_oxygen_grid)
        """
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes off due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
        """
        bathymetry,list_of_bathymetrie_indices = results[0]
        #BBL,list_of_BBL_indices = results[1] #not needed here
        BBL_range,list_of_BBL_range_indices = results[2]

        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0
        

        #compare the pressure value of the lowest valid data point per profile th the neighbouring profiles to determine outliers
        list_of_short_profiles = thesis.get_list_of_short_profiles(number_of_profiles,bathymetry,acceptable_slope = 2)
        
        #hard coded to accout for one short profile, that my condtions doesn't recognize
        #if the boolean array contains at least one true
        if np.any(lon == 20.466176666666666):
            list_of_short_profiles.append(np.argmax(lon == 20.466176666666666))
        
        for profile in range(number_of_profiles):
        
            #if the current profile is too short, skip it
            if profile in list_of_short_profiles:
                print(str(lon[profile])+": short profile")
                continue
                
            
            if np.nanmean(eps_oxygen_sat_grid[profile]) < 0:
                print(cruisename,transect_name,"negative oxygen values")
                continue
                
            
            #use the indices that are {height_above_ground} meters above the bathymetry
            else:     
                from_index = int(list_of_BBL_range_indices[profile]) 
                to_index = int(list_of_bathymetrie_indices[profile])
                
                
            #if the profile contains only nan values, profile is skipped
            if np.all(np.isnan(eps_oxygen_sat_grid[profile,:])): #from_index:to_index
                print("NaN profile")
                continue

                               
            #right now the criterion is only valid for emb217
            if cruisename == "emb217":
            #check for an outlier profile, ergo too high dissipation rates compared with the surrounding
                if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+2*spread_of_profile_medians):      
                    #print("\toutlier")
                    outlier_count += 1
                    continue
          
            if cruisename == "emb177":
                #index for a depth of 50db
                test_index = np.nanargmin(np.abs(eps_pressure-50))
                #print(eps_pressure[test_index],eps_oxygen_sat_grid[profile,test_index])
                
                #test if the saturation at that depth is under a certain level
                if eps_oxygen_sat_grid[profile,test_index] < 50:
                    print("Halocline is too high!")
                    outlier_count += 1
                    continue
            
            
            #if bathymetry[profile] <= 1/2*(search_depth_range[0]+search_depth_range[0]):
            #    continue
            
            #to confine the search for the maximum gradient to the area around the halocline
            upper_boundary = np.argmin(np.abs(interp_pressure-search_depth_range[0]))
            lower_boundary = np.argmin(np.abs(interp_pressure-search_depth_range[1]))
            
            data = oxygen_sat_grid[profile,upper_boundary:lower_boundary]
            if np.count_nonzero(np.isnan(data)) >= 0.8 * len(data): 
                #print(oxygen_sat_grid[profile,upper_boundary:lower_boundary])
                #print("\nWhy the hell?\n")
                continue
            
            #print(upper_boundary,interp_pressure[upper_boundary])
            #print(lower_boundary,interp_pressure[lower_boundary])
                        
            halocline_positions = []
            
            halocline_positions.append(interp_pressure[upper_boundary+np.nanargmin(thesis.central_differences(oxygen_sat_grid[profile,upper_boundary:lower_boundary]))]) #used argmin because oxygen is decreasing with depth
            halocline_positions.append(interp_pressure[upper_boundary+np.nanargmax(thesis.central_differences(salinity_grid[profile,upper_boundary:lower_boundary]))])
            halocline_positions.append(interp_pressure[upper_boundary+np.nanargmax(thesis.central_differences(consv_temperature_grid[profile,upper_boundary:lower_boundary]))])
            #halocline_positions.append(interp_pressure[upper_boundary+np.nanargmax(thesis.central_differences(density_grid[profile,upper_boundary:lower_boundary]))])
                        
            #halocline_positions.append(interp_pressure[upper_boundary+np.argmin(np.abs(oxygen_sat_grid[profile,upper_boundary:lower_boundary])-50)])
    
            #print(np.mean(halocline_positions),np.std(halocline_positions),halocline_positions)
            
            
            #deletion = np.argmin(np.abs(halocline_positions-np.median(halocline_positions)))
            
            transect_halocline.append(np.median(halocline_positions))
            
            """
            taxis.plot(oxygen_sat_grid[profile],interp_pressure)
            taxis.plot(thesis.central_differences(oxygen_sat_grid[profile]),interp_pressure,"k")
            taxis.hlines(halocline_positions,taxis.get_xlim()[0],taxis.get_xlim()[1], alpha = 0.2)
            taxis.hlines(search_depth_range,taxis.get_xlim()[0],taxis.get_xlim()[1], "tab:red", alpha = 0.2)
            """
            #plt.show()

        mean_halocline = np.mean(transect_halocline)
        std_halocline = np.std(transect_halocline)      
    
        cruise_halocline.append(mean_halocline)
        cruise_halocline_std.append(std_halocline)
        cruise_time_delta.append(transect_point_since_cruise_start / dt.timedelta(days=1))  #convert time delta objects into a float of days
        cruise_time.append(transect_point_in_time)
    
    for color,labelname,marker,cruise in zip(["tab:red","tab:blue","tab:green"],["summer cruise emb217","winter cruise emb177","autumn emb169"],["D","x","o"],["emb217","emb177","emb169"]):
        if cruisename == cruise:
            break
    
    axis.errorbar(cruise_time_delta,cruise_halocline,cruise_halocline_std, color = color, fmt = marker, capsize = 2, label = labelname)

    float_dates = []
    for date in cruise_time:
        float_dates.append(date.month+date.day/31)

    year_axis.errorbar(cruise_start.month+cruise_start.day/31,np.mean(cruise_halocline),np.std(cruise_halocline), color = color, fmt = marker, capsize = 2, label = labelname)
    #year_axis.plot(float_dates,cruise_halocline,"-",color = color,)
    print("\n\n\n",cruisename,np.mean(cruise_halocline),"\n\n\n")
    end_results.append([cruisename,np.mean(cruise_halocline)])

print("###########################")
print(end_results)
print("###########################")


width = 8 #4.7747
height = width / 1.618

axis.set_title("Halocline depth")
axis.set_xlabel("days since cruise start")
axis.set_ylabel("pressure [dbar]")

axis.legend() 
f.set_size_inches(width, height)
f.tight_layout()

year_axis.set_title("Halocline depth")
year_axis.set_xlim(0.5,12.5)
year_axis.legend()    
year_axis.set_xlabel("months")
year_axis.set_ylabel("pressure [dbar]")

y.tight_layout()
y.set_size_inches(width, height)

f.savefig("/home/ole/Thesis/Analysis/mss/pictures/halocline_depth_days",dpi=300)
y.savefig("/home/ole/Thesis/Analysis/mss/pictures/halocline_depth_months",dpi=300)
 
plt.show()    
