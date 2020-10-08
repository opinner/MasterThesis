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
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

import geopy.distance as geo
import gsw.conversions as gsw
import pathlib
import mss_functions as thesis
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb217","/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb169"]
LIST_OF_MSS_FOLDERS = ["/home/ole/windows/processed_mss/emb169","/home/ole/windows/processed_mss/emb177","/home/ole/windows/processed_mss/emb217"]

height_above_ground = 5 #Size of the averaging interval above ground for the BBL, has no meaning if (flux_through_halocline == True)
acceptable_slope = 2 #float('Inf') #acceptable bathymetrie difference in dbar between two neighboring data points. 

list_of_bad_profiles = []
       
for FOLDERNAME in LIST_OF_MSS_FOLDERS:
    number_of_fluxes_over_the_threshold = 0
    number_of_transects = 0
    total_number_of_fluxes = 0
    number_of_zero_flux = 0
    amount_of_missing_values = 0
    total_number_of_valid_profiles = 0    
    total_number_of_correct_profiles = 0
    
    
    
    path = pathlib.Path(FOLDERNAME)
    DATAFILENAMES = []

    splitted_foldername = FOLDERNAME.split("/")
    
    #get the cruise name from the folder name
    cruisename = splitted_foldername[-1]
    
    print(cruisename)
    
    
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
    
        #skip the short "S206" transects
        if transect_name[0:4] == "S106":
            print(transect_name,"skipped")
            list_of_bad_profiles.append(["_".join((cruisename,transect_name)),"\t\t\tshort_transect"])
            continue
            
        #something is not correct with this measurement
        if cruisename == "emb169" and transect_name[0:4] == "TS13":
            print(transect_name,"skipped, measurement looks wrong")
            list_of_bad_profiles.append(["_".join((cruisename,transect_name)),"\t\t\tmeasurement_looks_wrong"])
            continue
                
        print("\n",transect_name)
            
        
        data = np.load(datafile_path)
        
        try:
            number_of_profiles = data["number_of_profiles"] #
            lat = data["lat"] #Latitude of the profiles
            lon = data["lon"] #Longitude of the profiles
            distance = data["distance"] #distance from the starting profile (monotonically increasing)
            
            eps_pressure = data["bin_pressure"]
            eps_grid = data["bin_eps_grid"]
            eps_oxygen_sat_grid = data["bin_oxygen_sat_grid"]
            eps_oxygen_grid = data["bin_oxygen_grid"] 
            eps_pot_density_grid = data["bin_pot_density_grid"]

        
        except KeyError:
            print(transect_name," is skipped, Error during loading data")
            list_of_bad_profiles.append(["_".join((cruisename,transect_name,)),"\t\t\tdata_is_not_complete"])
            continue
        
        number_of_transects+=1 
            
        print("Number of profiles:",number_of_profiles)
        
        #print(min(eps_pressure),max(eps_pressure),len(eps_pressure))
        
        #calculate the idices of the bottom and some meters above that
        results = thesis.find_bottom_and_bottom_currents(number_of_profiles,eps_pressure,eps_pot_density_grid,eps_oxygen_grid,height_above_ground = height_above_ground)
        """
        bathymetrie                     pressure values of the first NaN value (in most cases this corresponds to the bottom, but is sometimes wrong due to missing data
        list_of_bathymetrie_indices     corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL                             pressure values of the calculated Bottom Boundary Layer (exact position depends on the criteria)
        list_of_BBL_indices             corresponding index (eg for interp_pressure or other arrays of the same size)
        BBL_range                       pressure values of "height_above_ground" meters. Follows therefore the batyhmetrie. 
        list_of_BBL_range_indices       corresponding index (eg for interp_pressure or other arrays of the same size)
        """
        bathymetrie,list_of_bathymetrie_indices = results[0]
        #BBL,list_of_BBL_indices = results[1] #not needed here
        BBL_range,list_of_BBL_range_indices = results[2]
               
        
        spread_of_profile_medians = np.nanstd(np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = 1))
        transect_median = np.nanmedian(np.log10(eps_grid[:,30:-30]),axis = None)
        outlier_count = 0

        

        #compare the pressure value of the lowest valid data point per profile th the neighbouring profiles to determine outliers
        list_of_short_profiles = thesis.get_list_of_short_profiles(number_of_profiles,bathymetrie,acceptable_slope)
        
        #hard coded to accout for one short profile, that my condtions doesn't recognize
        #if the boolean array contains at least one true
        if np.any(lon == 20.466176666666666):
            list_of_short_profiles.append(np.argmax(lon == 20.466176666666666))
        
        for profile in range(number_of_profiles):

            #if the current profile is too short, skip it
            if profile in list_of_short_profiles:
                print(str(lon[profile])+": short profile")
                list_of_bad_profiles.append(["_".join((cruisename,transect_name,str(profile))),"\t\t\tshort_profile"])
                continue
                
            
            if np.nanmean(eps_oxygen_sat_grid[profile]) < 0:
                print(cruisename,transect_name,"negative oxygen values")
                list_of_bad_profiles.append(["_".join((cruisename,transect_name,str(profile))),"\t\t\tnegative_mean_oxygen_values"])
                continue
                               
                
            #if the profile contains only nan values, profile is skipped
            if np.all(np.isnan(eps_grid[profile,:])): #from_index:to_index
                print("NaN profile")
                list_of_bad_profiles.append(["_".join((cruisename,transect_name,str(profile))),"\t\t\tprofile_contains_only_nans"])
                continue

                               
            #right now the criterion is only valid for emb217
            if cruisename == "emb217":
            #check for an outlier profile, ergo too high dissipation rates compared with the surrounding
                if np.nanmedian(np.log10(eps_grid[profile,30:-30])) > (transect_median+2*spread_of_profile_medians):      
                    #print("\toutlier")
                    outlier_count += 1
                    list_of_bad_profiles.append(["_".join((cruisename,transect_name,str(profile))),"\t\t\teps_is_too_high"])
                    continue
          
            if cruisename == "emb177":
                #index for a depth of 50db
                test_index = np.nanargmin(np.abs(eps_pressure-50))
                
                #test if the saturation at that depth is under a certain level
                if eps_oxygen_sat_grid[profile,test_index] < 50:
                    print("Oxycline drops too fast!",eps_oxygen_sat_grid[profile,test_index])
                    outlier_count += 1
                    list_of_bad_profiles.append(["_".join((cruisename,transect_name,str(profile))),"\t\t\toxygen_drops_too_fast"])
                    continue
                 
                try:    
                    #test if the saturation at that depth is under a certain level
                    if eps_pressure[np.nanargmin(thesis.central_differences(eps_oxygen_sat_grid[profile,:]))] < 50:
                        print("Oxycline is too high!",eps_pressure[np.nanargmin(thesis.central_differences(eps_oxygen_sat_grid[profile,:]))])
                        outlier_count += 1
                        #list_of_bad_profiles.append(["_".join((cruisename,transect_name,str(profile))),"oxycline_is_too_shallow"])
                        #continue
                        
                except ValueError:
                    #if np.all(np.isnan(eps_oxygen_sat_grid[profile,:])):
                    print("All NaN profile")
                    outlier_count += 1
                    list_of_bad_profiles.append(["_".join((cruisename,transect_name,str(profile))),"\t\t\tNaN_profile"])
                    continue

                                   
            total_number_of_correct_profiles+=1


    
    
np.savetxt("./data/"+"list_of_bad_profiles.txt", list_of_bad_profiles, fmt="%s", header = "profile\t\t\t\treason")    
    
    
