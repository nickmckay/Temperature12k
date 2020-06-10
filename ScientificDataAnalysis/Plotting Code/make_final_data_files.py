#==============================================================================
# This script opens the files for 12k reconstructions, restandardizes the
# values, and makes output files.
#    author: Michael P. Erb
#    date  : 4/17/2020
#==============================================================================

import numpy as np
import xarray as xr
from netCDF4 import Dataset
import os

reference_ages_12k_ens     = [0,12000]
reference_ages_12k_medians = [50,150]


### FUNCTIONS

# Calculate the area of a latitude band relative to the area of a sphere (using a sphere for Earth's shape)
def latband_weight(lat_lower,lat_upper):
    lat_lower_rad = lat_lower*(np.pi/180)
    lat_upper_rad = lat_upper*(np.pi/180)
    weight_relative = np.sin(lat_upper_rad) - np.sin(lat_lower_rad)
    weight_relative = weight_relative/2 # This makes the weight relative to the whole Earth instead of just one hemisphere
    return weight_relative

# Get the weights of the six 30 degree latitude bands
def six_latband_weights():
    lat_bounds_for_globalmean = np.array([[60,90],[30,60],[0,30],[-30,0],[-60,-30],[-90,-60]])
    lat_band_weights = np.zeros((lat_bounds_for_globalmean.shape[0])); lat_band_weights[:] = np.nan
    for i in range(lat_bounds_for_globalmean.shape[0]):
        lat_band_weights[i] = latband_weight(lat_bounds_for_globalmean[i,0],lat_bounds_for_globalmean[i,1])
    #
    return lat_band_weights


### LOAD DATA
data_dir = '/home/mpe32/analysis/14_Holocene_proxies/GMST_paper/data/'

# Load the global-means

# Load the standard calibrated composite (SCC)
gmst_scc_ensemble = np.genfromtxt(data_dir+'SCC/temp12k_SCC_globalmean.csv',delimiter=',')[:,1:]
gmst_scc_ages     = np.genfromtxt(data_dir+'SCC/temp12k_SCC_globalmean.csv',delimiter=',')[:,0]

# Load the dynamic calibrated composite (DCC)
gmst_dcc_ensemble = np.genfromtxt(data_dir+'DCC/globalMean500-6bands.csv',delimiter=',')[:,1:]
gmst_dcc_ages     = np.genfromtxt(data_dir+'DCC/globalMean500-6bands.csv',delimiter=',')[:,0]

# Load the CPS composite
gmst_cps_ensemble = np.genfromtxt(data_dir+'CPS/globalMean500-6bands-100yr2kwindow.csv',delimiter=',',skip_header=1)
gmst_cps_ages     = np.arange(0,12001,100)

# Load the paiCo composite
gmst_pai_ensemble = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled.csv',delimiter=',',skip_header=1)[:,1:],axis=0)
gmst_pai_ages     = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled.csv',delimiter=',',skip_header=1)[:,0],axis=0)


# Load the lat-band means

# Load the standard calibrated composite (SCC)
temp_scc_latbands = {}
temp_scc_latbands['60N_to_90N'] = np.genfromtxt(data_dir+'SCC/temp12k_SCC_latband_1_60N_to_90N.csv',delimiter=',')[:,1:]
temp_scc_latbands['30N_to_60N'] = np.genfromtxt(data_dir+'SCC/temp12k_SCC_latband_2_30N_to_60N.csv',delimiter=',')[:,1:]
temp_scc_latbands['0N_to_30N']  = np.genfromtxt(data_dir+'SCC/temp12k_SCC_latband_3_0N_to_30N.csv', delimiter=',')[:,1:]
temp_scc_latbands['30S_to_0S']  = np.genfromtxt(data_dir+'SCC/temp12k_SCC_latband_4_30S_to_0S.csv', delimiter=',')[:,1:]
temp_scc_latbands['60S_to_30S'] = np.genfromtxt(data_dir+'SCC/temp12k_SCC_latband_5_60S_to_30S.csv',delimiter=',')[:,1:]
temp_scc_latbands['90S_to_60S'] = np.genfromtxt(data_dir+'SCC/temp12k_SCC_latband_6_90S_to_60S.csv',delimiter=',')[:,1:]
temp_scc_latbands['ages']       = np.genfromtxt(data_dir+'SCC/temp12k_SCC_latband_1_60N_to_90N.csv',delimiter=',')[:,0]

# Load the dynamic calibrated composite (DCC)
temp_dcc_latbands = {}
temp_dcc_latbands['60N_to_90N'] = np.genfromtxt(data_dir+'DCC/60to90500-6bands.csv',  delimiter=',',skip_header=1)[:,2:]
temp_dcc_latbands['30N_to_60N'] = np.genfromtxt(data_dir+'DCC/30to60500-6bands.csv',  delimiter=',',skip_header=1)[:,2:]
temp_dcc_latbands['0N_to_30N']  = np.genfromtxt(data_dir+'DCC/0to30500-6bands.csv',   delimiter=',',skip_header=1)[:,2:]
temp_dcc_latbands['30S_to_0S']  = np.genfromtxt(data_dir+'DCC/-30to0500-6bands.csv',  delimiter=',',skip_header=1)[:,2:]
temp_dcc_latbands['60S_to_30S'] = np.genfromtxt(data_dir+'DCC/-60to-30500-6bands.csv',delimiter=',',skip_header=1)[:,2:]
temp_dcc_latbands['90S_to_60S'] = np.genfromtxt(data_dir+'DCC/-90to-60500-6bands.csv',delimiter=',',skip_header=1)[:,2:]
temp_dcc_latbands['ages']       = np.genfromtxt(data_dir+'DCC/60to90500-6bands.csv',  delimiter=',',skip_header=1)[:,1]

# Load the CPS composite
temp_cps_latbands = {}
temp_cps_latbands['60N_to_90N'] = np.genfromtxt(data_dir+'CPS/60to90500-6bands-100yr2kwindow.csv',  delimiter=',',skip_header=1)[:,2:]
temp_cps_latbands['30N_to_60N'] = np.genfromtxt(data_dir+'CPS/30to60500-6bands-100yr2kwindow.csv',  delimiter=',',skip_header=1)[:,2:]
temp_cps_latbands['0N_to_30N']  = np.genfromtxt(data_dir+'CPS/0to30500-6bands-100yr2kwindow.csv',   delimiter=',',skip_header=1)[:,2:]
temp_cps_latbands['30S_to_0S']  = np.genfromtxt(data_dir+'CPS/-30to0500-6bands-100yr2kwindow.csv',  delimiter=',',skip_header=1)[:,2:]
temp_cps_latbands['60S_to_30S'] = np.genfromtxt(data_dir+'CPS/-60to-30500-6bands-100yr2kwindow.csv',delimiter=',',skip_header=1)[:,2:]
temp_cps_latbands['90S_to_60S'] = np.genfromtxt(data_dir+'CPS/-90to-60500-6bands-100yr2kwindow.csv',delimiter=',',skip_header=1)[:,2:]
temp_cps_latbands['ages']       = np.genfromtxt(data_dir+'CPS/60to90500-6bands-100yr2kwindow.csv',  delimiter=',',skip_header=1)[:,1]

# Load the paiCo composite
temp_pai_latbands = {}
temp_pai_latbands['60N_to_90N'] = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled60to90.csv',  delimiter=',',skip_header=1)[:,1:],axis=0)
temp_pai_latbands['30N_to_60N'] = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled30to60.csv',  delimiter=',',skip_header=1)[:,1:],axis=0)
temp_pai_latbands['0N_to_30N']  = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled0to30.csv',   delimiter=',',skip_header=1)[:,1:],axis=0)
temp_pai_latbands['30S_to_0S']  = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled-30to0.csv',  delimiter=',',skip_header=1)[:,1:],axis=0)
temp_pai_latbands['60S_to_30S'] = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled-60to-30.csv',delimiter=',',skip_header=1)[:,1:],axis=0)
temp_pai_latbands['90S_to_60S'] = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled-90to-60.csv',delimiter=',',skip_header=1)[:,1:],axis=0)
temp_pai_latbands['ages']       = np.flip(np.genfromtxt(data_dir+'PAI/globalPaico500Rescaled60to90.csv',  delimiter=',',skip_header=1)[:,0], axis=0)

# Load the GAM composite
handle_gam_original = xr.open_dataset(data_dir+'GAM/annual-with-single-seas-gam-band-means.nc',decode_times=False)
temperature_gam_original         = np.flip(handle_gam_original['temperature'].values,axis=1)
temperature_samples_gam_original = np.flip(handle_gam_original['temperature_samples'].values,axis=1)
age_gam_original                 = np.flip(-1*handle_gam_original['time'].values,axis=0)
lat_gam_original                 = handle_gam_original['lat'].values
handle_gam_original.close()

# Compute global-means
lat_band_weights = six_latband_weights()
gmst_gam_ensemble = np.average(temperature_samples_gam_original,axis=0,weights=lat_band_weights)
gmst_gam_ages     = age_gam_original

# Organize latbands
temp_gam_latbands = {}
temp_gam_latbands['ages'] = age_gam_original
latbands_gam = ['90S_to_60S','60S_to_30S','30S_to_0S','0N_to_30N','30N_to_60N','60N_to_90N']
for i,latband_txt in enumerate(latbands_gam):
    print(lat_gam_original[i],latband_txt)
    temp_gam_latbands[latband_txt] = temperature_samples_gam_original[i,:,:]


# Examine the data

print('Ensemble shapes, global-means')
print('SCC: ',str(gmst_scc_ensemble.shape))
print('DCC: ',str(gmst_dcc_ensemble.shape))
print('CPS: ',str(gmst_cps_ensemble.shape))
print('PAI: ',str(gmst_pai_ensemble.shape))
print('GAM: ',str(gmst_gam_ensemble.shape))

print('Ensemble shapes, lat-band means')
print('SCC: ',str(temp_scc_latbands['60N_to_90N'].shape))
print('DCC: ',str(temp_dcc_latbands['60N_to_90N'].shape))
print('CPS: ',str(temp_cps_latbands['60N_to_90N'].shape))
print('PAI: ',str(temp_pai_latbands['60N_to_90N'].shape))
print('GAM: ',str(temp_gam_latbands['60N_to_90N'].shape))

print('Number of ages without NaNs in reconstructions, global-means')
print('SCC: ',str(np.sum(np.isfinite(np.mean(gmst_scc_ensemble,axis=1)))))
print('DCC: ',str(np.sum(np.isfinite(np.mean(gmst_dcc_ensemble,axis=1)))))
print('CPS: ',str(np.sum(np.isfinite(np.mean(gmst_cps_ensemble,axis=1)))))
print('PAI: ',str(np.sum(np.isfinite(np.mean(gmst_pai_ensemble,axis=1)))))
print('GAM: ',str(np.sum(np.isfinite(np.mean(gmst_gam_ensemble,axis=1)))))

print('Number of ages without NaNs in reconstructions, lat-band means')
print('SCC: ',str(np.sum(np.isfinite(np.mean(temp_scc_latbands['60N_to_90N'],axis=1)))))
print('DCC: ',str(np.sum(np.isfinite(np.mean(temp_dcc_latbands['60N_to_90N'],axis=1)))))
print('CPS: ',str(np.sum(np.isfinite(np.mean(temp_cps_latbands['60N_to_90N'],axis=1)))))
print('PAI: ',str(np.sum(np.isfinite(np.mean(temp_pai_latbands['60N_to_90N'],axis=1)))))
print('GAM: ',str(np.sum(np.isfinite(np.mean(temp_gam_latbands['60N_to_90N'],axis=1)))))


### CALCULATIONS

latbands_all = ['60N_to_90N','30N_to_60N','0N_to_30N','30S_to_0S','60S_to_30S','90S_to_60S']

# Check that ages are the same for each method
if ((gmst_scc_ages == gmst_dcc_ages).all() and (gmst_scc_ages == gmst_cps_ages).all() and (gmst_scc_ages == gmst_pai_ages).all() and (gmst_scc_ages == gmst_gam_ages).all() and (gmst_scc_ages == temp_scc_latbands['ages']).all() and (gmst_scc_ages == temp_dcc_latbands['ages']).all() and (gmst_scc_ages == temp_cps_latbands['ages']).all() and (gmst_scc_ages == temp_pai_latbands['ages']).all() and (gmst_scc_ages == temp_gam_latbands['ages']).all()):
    ages_12k = gmst_scc_ages
    print('Good.  Ages are the same.')
else:
    print('WARNING! AGES ARE DIFFERENT IN 12K COMPOSITES')

# Remove a reference period from each ensemble member individually
indices_ref_12k_ens = np.where((ages_12k >= reference_ages_12k_ens[0]) & (ages_12k <= reference_ages_12k_ens[1]))[0]
gmst_scc_ensemble = gmst_scc_ensemble - np.nanmean(gmst_scc_ensemble[indices_ref_12k_ens,:],axis=0)[None,:]
gmst_dcc_ensemble = gmst_dcc_ensemble - np.nanmean(gmst_dcc_ensemble[indices_ref_12k_ens,:],axis=0)[None,:]
gmst_cps_ensemble = gmst_cps_ensemble - np.nanmean(gmst_cps_ensemble[indices_ref_12k_ens,:],axis=0)[None,:]
gmst_pai_ensemble = gmst_pai_ensemble - np.nanmean(gmst_pai_ensemble[indices_ref_12k_ens,:],axis=0)[None,:]
gmst_gam_ensemble = gmst_gam_ensemble - np.nanmean(gmst_gam_ensemble[indices_ref_12k_ens,:],axis=0)[None,:]
for latband_txt in latbands_all:
    temp_scc_latbands[latband_txt] = temp_scc_latbands[latband_txt] - np.nanmean(temp_scc_latbands[latband_txt][indices_ref_12k_ens,:],axis=0)[None,:]
    temp_dcc_latbands[latband_txt] = temp_dcc_latbands[latband_txt] - np.nanmean(temp_dcc_latbands[latband_txt][indices_ref_12k_ens,:],axis=0)[None,:]
    temp_cps_latbands[latband_txt] = temp_cps_latbands[latband_txt] - np.nanmean(temp_cps_latbands[latband_txt][indices_ref_12k_ens,:],axis=0)[None,:]
    temp_pai_latbands[latband_txt] = temp_pai_latbands[latband_txt] - np.nanmean(temp_pai_latbands[latband_txt][indices_ref_12k_ens,:],axis=0)[None,:]
    temp_gam_latbands[latband_txt] = temp_gam_latbands[latband_txt] - np.nanmean(temp_gam_latbands[latband_txt][indices_ref_12k_ens,:],axis=0)[None,:]

# For each method, remove a shorter reference period from the median of all ensemble members
indices_ref_12k_medians = np.where((ages_12k >= reference_ages_12k_medians[0]) & (ages_12k <= reference_ages_12k_medians[1]))[0]
gmst_scc_ensemble = gmst_scc_ensemble - np.nanmean(np.median(gmst_scc_ensemble[indices_ref_12k_medians,:],axis=1),axis=0)
gmst_dcc_ensemble = gmst_dcc_ensemble - np.nanmean(np.median(gmst_dcc_ensemble[indices_ref_12k_medians,:],axis=1),axis=0)
gmst_cps_ensemble = gmst_cps_ensemble - np.nanmean(np.median(gmst_cps_ensemble[indices_ref_12k_medians,:],axis=1),axis=0)
gmst_pai_ensemble = gmst_pai_ensemble - np.nanmean(np.median(gmst_pai_ensemble[indices_ref_12k_medians,:],axis=1),axis=0)
gmst_gam_ensemble = gmst_gam_ensemble - np.nanmean(np.median(gmst_gam_ensemble[indices_ref_12k_medians,:],axis=1),axis=0)
for latband_txt in latbands_all:
    temp_scc_latbands[latband_txt] = temp_scc_latbands[latband_txt] - np.nanmean(np.median(temp_scc_latbands[latband_txt][indices_ref_12k_medians,:],axis=1),axis=0)
    temp_dcc_latbands[latband_txt] = temp_dcc_latbands[latband_txt] - np.nanmean(np.median(temp_dcc_latbands[latband_txt][indices_ref_12k_medians,:],axis=1),axis=0)
    temp_cps_latbands[latband_txt] = temp_cps_latbands[latband_txt] - np.nanmean(np.median(temp_cps_latbands[latband_txt][indices_ref_12k_medians,:],axis=1),axis=0)
    temp_pai_latbands[latband_txt] = temp_pai_latbands[latband_txt] - np.nanmean(np.median(temp_pai_latbands[latband_txt][indices_ref_12k_medians,:],axis=1),axis=0)
    temp_gam_latbands[latband_txt] = temp_gam_latbands[latband_txt] - np.nanmean(np.median(temp_gam_latbands[latband_txt][indices_ref_12k_medians,:],axis=1),axis=0)

# Combine latbands into arrays
def make_latband_array(latband_data):
    #
    n_time = latband_data['60N_to_90N'].shape[0]
    n_ens  = latband_data['60N_to_90N'].shape[1]
    latband_data_array = np.zeros((6,n_time,n_ens)); latband_data_array[:] = np.nan
    #
    for i,latband_txt in enumerate(latbands_all):
        print(i,latband_txt)
        latband_data_array[i,:,:] = latband_data[latband_txt]
    #
    return latband_data_array

temp_scc_latbands_array = make_latband_array(temp_scc_latbands)
temp_dcc_latbands_array = make_latband_array(temp_dcc_latbands)
temp_cps_latbands_array = make_latband_array(temp_cps_latbands)
temp_pai_latbands_array = make_latband_array(temp_pai_latbands)
temp_gam_latbands_array = make_latband_array(temp_gam_latbands)

# Combine all 12k methods into one array
gmst_all_ensemble = np.concatenate((gmst_scc_ensemble,gmst_dcc_ensemble,gmst_gam_ensemble,gmst_cps_ensemble,gmst_pai_ensemble),axis=1)
latband_all = {}
for latband_txt in latbands_all:
    latband_all[latband_txt] = np.concatenate((temp_scc_latbands[latband_txt],temp_dcc_latbands[latband_txt],temp_gam_latbands[latband_txt],temp_cps_latbands[latband_txt],temp_pai_latbands[latband_txt]),axis=1)



### OUTPUT

# If the output directory doesn't exist, create it
if not os.path.exists(data_dir+'final_data'): os.makedirs(data_dir+'final_data')

# Output all methods to a netcdf file
outputfile = Dataset(data_dir+'final_data/temp12k_alldata.nc','w',format='NETCDF4')
outputfile.createDimension('ages',len(ages_12k))
outputfile.createDimension('ens',500)
outputfile.createDimension('latbands',6)

output_age             = outputfile.createVariable('age','f8',('ages',))
output_latband_ranges  = outputfile.createVariable('latband_ranges','S10',('latbands',))
output_latband_weights = outputfile.createVariable('latband_weights','f8',('latbands',))
output_scc_globalmean  = outputfile.createVariable('scc_globalmean','f8',('ages','ens',))
output_dcc_globalmean  = outputfile.createVariable('dcc_globalmean','f8',('ages','ens',))
output_gam_globalmean  = outputfile.createVariable('gam_globalmean','f8',('ages','ens',))
output_cps_globalmean  = outputfile.createVariable('cps_globalmean','f8',('ages','ens',))
output_pai_globalmean  = outputfile.createVariable('pai_globalmean','f8',('ages','ens',))
output_scc_latbands    = outputfile.createVariable('scc_latbands','f8',('latbands','ages','ens',))
output_dcc_latbands    = outputfile.createVariable('dcc_latbands','f8',('latbands','ages','ens',))
output_gam_latbands    = outputfile.createVariable('gam_latbands','f8',('latbands','ages','ens',))
output_cps_latbands    = outputfile.createVariable('cps_latbands','f8',('latbands','ages','ens',))
output_pai_latbands    = outputfile.createVariable('pai_latbands','f8',('latbands','ages','ens',))

output_age[:]             = ages_12k
output_latband_ranges[:]  = np.array(latbands_all)
output_latband_weights[:] = lat_band_weights
output_scc_globalmean[:]  = gmst_scc_ensemble
output_dcc_globalmean[:]  = gmst_dcc_ensemble
output_gam_globalmean[:]  = gmst_gam_ensemble
output_cps_globalmean[:]  = gmst_cps_ensemble
output_pai_globalmean[:]  = gmst_pai_ensemble
output_scc_latbands[:]    = temp_scc_latbands_array
output_dcc_latbands[:]    = temp_dcc_latbands_array
output_gam_latbands[:]    = temp_gam_latbands_array
output_cps_latbands[:]    = temp_cps_latbands_array
output_pai_latbands[:]    = temp_pai_latbands_array

outputfile.close()


# Function to output data to CSV files.
def output_csv(gmst_data,latband_data,method):
    #
    # Save global-mean values
    values_to_save_gmst = np.concatenate((ages_12k[:,None],gmst_data),axis=1)
    np.savetxt(data_dir+'final_data/temp12k_'+method+'_globalmean.csv',values_to_save_gmst,delimiter=',')
    #
    for i,latband_txt in enumerate(latbands_all):
        values_to_save = np.concatenate((ages_12k[:,None],latband_data[latband_txt]),axis=1)
        np.savetxt(data_dir+'final_data/temp12k_'+method+'_latband_'+str(i+1)+'_'+latband_txt+'.csv',values_to_save,delimiter=',')
    #
    print('Outputted data for '+method)

output_csv(gmst_scc_ensemble,temp_scc_latbands,'SCC')
output_csv(gmst_dcc_ensemble,temp_dcc_latbands,'DCC')
output_csv(gmst_cps_ensemble,temp_cps_latbands,'CPS')
output_csv(gmst_pai_ensemble,temp_pai_latbands,'PAI')
output_csv(gmst_gam_ensemble,temp_gam_latbands,'GAM')


# Calculate all percentiles
global_5      = np.nanpercentile(gmst_all_ensemble,5,axis=1)
global_median = np.nanmedian(gmst_all_ensemble,axis=1)
global_95     = np.nanpercentile(gmst_all_ensemble,95,axis=1)
latband_all_5      = {}
latband_all_median = {}
latband_all_95     = {}
for latband_txt in latbands_all:
    latband_all_5[latband_txt]      = np.nanpercentile(latband_all[latband_txt],5,axis=1)
    latband_all_median[latband_txt] = np.nanmedian(latband_all[latband_txt],axis=1)
    latband_all_95[latband_txt]     = np.nanpercentile(latband_all[latband_txt],95,axis=1)

percentiles_to_save = np.concatenate((ages_12k[:,None],global_5[:,None],global_median[:,None],global_95[:,None],
                                      latband_all_5['60N_to_90N'][:,None],latband_all_median['60N_to_90N'][:,None],latband_all_95['60N_to_90N'][:,None],
                                      latband_all_5['30N_to_60N'][:,None],latband_all_median['30N_to_60N'][:,None],latband_all_95['30N_to_60N'][:,None],
                                      latband_all_5['0N_to_30N'][:,None], latband_all_median['0N_to_30N'][:,None], latband_all_95['0N_to_30N'][:,None],
                                      latband_all_5['30S_to_0S'][:,None], latband_all_median['30S_to_0S'][:,None], latband_all_95['30S_to_0S'][:,None],
                                      latband_all_5['60S_to_30S'][:,None],latband_all_median['60S_to_30S'][:,None],latband_all_95['60S_to_30S'][:,None],
                                      latband_all_5['90S_to_60S'][:,None],latband_all_median['90S_to_60S'][:,None],latband_all_95['90S_to_60S'][:,None]),axis=1)

percentiles_header = 'ages, global_5, global_median, global_95, 60N_to_90N_5, 60N_to_90N_median, 60N_to_90N_95, 30N_to_60N_5, 30N_to_60N_median, 30N_to_60N_95, 0N_to_30N_5, 0N_to_30N_median, 0N_to_30N_95, 30S_to_0S_5, 30S_to_0S_median, 30S_to_0S_95, 60S_to_30S_5, 60S_to_30S_median, 60S_to_30S_95, 90S_to_60S_5, 90S_to_60S_median, 90S_to_60S_95'

# Save a file with percentiles
np.savetxt(data_dir+'final_data/temp12k_allmethods_percentiles.csv',percentiles_to_save,delimiter=',',header=percentiles_header,comments='')


"""
# Check the new files
handle_new = xr.open_dataset(data_dir+'final_data/temp12k_alldata.nc',decode_times=False)
new_age             = handle_new['age'].values
new_latband_ranges  = handle_new['latband_ranges'].values
new_latband_weights = handle_new['latband_weights'].values
new_scc_globalmean  = handle_new['scc_globalmean'].values
new_scc_latbands    = handle_new['scc_latbands'].values
handle_new.close()

new_scc_globalmean_csv = np.loadtxt(data_dir+'final_data/temp12k_SCC_globalmean.csv',delimiter=',')
"""

