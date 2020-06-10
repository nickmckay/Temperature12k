#=============================================================================
# This script makes some figures of the Temp12k GMST data, as well as a
# comparison between different Holocene composites.
#    author: Michael P. Erb
#    date  : 6/8/2020
#=============================================================================

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import copy
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from mpl_toolkits.basemap import Basemap

save_instead_of_plot = True
reference_period     = [1800,1899]

### LOAD DATA
data_dir = '/home/mpe32/analysis/14_Holocene_proxies/GMST_paper/data/'

# Load 12k data
handle = xr.open_dataset(data_dir+'final_data/temp12k_alldata.nc',decode_times=False)
gmst_12k_ages     = handle['age'].values
gmst_scc_ensemble = handle['scc_globalmean'].values
gmst_dcc_ensemble = handle['dcc_globalmean'].values
gmst_cps_ensemble = handle['cps_globalmean'].values
gmst_pai_ensemble = handle['pai_globalmean'].values
gmst_gam_ensemble = handle['gam_globalmean'].values
handle.close()

# Load 2k GMST data
data_dir_2k = '/projects/pd_lab/data/paleoclimate_reconstructions/GMST_2k_reconstructions/final_recon_results/recons/'
filenames_2k = ['BHM','CPS','DA','M08','OIE','PAI','PCR']

years_2k = np.genfromtxt(data_dir_2k+'BHM.txt',skip_header=1)[:,0]
ages_2k = 1950-years_2k

nfiles = len(filenames_2k)
ntime = 2000
niter = 1000
data_gmst_all = np.zeros((nfiles,ntime,niter)); data_gmst_all[:] = np.nan
for i,filename in enumerate(filenames_2k):
    print(filename)
    data_gmst_all[i,:,:] = np.genfromtxt(data_dir_2k+filename+'.txt',skip_header=1)[:,1:]

# Load the ERA-20C reanalysis data
data_dir_era20c = '/projects/pd_lab/data/modern_datasets/ERA20C/'
handle_era20c = xr.open_dataset(data_dir_era20c+'t2m_annual_1900_to_2010_era20c.nc',decode_times=False)
tas_era20c   = handle_era20c['t2m'].values
lon_era20c   = handle_era20c['lon'].values
lat_era20c   = handle_era20c['lat'].values
years_era20c = handle_era20c['years'].values
handle_era20c.close()
ages_era20c = 1950-years_era20c

# Load other Holocene reconstructions
data_dir_recons = '/projects/pd_lab/data/paleoclimate_reconstructions/Holocene_reconstructions/'

data_shakun = pd.ExcelFile(data_dir_recons+'Shakun_etal_2012/41586_2012_BFnature10915_MOESM60_ESM.xls').parse('TEMPERATURE STACKS').values
ages_shakun       = (data_shakun[:,0]*1000).astype(np.float)
globalmean_shakun = data_shakun[:,1].astype(np.float)
onesigma_shakun   = data_shakun[:,2].astype(np.float)

data_marcott = pd.ExcelFile(data_dir_recons+'Marcott_etal_2013/Marcott.SM.database.S1.xlsx').parse('TEMPERATURE STACKS').values
ages_marcott       = data_marcott[5:,2].astype(np.float)
globalmean_marcott = data_marcott[5:,3].astype(np.float)
onesigma_marcott   = data_marcott[5:,4].astype(np.float)

# Load the proxy locations

# Load the Temp12k proxy metadata
data_dir_proxies_temp12k = '/projects/pd_lab/data/proxies/Holocene/database_v1/'
metadata_temp12k = pd.ExcelFile(data_dir_proxies_temp12k+'Temp12k_v1_record_list.xlsx').parse('Suppl_Table_2').values
sitenames_all_temp12k = metadata_temp12k[2:798,0]
lat_all_temp12k       = metadata_temp12k[2:798,1].astype(np.float)
lon_all_temp12k       = metadata_temp12k[2:798,2].astype(np.float)

# Get the proxy lat and lon metadata for the other Holocene reconstructions
data_dir_proxies_recons = '/projects/pd_lab/data/paleoclimate_reconstructions/Holocene_reconstructions/'
metadata_shakun = pd.ExcelFile(data_dir_proxies_recons+'Shakun_etal_2012/41586_2012_BFnature10915_MOESM60_ESM.xls').parse('METADATA').values
sitenames_all_shakun = metadata_shakun[:,1]
lat_all_shakun       = metadata_shakun[:,4].astype(np.float)
lon_all_shakun       = metadata_shakun[:,5].astype(np.float)

metadata_marcott = pd.ExcelFile(data_dir_proxies_recons+'Marcott_etal_2013/Marcott.SM.database.S1.xlsx').parse('METADATA').values
sitenames_all_marcott = metadata_marcott[2:78,1]
lat_all_marcott       = metadata_marcott[2:78,4]
lon_all_marcott       = metadata_marcott[2:78,5]

# One of Marcott's proxies uses two locations.  Plot both locations
index_proxy = np.where(lat_all_marcott == '71.3/ 81.0')[0]
lat_all_marcott[index_proxy] = 71.3
lon_all_marcott[index_proxy] = 26.7
lat_all_marcott = np.insert(lat_all_marcott,index_proxy+1,81.0)
lon_all_marcott = np.insert(lon_all_marcott,index_proxy+1,-71)
sitenames_all_marcott = np.insert(sitenames_all_marcott,index_proxy+1,'Agassiz & Renland, other location')
lat_all_marcott = lat_all_marcott.astype(np.float)
lon_all_marcott = lon_all_marcott.astype(np.float)

# A function to remove multiple proxy entrees at the same location
def proxy_locations(sitenames_all,lat_all,lon_all):
    #
    n_rows = len(sitenames_all)
    sitenames_chosen = []
    for i in range(n_rows):
        if (sitenames_all[i] not in sitenames_chosen) or (sitenames_all[i] == '-'):
            sitenames_chosen.append(sitenames_all[i])
        else:
            lat_all[i] = np.nan
            lon_all[i] = np.nan
    #
    # Get rid of nans
    valid_data = np.isfinite(lat_all) & np.isfinite(lon_all)
    sitenames_all = sitenames_all[valid_data]
    lat_all       = lat_all[valid_data]
    lon_all       = lon_all[valid_data]
    #
    return sitenames_all,lat_all,lon_all

# Remove multiple proxy entrees at the same location
sitenames_all_temp12k,lat_all_temp12k,lon_all_temp12k = proxy_locations(sitenames_all_temp12k,lat_all_temp12k,lon_all_temp12k)
sitenames_all_shakun, lat_all_shakun, lon_all_shakun  = proxy_locations(sitenames_all_shakun, lat_all_shakun, lon_all_shakun)
sitenames_all_marcott,lat_all_marcott,lon_all_marcott = proxy_locations(sitenames_all_marcott,lat_all_marcott,lon_all_marcott)


### CALCULATIONS

# Remove a reference period from each 2k method
indices_ref_2k = np.where((years_2k >= reference_period[0]) & (years_2k <= reference_period[1]))[0]
for i in range(data_gmst_all.shape[0]):
    value_ref = np.mean(np.median(data_gmst_all[i,indices_ref_2k,:],axis=1),axis=0)
    data_gmst_all[i,:,:] = data_gmst_all[i,:,:]-value_ref

# Put all methods on the same axis
data_gmst_all_swap = np.swapaxes(data_gmst_all,0,1)
data_gmst_all_2d = np.reshape(data_gmst_all_swap,(ntime,nfiles*niter))

# This function takes a time-lat-lon variable and computes the global-mean.
def global_mean(variable,lats):
    variable_zonal = np.nanmean(variable,axis=2)
    lat_weights = np.cos(np.radians(lats))
    variable_global = np.average(variable_zonal,axis=1,weights=lat_weights)
    return variable_global

# Compute the annual-mean of the ERA-20C data
tas_era20c_globalmean = global_mean(tas_era20c,lat_era20c)

# Find the difference in means between two time series during their period of overlap
def mean_of_overlap(ts1,ages1,ts2,ages2):
    overlap_age_min = np.max([np.min(ages1),np.min(ages2)])
    overlap_age_max = np.min([np.max(ages1),np.max(ages2)])
    ts1_mean_overlap_period = np.mean(ts1[np.where((ages1 >= overlap_age_min) & (ages1 <= overlap_age_max))[0]])
    ts2_mean_overlap_period = np.mean(ts2[np.where((ages2 >= overlap_age_min) & (ages2 <= overlap_age_max))[0]])
    difference_in_means = ts1_mean_overlap_period - ts2_mean_overlap_period
    return difference_in_means

# Scale the mean of the ERA20C value to the 2k composite using their overlapping segments
difference_in_means_era20c = mean_of_overlap(tas_era20c_globalmean,ages_era20c,np.median(data_gmst_all_2d,axis=1),ages_2k)
tas_era20c_globalmean = tas_era20c_globalmean - difference_in_means_era20c

# Combine all 12k methods into one array
gmst_all_ensemble   = np.concatenate((gmst_scc_ensemble,gmst_dcc_ensemble,gmst_cps_ensemble,gmst_pai_ensemble,gmst_gam_ensemble),axis=1)
gmst_all_except_cps = np.concatenate((gmst_scc_ensemble,gmst_dcc_ensemble,gmst_pai_ensemble,gmst_gam_ensemble),axis=1)

# Test for NaNs
nantest_ensemble = copy.deepcopy(gmst_all_ensemble)
nantest_ensemble[~np.isnan(nantest_ensemble)] = 0
nantest_ensemble[np.isnan(nantest_ensemble)]  = 1
nan_counts = np.sum(nantest_ensemble,axis=0)
print('Extra NaNs at indices:',str(np.where(nan_counts > 2)[0]))  # GAM has NaN for 2 values.  This looks for other NaNs.
#plt.pcolormesh(nantest_ensemble)


# Compute the values at 6ka vs 1800-1900
indices_6k        = np.where((gmst_12k_ages >= 5500) & (gmst_12k_ages <= 6500))[0]  # 11 values
indices_1800_1900 = np.where((gmst_12k_ages >= 50)   & (gmst_12k_ages <= 150))[0]   # 1 value
values_6ka        = np.mean(gmst_all_ensemble[indices_6k,:],  axis=0) - np.mean(gmst_all_ensemble[indices_1800_1900,:],  axis=0)
values_6ka_scc    = np.mean(gmst_scc_ensemble[indices_6k,:],  axis=0) - np.mean(gmst_scc_ensemble[indices_1800_1900,:],  axis=0)
values_6ka_dcc    = np.mean(gmst_dcc_ensemble[indices_6k,:],  axis=0) - np.mean(gmst_dcc_ensemble[indices_1800_1900,:],  axis=0)
values_6ka_cps    = np.mean(gmst_cps_ensemble[indices_6k,:],  axis=0) - np.mean(gmst_cps_ensemble[indices_1800_1900,:],  axis=0)
values_6ka_pai    = np.mean(gmst_pai_ensemble[indices_6k,:],  axis=0) - np.mean(gmst_pai_ensemble[indices_1800_1900,:],  axis=0)
values_6ka_gam    = np.mean(gmst_gam_ensemble[indices_6k,:],  axis=0) - np.mean(gmst_gam_ensemble[indices_1800_1900,:],  axis=0)
values_6ka_no_cps = np.mean(gmst_all_except_cps[indices_6k,:],axis=0) - np.mean(gmst_all_except_cps[indices_1800_1900,:],axis=0)

# List the median and 90% ranges
print('6k - present, median and 90% ranges')
print('SCC,         median: '+str('{:.2f}'.format(np.median(values_6ka_scc)))+', range: '+str('{:.2f}'.format(np.percentile(values_6ka_scc,5)))+' - '+str('{:.2f}'.format(np.percentile(values_6ka_scc,95))))
print('DCC,         median: '+str('{:.2f}'.format(np.median(values_6ka_dcc)))+', range: '+str('{:.2f}'.format(np.percentile(values_6ka_dcc,5)))+' - '+str('{:.2f}'.format(np.percentile(values_6ka_dcc,95))))
print('GAM,         median: '+str('{:.2f}'.format(np.median(values_6ka_gam)))+', range: '+str('{:.2f}'.format(np.percentile(values_6ka_gam,5)))+' - '+str('{:.2f}'.format(np.percentile(values_6ka_gam,95))))
print('CPS,         median: '+str('{:.2f}'.format(np.median(values_6ka_cps)))+', range: '+str('{:.2f}'.format(np.percentile(values_6ka_cps,5)))+' - '+str('{:.2f}'.format(np.percentile(values_6ka_cps,95))))
print('PAI,         median: '+str('{:.2f}'.format(np.median(values_6ka_pai)))+', range: '+str('{:.2f}'.format(np.percentile(values_6ka_pai,5)))+' - '+str('{:.2f}'.format(np.percentile(values_6ka_pai,95))))
print('---')
print('All methods, median: '+str('{:.2f}'.format(np.median(values_6ka)))       +', range: '+str('{:.2f}'.format(np.percentile(values_6ka,5)))       +' - '+str('{:.2f}'.format(np.percentile(values_6ka,95))))
print('No CPS,      median: '+str('{:.2f}'.format(np.median(values_6ka_no_cps)))+', range: '+str('{:.2f}'.format(np.percentile(values_6ka_no_cps,5)))+' - '+str('{:.2f}'.format(np.percentile(values_6ka_no_cps,95))))


# Compute the timing and magnitude of the warmest 200 and 1000 year segments
n_ens = gmst_all_ensemble.shape[1]
max_1000yr_values = np.zeros((n_ens)); max_1000yr_values[:] = np.nan
max_1000yr_ages   = np.zeros((n_ens)); max_1000yr_ages[:]   = np.nan
max_200yr_values  = np.zeros((n_ens)); max_200yr_values[:]  = np.nan
max_200yr_ages    = np.zeros((n_ens)); max_200yr_ages[:]    = np.nan
gmst_12k_ages_1000yr_means = np.convolve(gmst_12k_ages,np.ones((10))/10,mode='valid')
gmst_12k_ages_200yr_means  = np.convolve(gmst_12k_ages,np.ones((2))/2,  mode='valid')
for i in range(n_ens):
    gmst_1000yr_means = np.convolve(gmst_all_ensemble[:,i],np.ones((10))/10,mode='valid') - np.mean(gmst_all_ensemble[indices_1800_1900,i],axis=0)
    gmst_200yr_means  = np.convolve(gmst_all_ensemble[:,i],np.ones((2))/2,  mode='valid') - np.mean(gmst_all_ensemble[indices_1800_1900,i],axis=0)
    max_1000yr_values[i] = np.nanmax(gmst_1000yr_means)
    max_200yr_values[i]  = np.nanmax(gmst_200yr_means)
    max_1000yr_ages[i] = gmst_12k_ages_1000yr_means[np.nanargmax(gmst_1000yr_means)]
    max_200yr_ages[i]  = gmst_12k_ages_200yr_means[np.nanargmax(gmst_200yr_means)]

# Determine how many of the 12k ensemble members ever exceed the 2001-2010 decade in ERA-20C
#recent_decade_indices = np.where((years_era20c >= 2001) & (years_era20c <= 2010))[0]
#tas_recent_decade = np.mean(tas_era20c_globalmean[recent_decade_indices]); label_txt = 'ERA-20C'
tas_recent_decade = 1. + 0.03; label_txt = 'WHO'  # This is a value from WHO (1.) plus a 0.03 differential to account for the different reference periods

gmst_all_ensemble_anom_to_1800_1900 = gmst_all_ensemble - np.mean(gmst_all_ensemble[indices_1800_1900,:],axis=0)[None,:]
gmst_all_ensemble_maxes = np.nanmax(gmst_all_ensemble_anom_to_1800_1900,axis=0)
ens_warmer_than_recent_100yr = sum(gmst_all_ensemble_maxes > tas_recent_decade)
ens_warmer_than_recent_200yr = sum(max_200yr_values        > tas_recent_decade)
percent_ens_warmer_than_recent_100yr = 100*(ens_warmer_than_recent_100yr/n_ens)
percent_ens_warmer_than_recent_200yr = 100*(ens_warmer_than_recent_200yr/n_ens)

print('===')
print('Percent of ensemble members with segments warmer than the 2001-2010 decade (recent value based on '+label_txt+')')
print('  For 100yr segments:',percent_ens_warmer_than_recent_100yr)
print('  For 200yr segments:',percent_ens_warmer_than_recent_200yr)
print('===')

# Determine which methods contain the warmer values
index_of_warmer = np.where(max_200yr_values > tas_recent_decade)[0]
print('Number of ensemble members with 200yr segments warmer than the 2001-2010 decade (recent value based on '+label_txt+')')
print('SCC: ',index_of_warmer[(index_of_warmer >= 0)    & (index_of_warmer < 500)].shape[0])
print('DCC: ',index_of_warmer[(index_of_warmer >= 500)  & (index_of_warmer < 1000)].shape[0])
print('CPS: ',index_of_warmer[(index_of_warmer >= 1000) & (index_of_warmer < 1500)].shape[0])
print('PAI: ',index_of_warmer[(index_of_warmer >= 1500) & (index_of_warmer < 2000)].shape[0])
print('GAM: ',index_of_warmer[(index_of_warmer >= 2000) & (index_of_warmer < 2500)].shape[0])

# Relative warmth of recent decade
print('Relative warmth of recent decade in',label_txt,':',tas_recent_decade)

# Calculate the change between ~12ka and 10ka for each method
index_10000 = np.where(gmst_12k_ages == 10000)[0][0]
index_12000 = np.where(gmst_12k_ages == 12000)[0][0]
index_11800 = np.where(gmst_12k_ages == 11800)[0][0]
print('Temperature change between ~12ka and 10ka for each method')
print('SSC: ',str(np.nanmedian(gmst_scc_ensemble[index_10000,:],axis=0) - np.nanmedian(gmst_scc_ensemble[index_12000,:],axis=0)))
print('DCC: ',str(np.nanmedian(gmst_dcc_ensemble[index_10000,:],axis=0) - np.nanmedian(gmst_dcc_ensemble[index_12000,:],axis=0)))
print('CPS: ',str(np.nanmedian(gmst_cps_ensemble[index_10000,:],axis=0) - np.nanmedian(gmst_cps_ensemble[index_12000,:],axis=0)))
print('PAI: ',str(np.nanmedian(gmst_pai_ensemble[index_10000,:],axis=0) - np.nanmedian(gmst_pai_ensemble[index_12000,:],axis=0)))
print('GAM: ',str(np.nanmedian(gmst_gam_ensemble[index_10000,:],axis=0) - np.nanmedian(gmst_gam_ensemble[index_11800,:],axis=0)))  # GAM doesn't have values all the way back to 12ka

# Average the Marcott reconstruction into 120 year bins, centered on the same ages at the Temp12k reconstruction
# Note: 120 year means are used rather than 100 year means because of the original ages of the Marcott reconstruction.
ages_marcott_central = np.arange(0,11301,100)
globalmean_marcott_120yrmeans = np.zeros((len(ages_marcott_central))); globalmean_marcott_120yrmeans[:] = np.nan
onesigma_marcott_120yrmeans   = np.zeros((len(ages_marcott_central))); onesigma_marcott_120yrmeans[:]   = np.nan
for i,central_age in enumerate(ages_marcott_central):
    age_indices_selected = np.where((ages_marcott >= (central_age-50)) & (ages_marcott <= (central_age+50)))[0]
    if len(age_indices_selected) != 6: print('Marcott average centered on age '+str(central_age)+' is using '+str(len(age_indices_selected))+' points')
    globalmean_marcott_120yrmeans[i] = np.mean(globalmean_marcott[age_indices_selected])
    onesigma_marcott_120yrmeans[i]   = np.mean(onesigma_marcott[age_indices_selected])

ages_marcott_120yrmeans = ages_marcott_central

# Remove a reference period from the Marcott reconstruction
reference_ages_12k_medians = [50,150]
indices_ref_marcott = np.where((ages_marcott_120yrmeans >= reference_ages_12k_medians[0]) & (ages_marcott_120yrmeans <= reference_ages_12k_medians[1]))[0]
globalmean_marcott_20yrmean   = globalmean_marcott            - np.mean(globalmean_marcott_120yrmeans[indices_ref_marcott])
globalmean_marcott_120yrmeans = globalmean_marcott_120yrmeans - np.mean(globalmean_marcott_120yrmeans[indices_ref_marcott])

# The Shukun reconstruction stops at 6.5 ka.  Align it to the Marcott reconstruction using the period of overlap
difference_in_means_shakun = mean_of_overlap(globalmean_shakun,ages_shakun,globalmean_marcott_120yrmeans,ages_marcott_120yrmeans)
globalmean_shakun = globalmean_shakun - difference_in_means_shakun

# In Marcott's reconstruction, find the largest 200-year mean value
marcott_200yr_means  = np.convolve(globalmean_marcott_20yrmean,np.ones((10))/10,mode='valid')
temp12k_200yr_means  = np.convolve(np.nanmedian(gmst_all_ensemble,axis=1),np.ones((2))/2,mode='valid')
print('Maximum 200-year mean in the Marcott reconstruction: ',str(max(marcott_200yr_means)))
print('Maximum 200-year mean in the Temp12k reconstruction: ',str(max(temp12k_200yr_means)))
print('Difference: ',str(max(marcott_200yr_means)-max(temp12k_200yr_means)))

# Calculate the 5.5-6.5ka mean from Marcott.
indices_6ka_marcott = np.where((ages_marcott >= 5500) & (ages_marcott <= 6500))[0]
print('Mean of Marcott for 5.5-6.5 ka:',str(np.mean(globalmean_marcott_20yrmean[indices_6ka_marcott])))


# Make sure that Holocene proxy values go from -180 to 180
print(np.min(lon_all_temp12k),np.max(lon_all_temp12k))
print(np.min(lon_all_shakun),np.max(lon_all_shakun))
print(np.min(lon_all_marcott),np.max(lon_all_marcott))

# Count the proxies
nproxies_temp12k = len(lat_all_temp12k)
nproxies_shakun  = len(lat_all_shakun)
nproxies_marcott = len(lat_all_marcott)
print('Number of proxies:',nproxies_temp12k,nproxies_shakun,nproxies_marcott)



### FIGURES
plt.style.use('ggplot')

# Plot time series of the reconstructions
f, ax1 = plt.subplots(1,1,figsize=(16,14),sharex=True)

line_era20c, = ax1.plot(ages_era20c,tas_era20c_globalmean,             color='k',   linewidth=1)
line_2k,     = ax1.plot(ages_2k,    np.median(data_gmst_all_2d,axis=1),color='navy',linewidth=.5)
for i in range(1,10):
    percentiles_12k = ax1.fill_between(gmst_12k_ages,np.nanpercentile(gmst_all_ensemble,i*5,axis=1),np.nanpercentile(gmst_all_ensemble,100-(i*5),axis=1),color='k',alpha=0.12)

outer_percentiles_12k, = ax1.plot(gmst_12k_ages,np.nanpercentile(gmst_all_ensemble,5,axis=1), '--',color='gray',linewidth=1)
outer_percentiles_12k, = ax1.plot(gmst_12k_ages,np.nanpercentile(gmst_all_ensemble,95,axis=1),'--',color='gray',linewidth=1)
line_scc_12k, = ax1.plot(gmst_12k_ages,np.nanmedian(gmst_scc_ensemble,axis=1),color='tab:blue',  linewidth=3)
line_dcc_12k, = ax1.plot(gmst_12k_ages,np.nanmedian(gmst_dcc_ensemble,axis=1),color='tab:cyan',  linewidth=3)
line_gam_12k, = ax1.plot(gmst_12k_ages,np.nanmedian(gmst_gam_ensemble,axis=1),color='tab:olive', linewidth=3)
line_cps_12k, = ax1.plot(gmst_12k_ages,np.nanmedian(gmst_cps_ensemble,axis=1),color='tab:orange',linewidth=3)
line_pai_12k, = ax1.plot(gmst_12k_ages,np.nanmedian(gmst_pai_ensemble,axis=1),color='tab:red',   linewidth=3)
ax1.axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)

ax1.legend((line_scc_12k,line_dcc_12k,line_gam_12k,line_cps_12k,line_pai_12k,line_era20c,line_2k,outer_percentiles_12k,percentiles_12k),('SCC','DCC','GAM','CPS','PAI','ERA-20C','2k median','All 12k methods, 5 - 95th percentiles','All 12k methods, every 5th percentile'),loc=7,ncol=2,fontsize=16)
ax1.set_xlabel('Age (yr BP)',fontsize=20)
ax1.set_ylabel('$\Delta$ Global-mean temperature ($^\circ$C)',fontsize=20)
ax1.set_xlim(12000,-60)
ax1.set_ylim(-5,2)
ax1.set_title('Global-mean temperature composites',fontsize=26)
ax1.tick_params(axis='both',which='major',labelsize=20)

# Create inset plot
ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(ax1,[.25,.05,.73,.36])
ax2.set_axes_locator(ip)

line_era20c, = ax2.plot(ages_era20c,tas_era20c_globalmean,             color='k',   linewidth=1)
line_2k,     = ax2.plot(ages_2k,    np.median(data_gmst_all_2d,axis=1),color='navy',linewidth=.5)
for i in range(1,10):
    percentiles_12k = ax2.fill_between(gmst_12k_ages,np.nanpercentile(gmst_all_ensemble,i*5,axis=1),np.nanpercentile(gmst_all_ensemble,100-(i*5),axis=1),color='k',alpha=0.12)

outer_percentiles_12k, = ax2.plot(gmst_12k_ages,np.nanpercentile(gmst_all_ensemble,5,axis=1), '--',color='gray',linewidth=1)
outer_percentiles_12k, = ax2.plot(gmst_12k_ages,np.nanpercentile(gmst_all_ensemble,95,axis=1),'--',color='gray',linewidth=1)
line_scc_12k, = ax2.plot(gmst_12k_ages,np.nanmedian(gmst_scc_ensemble,axis=1),color='tab:blue',  linewidth=3)
line_dcc_12k, = ax2.plot(gmst_12k_ages,np.nanmedian(gmst_dcc_ensemble,axis=1),color='tab:cyan',  linewidth=3)
line_gam_12k, = ax2.plot(gmst_12k_ages,np.nanmedian(gmst_gam_ensemble,axis=1),color='tab:olive', linewidth=3)
line_cps_12k, = ax2.plot(gmst_12k_ages,np.nanmedian(gmst_cps_ensemble,axis=1),color='tab:orange',linewidth=3)
line_pai_12k, = ax2.plot(gmst_12k_ages,np.nanmedian(gmst_pai_ensemble,axis=1),color='tab:red',   linewidth=3)
ax2.axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)

ax2.set_xlabel('Age (yr BP)',fontsize=12)
ax2.set_ylabel('$\Delta$ Global-mean temperature ($^\circ$C)',fontsize=12)
ax2.set_xlim(2000,-60)
ax2.set_ylim(-.4,1.3)
ax2.set_title('  Close-up of values since 2 ka',fontsize=12,pad=-20,loc='left',color='gray')
ax2.tick_params(axis='both',which='major',labelsize=12)

if save_instead_of_plot:
    plt.savefig('figures/Fig3_gmst_composite_12k.png',dpi=300,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()


print(np.min([np.min(values_6ka),np.min(max_1000yr_values),np.min(max_200yr_values)]))
print(np.max([np.max(values_6ka),np.max(max_1000yr_values),np.max(max_200yr_values)]))

# Plot histograms of values
f, ax = plt.subplots(2,1,figsize=(12,8),sharex=False,sharey=False)
ax = ax.ravel()
bins_values = np.arange(-.5,3.51,.1)
bins_ages   = np.arange(0,12001,100)

ax[0].hist(values_6ka,       bins=bins_values,color='k',histtype='step',linewidth=1.5,label='6ka;                    median: '+str('{:.2f}'.format(np.median(values_6ka)))+'$^\circ$C,   range: '+str('{:.2f}'.format(np.percentile(values_6ka,5)))+' - '+str('{:.2f}'.format(np.percentile(values_6ka,95)))+'$^\circ$C',zorder=2)
ax[0].hist(max_1000yr_values,bins=bins_values,color='tab:red', alpha=0.5,label='1000yr means;   median: '+str('{:.2f}'.format(np.median(max_1000yr_values)))+'$^\circ$C,   range: '+str('{:.2f}'.format(np.percentile(max_1000yr_values,5)))+' - '+str('{:.2f}'.format(np.percentile(max_1000yr_values,95)))+'$^\circ$C')
ax[0].hist(max_200yr_values, bins=bins_values,color='tab:blue',alpha=0.5,label='200yr means;     median: '+str('{:.2f}'.format(np.median(max_200yr_values)))+'$^\circ$C,   range: '+str('{:.2f}'.format(np.percentile(max_200yr_values,5)))+' - '+str('{:.2f}'.format(np.percentile(max_200yr_values,95)))+'$^\circ$C')
ax[0].axvline(x=np.median(values_6ka),       color='k',        alpha=1,linestyle=':',linewidth=2)
ax[0].axvline(x=np.median(max_1000yr_values),color='tab:red',  alpha=1,linestyle=':',linewidth=2)
ax[0].axvline(x=np.median(max_200yr_values), color='tab:blue', alpha=1,linestyle=':',linewidth=2)
ax[0].set_xlim(-.5,3.5)
ax[0].set_title('a) Warmest values',loc='left',fontsize=16)
ax[0].set_xlabel('$\Delta$T ($^\circ$C)',fontsize=16)
ax[0].axvline(x=0,color='k',linestyle='--')
ax[0].legend(loc=1)

ax[1].hist(max_1000yr_ages,bins=bins_ages,color='tab:red',  alpha=0.5,label='Ages of 1000yr means;   median: '+str('{:.0f}'.format(np.median(max_1000yr_ages)))+' yr BP,   range: '+str('{:.0f}'.format(np.percentile(max_1000yr_ages,5)))+' - '+str('{:.0f}'.format(np.percentile(max_1000yr_ages,95)))+' yr BP')
ax[1].hist(max_200yr_ages, bins=bins_ages,color='tab:blue', alpha=0.5,label='Ages of 200yr means;     median: '+str('{:.0f}'.format(np.median(max_200yr_ages)))+' yr BP,   range: '+str('{:.0f}'.format(np.percentile(max_200yr_ages,5)))+' - '+str('{:.0f}'.format(np.percentile(max_200yr_ages,95)))+' yr BP')
ax[1].axvline(x=np.median(max_1000yr_ages),color='tab:red', alpha=1,linestyle=':',linewidth=2)
ax[1].axvline(x=np.median(max_200yr_ages), color='tab:blue',alpha=1,linestyle=':',linewidth=2)
ax[1].set_xlim(12000,0)
ax[1].set_title('b) Ages of warmest values',loc='left',fontsize=16)
ax[1].set_xlabel('Age (yr BP)',fontsize=16)
ax[1].legend(loc=1)

for i in range(2):
    ax[i].set_ylabel('Frequency',fontsize=16)
    ax[i].tick_params(axis='both',which='major',labelsize=14)

plt.suptitle('Values and ages of maximum warmth',fontsize=20)
f.tight_layout()
f.subplots_adjust(top=.9)

if save_instead_of_plot:
    plt.savefig('figures/Fig4_max_warmth_histograms.png',dpi=200,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()


# Compare the new Temp12k reconstructions against older Holocene reconstructions
f, ax = plt.subplots(2,1,figsize=(16,18),sharex=False,sharey=False)
ax = ax.ravel()

range_temp12k = ax[0].fill_between(gmst_12k_ages,          np.nanpercentile(gmst_all_ensemble,16,axis=1),            np.nanpercentile(gmst_all_ensemble,84,axis=1),            color='darkgreen',  alpha=0.5)
range_shakun  = ax[0].fill_between(ages_shakun,            globalmean_shakun-onesigma_shakun,                        globalmean_shakun+onesigma_shakun,                        color='tab:purple', alpha=0.5)
range_marcott = ax[0].fill_between(ages_marcott_120yrmeans,globalmean_marcott_120yrmeans-onesigma_marcott_120yrmeans,globalmean_marcott_120yrmeans+onesigma_marcott_120yrmeans,color='deepskyblue',alpha=0.5)

line_temp12k, = ax[0].plot(gmst_12k_ages,          np.nanmedian(gmst_all_ensemble,axis=1),color='darkgreen',  linewidth=3)
line_shakun,  = ax[0].plot(ages_shakun,            globalmean_shakun,                     color='tab:purple', linewidth=3)
line_marcott, = ax[0].plot(ages_marcott_120yrmeans,globalmean_marcott_120yrmeans,         color='deepskyblue',linewidth=3)
ax[0].axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)

ax[0].legend((line_temp12k,line_shakun,line_marcott),('Temp12k reconstruction (this study)','Shakun et al. 2012','Marcott et al. 2013 (120yr means)'),loc='lower right',ncol=1,fontsize=16)
ax[0].set_xlabel('Age (yr BP)',fontsize=16)
ax[0].set_ylabel('$\Delta$ Global-mean temperature ($^\circ$C)',fontsize=16)
ax[0].set_xlim(12000,0)
ax[0].set_ylim(-1.1,1.1)
ax[0].set_title('(a) Comparison of Holocene global-mean temperature composites',fontsize=20)
ax[0].tick_params(axis='both',which='major',labelsize=16)

# Plot the proxy locations.
m = Basemap(projection='robin',lon_0=0,resolution='c',ax=ax[1])
x_temp12k, y_temp12k = m(lon_all_temp12k, lat_all_temp12k)
x_shakun,  y_shakun  = m(lon_all_shakun,  lat_all_shakun)
x_marcott, y_marcott = m(lon_all_marcott, lat_all_marcott)

m.drawcoastlines()
m.fillcontinents(color='whitesmoke',lake_color='white',zorder=0)
scatter_shakun  = m.scatter(x_shakun, y_shakun, 300,color='tab:purple')
scatter_marcott = m.scatter(x_marcott,y_marcott,200,color='deepskyblue')
scatter_temp12k = m.scatter(x_temp12k,y_temp12k, 50,color='darkgreen')
legend = plt.legend([scatter_temp12k,scatter_shakun,scatter_marcott],['Temp12k ($N_{sites}$='+str(nproxies_temp12k)+')','Shakun ($N_{sites}$='+str(nproxies_shakun)+')','Marcott ($N_{sites}$='+str(nproxies_marcott)+')'],scatterpoints=1,loc=1,ncol=1,bbox_to_anchor=(.23,.25),prop={'size':14})

m.drawparallels(np.arange(-90,91,30),  labels=[True, False,False,False],fontsize=12)
m.drawmeridians(np.arange(-180,181,30),labels=[False,False,False,True ],fontsize=12)
ax[1].set_title('b) Proxy locations',loc='center',fontsize=20)

if save_instead_of_plot:
    plt.savefig('figures/Fig6_gmst_composite_12k_comparison_and_map.png',dpi=300,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()

