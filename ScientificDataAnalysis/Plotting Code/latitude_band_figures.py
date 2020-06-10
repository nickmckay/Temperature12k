#=============================================================================
# This script makes some figures of the Temp12k latitude-band data.
#    author: Michael P. Erb
#    date  : 6/8/2020
#=============================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import xarray as xr

save_instead_of_plot = True
reference_period = [1800,1899]


### LOAD DATA
data_dir = '/home/mpe32/analysis/14_Holocene_proxies/GMST_paper/data/'

# Load 12k data
handle = xr.open_dataset(data_dir+'final_data/temp12k_alldata.nc',decode_times=False)
temp_12k_ages = handle['age'].values
scc_latbands = handle['scc_latbands'].values
dcc_latbands = handle['dcc_latbands'].values
cps_latbands = handle['cps_latbands'].values
pai_latbands = handle['pai_latbands'].values
gam_latbands = handle['gam_latbands'].values
latbands_all    = handle['latband_ranges'].values
latband_weights = handle['latband_weights'].values
gmst_scc_ensemble = handle['scc_globalmean'].values
gmst_dcc_ensemble = handle['dcc_globalmean'].values
gmst_cps_ensemble = handle['cps_globalmean'].values
gmst_pai_ensemble = handle['pai_globalmean'].values
gmst_gam_ensemble = handle['gam_globalmean'].values
handle.close()

# Load latitude band data for 2k reconstructions and ERA-20C
data_dir_lat_bands = '/home/mpe32/analysis/14_Holocene_proxies/GMST_paper/latitude_bands/data_lat_bands/'
temp_2k_latbands     = np.load(data_dir_lat_bands+'tas_lat_bands_2k.npy',allow_pickle=True).item()
temp_era20c_latbands = np.load(data_dir_lat_bands+'tas_lat_bands_era20c.npy',allow_pickle=True).item()



### CALCULATIONS

# Reorganize the 12k data
temp_scc_latbands = {}
temp_dcc_latbands = {}
temp_cps_latbands = {}
temp_pai_latbands = {}
temp_gam_latbands = {}
for i,latband_txt in enumerate(latbands_all):
    temp_scc_latbands[latband_txt] = scc_latbands[i,:,:]
    temp_dcc_latbands[latband_txt] = dcc_latbands[i,:,:]
    temp_cps_latbands[latband_txt] = cps_latbands[i,:,:]
    temp_pai_latbands[latband_txt] = pai_latbands[i,:,:]
    temp_gam_latbands[latband_txt] = gam_latbands[i,:,:]

# Compute years and ages
temp_12k_years               = 1950-temp_12k_ages
temp_2k_latbands['ages']     = 1950-temp_2k_latbands['years']
temp_era20c_latbands['ages'] = 1950-temp_era20c_latbands['years']

# Remove a reference period every ensemble member of the 2k data
indices_ref_2k = np.where((temp_2k_latbands['years'] >= reference_period[0]) & (temp_2k_latbands['years'] <= reference_period[1]))[0]
for latband_txt in latbands_all:
    temp_2k_latbands[latband_txt] = temp_2k_latbands[latband_txt] - np.mean(temp_2k_latbands[latband_txt][indices_ref_2k,:],axis=0)[None,:]

# Scale the mean of the ERA20C value to the 2k composite using their overlapping segments
def mean_of_overlap(ts1,ages1,ts2,ages2):
    overlap_age_min = np.max([np.min(ages1),np.min(ages2)])
    overlap_age_max = np.min([np.max(ages1),np.max(ages2)])
    ts1_mean_overlap_period = np.mean(ts1[np.where((ages1 >= overlap_age_min) & (ages1 <= overlap_age_max))[0]])
    ts2_mean_overlap_period = np.mean(ts2[np.where((ages2 >= overlap_age_min) & (ages2 <= overlap_age_max))[0]])
    difference_in_means = ts1_mean_overlap_period - ts2_mean_overlap_period
    return difference_in_means

for latband_txt in latbands_all:
    difference_in_means_era20c = mean_of_overlap(temp_era20c_latbands[latband_txt],temp_era20c_latbands['ages'],np.median(temp_2k_latbands[latband_txt],axis=1),temp_2k_latbands['ages'])
    temp_era20c_latbands[latband_txt] = temp_era20c_latbands[latband_txt] - difference_in_means_era20c

# Combine all 12k methods into one array
n_latbands = 6
n_methods  = 5
n_ens      = 500
n_ages     = len(temp_12k_ages)
temp_all_latbands = np.zeros((n_latbands,n_ages,n_methods*n_ens)); temp_all_latbands[:] = np.nan
for i,latband_txt in enumerate(latbands_all):
    print(i,latband_txt,latband_weights[i])
    temp_all_latbands[i,:,:] = np.concatenate((temp_scc_latbands[latband_txt],temp_dcc_latbands[latband_txt],temp_cps_latbands[latband_txt],temp_pai_latbands[latband_txt],temp_gam_latbands[latband_txt]),axis=1)

# Compute spatial means for the NH, SH, and the globe
temp_all_NH     = np.average(temp_all_latbands[0:3,:,:],axis=0,weights=latband_weights[0:3])
temp_all_SH     = np.average(temp_all_latbands[3:6,:,:],axis=0,weights=latband_weights[3:6])
temp_all_global = np.concatenate((gmst_scc_ensemble,gmst_dcc_ensemble,gmst_cps_ensemble,gmst_pai_ensemble,gmst_gam_ensemble),axis=1)

# Compute trends from 6ka to 0ka (excluding the warming signal of the recent past) for each case
indices_in_interval = np.where((temp_12k_ages >= 100) & (temp_12k_ages <= 6000))[0]
trends_global = np.zeros((n_methods*n_ens)); trends_global[:] = np.nan
trends_NH     = np.zeros((n_methods*n_ens)); trends_NH[:]     = np.nan
trends_SH     = np.zeros((n_methods*n_ens)); trends_SH[:]     = np.nan
for i in range(n_methods*n_ens):
    trends_global[i],_,_,_,_ = stats.linregress(-1*temp_12k_ages[indices_in_interval],temp_all_global[indices_in_interval,i])
    trends_NH[i],_,_,_,_     = stats.linregress(-1*temp_12k_ages[indices_in_interval],temp_all_NH[indices_in_interval,i])
    trends_SH[i],_,_,_,_     = stats.linregress(-1*temp_12k_ages[indices_in_interval],temp_all_SH[indices_in_interval,i])

# Calculate the percentile of the 0 value
percentile_0_global = stats.percentileofscore(trends_global,0,kind='strict')
percentile_0_NH     = stats.percentileofscore(trends_NH,    0,kind='strict')
percentile_0_SH     = stats.percentileofscore(trends_SH,    0,kind='strict')



### FIGURES
plt.style.use('ggplot')

latband_titles = ['60-90$^\circ$N','30-60$^\circ$N','0-30$^\circ$N','0-30$^\circ$S','30-60$^\circ$S','60-90$^\circ$S']

# Plot the latitude band value for NH and SH
f, ax = plt.subplots(3,2,figsize=(20,15),sharex=True,sharey=True)
ax = ax.ravel()

num_range = [0,5,1,4,2,3]
letters = ['a','d','b','e','c','f']
for i,latband_num in enumerate(num_range):
    #
    latband_txt   = latbands_all[latband_num]
    latband_title = latband_titles[latband_num]
    #
    line_era20c, = ax[i].plot(temp_era20c_latbands['ages'],temp_era20c_latbands[latband_txt],              color='k',   linewidth=1)
    line_2k,     = ax[i].plot(temp_2k_latbands['ages'],    np.median(temp_2k_latbands[latband_txt],axis=1),color='navy',linewidth=.5)
    for j in range(1,10):
        percentiles_12k = ax[i].fill_between(temp_12k_ages,np.nanpercentile(temp_all_latbands[latband_num,:,:],j*5,axis=1),np.nanpercentile(temp_all_latbands[latband_num,:,:],100-(j*5),axis=1),color='k',alpha=0.12)
    #
    outer_percentiles_12k, = ax[i].plot(temp_12k_ages,np.nanpercentile(temp_all_latbands[latband_num,:,:],5,axis=1), '--',color='gray',linewidth=1)
    outer_percentiles_12k, = ax[i].plot(temp_12k_ages,np.nanpercentile(temp_all_latbands[latband_num,:,:],95,axis=1),'--',color='gray',linewidth=1)
    line_scc_12k, = ax[i].plot(temp_12k_ages,np.nanmedian(temp_scc_latbands[latband_txt],axis=1),color='tab:blue',  linewidth=3)
    line_dcc_12k, = ax[i].plot(temp_12k_ages,np.nanmedian(temp_dcc_latbands[latband_txt],axis=1),color='tab:cyan',  linewidth=3)
    line_gam_12k, = ax[i].plot(temp_12k_ages,np.nanmedian(temp_gam_latbands[latband_txt],axis=1),color='tab:olive', linewidth=3)
    line_cps_12k, = ax[i].plot(temp_12k_ages,np.nanmedian(temp_cps_latbands[latband_txt],axis=1),color='tab:orange',linewidth=3)
    line_pai_12k, = ax[i].plot(temp_12k_ages,np.nanmedian(temp_pai_latbands[latband_txt],axis=1),color='tab:red',   linewidth=3)
    ax[i].axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)
    #
    ax[i].set_xlim(12000,-60)
    ax[i].set_ylim(-6,6)
    ax[i].set_title(letters[i]+') '+latband_title,fontsize=20,loc='left')
    ax[i].tick_params(axis='both',which='major',labelsize=20)

ax[0].legend((line_scc_12k,line_dcc_12k,line_gam_12k,line_cps_12k,line_pai_12k,line_era20c,line_2k,outer_percentiles_12k,percentiles_12k),('SCC','DCC','GAM','CPS','PAI','ERA-20C','2k median','All 12k methods, 5 - 95th percentiles','All 12k methods, every 5th percentile'),ncol=2,loc=8,fontsize=13)
ax[2].set_ylabel('$\Delta$ Temperature ($^\circ$C)',fontsize=20)
ax[4].set_xlabel('Age (yr BP)',fontsize=20)
ax[5].set_xlabel('Age (yr BP)',fontsize=20)
plt.suptitle('Temperature composites for latitude bands',fontsize=28)
f.tight_layout()
f.subplots_adjust(top=.91)
if save_instead_of_plot:
    plt.savefig('figures/Fig2_latitude_composites_12k.png',dpi=300,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()


print('Min:',min([min(trends_global),min(trends_NH),min(trends_SH)])*1000)
print('Max:',max([max(trends_global),max(trends_NH),max(trends_SH)])*1000)


# Plot histograms of trends
f, ax = plt.subplots(2,1,figsize=(12,8),sharex=True,sharey=True)
ax = ax.ravel()
bins = np.arange(-1,.21,.025)

ax[0].hist(trends_global*1000,bins=bins,color='k',alpha=0.5,label='Global-mean;   median: '+str('{:.2f}'.format(np.median(trends_global*1000)))+',   range: '+str('{:.2f}'.format(np.percentile(trends_global*1000,5)))+' - '+str('{:.2f}'.format(np.percentile(trends_global*1000,95))),zorder=2)
ax[0].axvline(x=np.median(trends_global*1000),color='k',alpha=0.75,linestyle=':',linewidth=2)
ax[0].set_title('a) Global-mean',loc='left',fontsize=16)

ax[1].hist(trends_NH*1000,bins=bins,color='darkgreen',alpha=0.5,label='NH-mean;   median: '+str('{:.2f}'.format(np.median(trends_NH*1000)))+',   range: '+str('{:.2f}'.format(np.percentile(trends_NH*1000,5)))+' - '+str('{:.2f}'.format(np.percentile(trends_NH*1000,95))),zorder=2)
ax[1].hist(trends_SH*1000,bins=bins,color='purple',   alpha=0.5,label='SH-mean;   median: '+str('{:.2f}'.format(np.median(trends_SH*1000)))+',   range: '+str('{:.2f}'.format(np.percentile(trends_SH*1000,5)))+' - '+str('{:.2f}'.format(np.percentile(trends_SH*1000,95))),zorder=2)
ax[1].axvline(x=np.median(trends_NH*1000),color='darkgreen',alpha=0.75,linestyle=':',linewidth=2)
ax[1].axvline(x=np.median(trends_SH*1000),color='purple',   alpha=0.75,linestyle=':',linewidth=2)
ax[1].set_title('b) Hemispheric means',loc='left',fontsize=16)
ax[1].set_xlabel('Rate of $\Delta$T ($^\circ$C per millennia)',fontsize=16)

for i in range(2):
    ax[i].axvline(x=0,color='k',linestyle='--')
    ax[i].set_ylabel('Frequency',fontsize=16)
    ax[i].legend(loc=2,fontsize=12)
    ax[i].set_xlim(-1,.2)
    ax[i].tick_params(axis='both',which='major',labelsize=14)

plt.suptitle('Temperature trends from 6 to 0.1 ka',fontsize=20)
f.tight_layout()
f.subplots_adjust(top=.9)

if save_instead_of_plot:
    plt.savefig('figures/Fig5_trend_6_0ka_histograms_nmethods_'+str(n_methods)+'.png',dpi=200,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()


print('Percent of ensemble members which show cooling')
print('  Global-mean: '+str('{:.2f}'.format(percentile_0_global)+'%'))
print('  NH:          '+str('{:.2f}'.format(percentile_0_NH)+'%'))
print('  SH:          '+str('{:.2f}'.format(percentile_0_SH)+'%'))


# An extra figure, for comparison: plot a latitude by latitude comparison of the different methods
f, ax = plt.subplots(6,2,figsize=(15,25))
ax = ax.ravel()

for i,latband_txt in enumerate(latbands_all):
    #
    ax[i*2].fill_between(temp_2k_latbands['years'],np.percentile(temp_2k_latbands[latband_txt],5,axis=1),np.percentile(temp_2k_latbands[latband_txt],95,axis=1),color='navy',linewidth=0,alpha=0.2)
    line_2k,      = ax[i*2].plot(temp_2k_latbands['years'],np.nanmedian(temp_2k_latbands[latband_txt],axis=1), color='navy',      linewidth=1)
    line_scc_12k, = ax[i*2].plot(temp_12k_years,           np.nanmedian(temp_scc_latbands[latband_txt],axis=1),color='tab:blue',  linewidth=3)
    line_dcc_12k, = ax[i*2].plot(temp_12k_years,           np.nanmedian(temp_dcc_latbands[latband_txt],axis=1),color='tab:gray',  linewidth=3)
    line_gam_12k, = ax[i*2].plot(temp_12k_years,           np.nanmedian(temp_gam_latbands[latband_txt],axis=1),color='tab:olive', linewidth=3)
    line_cps_12k, = ax[i*2].plot(temp_12k_years,           np.nanmedian(temp_cps_latbands[latband_txt],axis=1),color='tab:orange',linewidth=3)
    line_pai_12k, = ax[i*2].plot(temp_12k_years,           np.nanmedian(temp_pai_latbands[latband_txt],axis=1),color='tab:red',   linewidth=3)
    ax[i*2].axhline(y=0,color='gray',linewidth=1,linestyle='dashed')
    if i == 0: ax[i*2].legend((line_2k,line_scc_12k,line_dcc_12k,line_gam_12k,line_cps_12k,line_pai_12k),('2k Median','SCC','DCC','GAM','CPS','PAI'),loc=3,fontsize=10,ncol=2)
    #
    ax[i*2].set_ylabel('$\Delta$T ($^\circ$C)',fontsize=12)
    if i == 5: ax[i*2].set_xlabel('Year (CE)',fontsize=12)
    ax[i*2].set_xlim(0,2000)
    ax[i*2].set_ylim(-2,2)
    ax[i*2].set_title(latband_titles[i],fontsize=12)
    #
    ax[(i*2)+1].fill_between(temp_2k_latbands['years'],np.percentile(temp_2k_latbands[latband_txt],5,axis=1),np.percentile(temp_2k_latbands[latband_txt],95,axis=1),color='navy',linewidth=0,alpha=0.2)
    line_2k,     = ax[(i*2)+1].plot(temp_2k_latbands['years'],    np.nanmedian(temp_2k_latbands[latband_txt],axis=1),color='navy',linewidth=1.5)
    line_era20c, = ax[(i*2)+1].plot(temp_era20c_latbands['years'],temp_era20c_latbands[latband_txt],                 color='k',   linewidth=1.5)
    ax[(i*2)+1].axhline(y=0,color='gray',linewidth=1,linestyle='dashed')
    if i == 0: ax[(i*2)+1].legend((line_era20c,line_2k),('ERA-20C','2k Median'),loc=4,fontsize=10,ncol=2)
    #
    ax[(i*2)+1].set_ylabel('$\Delta$T ($^\circ$C)',fontsize=12)
    if i == 5: ax[(i*2)+1].set_xlabel('Year (CE)',fontsize=12)
    ax[(i*2)+1].set_xlim(1900,2010)
    ax[(i*2)+1].set_ylim(-2,4)
    ax[(i*2)+1].set_title(latband_titles[i],fontsize=12)

plt.suptitle('Comparison of latitude-band time series during periods of overlap\nshowing Temp12k reconstructions, 2k reconstructions, and ERA-20C',fontsize=20)
f.tight_layout()
f.subplots_adjust(top=.94)

if save_instead_of_plot:
    plt.savefig('figures/ts_overlap_comparison.png',dpi=200,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()

