import numpy as np
import numpy.matlib
import math
from scipy import signal
import xarray as xr
from scipy.stats import gaussian_kde, ttest_ind
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.rasterize_and_save as raster
import cartopy
import cartopy.crs as ccrs
from datetime import timedelta
from dateutil.relativedelta import relativedelta, MO
import datetime
import matplotlib.animation as animation
from matplotlib import rc
from dask.distributed import Client, LocalCluster
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import os
import glob
import gc
import multiprocessing as mp
from itertools import repeat
import time
import warnings

start = time.time()
# Set model, variables, resolutions
model_name =['ERA5_','HadGEM3-GC31_','CNRM-CM6_','EC-Earth3P_','ECMWF-IFS_','CMCC-CM2_', 'MPI-ESM1-2_','FOCI_']
ERA5_variable = ['sea_surface_temperature','t2m','z300','precipitation_total','surface_latent_heat_flux']
FOCI_variable = ['tsw_', 'temp2_','geopoth_pl_300hPa_','precip_','ahfl_']
variable = ['tos_', 'tas_','zg_','pr_','hfls_']
resolution = [[''],
              ['-LL', '-MM', '-HH'],
              ['-1' , '-1-HR'],
              ['', '-HR'],
              ['-LR', '-MR','-HR'],
              ['-HR4','-VHR4'],
              ['-HR', '-XR','-ER', '-HR-PP'],
              ['-Standard' ,'-VIKING10']]

var_name = [['var34', 'var167','z','var228','var147'],['tos', 'tas','zg23','pr','hfls'],['tos', 'tas','zg10','pr','hfls'] ,['tos', 'tas','zg','pr','hfls'] ,['tos', 'tas','zg','pr','hfls'] ,['tos', 'tas','zg','pr','hfls'],['tos', 'tas','zg','pr','hfls'],['tsw', 'temp2','geopoth','precip','ahfl']]# 'slp_','precip_', 'ahfl_', 'ahfs_']
unit_factor= [[1,1, (1/9.81), (1000*24), (-1)], # ERA5
              [1,1, 1,        (3600*24), (-1)],     # HadGEM
              [1,1, 1,        (3600*24), (-1)],     # CNRM
              [1,1, 1,        (3600*24), (-1)],     # EC-Earth
              [1,1, 1,        (3600*24), (-1)],     # ECMWF-IFS
              [1,1, 1,        (3600*24), (-1)],     # CMCC
              [1,1, 1,        (3600*24), (-1)],     # MPI-ESM
              [1,1, 1,        (3600*24), (-1)],     # FOCI
              ]
HighRes_model={}
for mod in range(len(model_name)):
    print(model_name[mod])
    for var in range(len(variable)):
            for res in range(len(resolution[mod])):
                if mod==0:
                    #Input of ERA5
                    directory='/work/bm1398/m301111/Reanalysis/ERA5/daily/ERA5__'; file_end= '__1979-2019.nc'
                    HighRes_model['ERA5_'+variable[var]]=xr.open_dataset(directory+ERA5_variable[var]+file_end)

                if mod>0 and mod<7:
                    #Input of HighResMIP models
                    directory ='/work/bm1398/m301111/Models/HighResMIP/'+model_name[mod][:-1]+'/'+variable[var][:-1]+'/'+variable[var]+'daily/'+variable[var]+'day_'+model_name[mod][:-1]; #file_end='_control-1950_r1i1p1f1_gn_19500101-20501231_2_5_deg_06_00_shifted_inverted.nc'
                    file=glob.glob(directory+resolution[mod][res]+'_control*.nc')
                    #print(file)
                    HighRes_model[model_name[mod]+variable[var]+resolution[mod][res]] = xr.open_dataset(file[0])

                if mod==7:
                    #Input of FOCI
                    period ='1950-2049'
                    path =     '/work/bm1398/m301111/Models/FOCI/FOCI1.10-TM020_preind_no_nest/'+period+'/'
                    path_nest= '/work/bm1398/m301111/Models/FOCI/FOCI1.10-TM026_preind_VIKING10_nest/'+period+'/'
                    if res==0:
                        FOCI_file       = path+'FOCI1.10-TM020_echam6_echam_'+period+'_'+FOCI_variable[var][:-1]+'.nc'
                    if res==1:
                        FOCI_file   = path_nest+'FOCI1.10-TM026_echam6_echam_'+period+'_'+FOCI_variable[var][:-1]+'.nc'
                    HighRes_model[model_name[mod]+variable[var]+resolution[mod][res]] = xr.open_dataset(FOCI_file)
print('Loading data sets took', time.time()-start, 'seconds.')


def bootstrap_composites_map (var, sst_ds_Atl_dt_JJA, sample_num, sample_length, rm_value=5):

    ### Apply Bootstrap method with "sample_num" samples of length "sample_length" for composite of cold SST events with a negative tendency ###

    # Obtain a single sample with length "sample_length"
    sample = sst_ds_Atl_dt_JJA.sel(time=np.random.choice(sst_ds_Atl_dt_JJA.time, size=sample_length, replace=False))
    # Select data of randomly chosen time points
    var_boot_sample = var.where(var.time == sample.time).assign_coords(time=sample.time )
    # Compute sample mean over dim "time" with length "sample_length"
    boot_mean = var_boot_sample.mean('time')
    # Assign number to the sample mean
    boot_mean.coords['sample_number'] = sample_num+1
    #print(var_boot_sample)
    return boot_mean


start = time.time()

boots_mean_quant_025_map_lag0 = {}
boots_mean_quant_975_map_lag0 = {}

for mod in range(len(model_name)):
    boots_map={}
    boots_mean_map={}
    print(model_name[mod])
    ##### Deseasonalize SST, calculate NA average, detrend, select JJA values #####
    HighRes_model_tos_ds = {}                                                                       ### Used for step 1
    HighRes_model_tos_ds_Atl = {}                                                                   ### Used for step 2
    weights = np.cos(np.deg2rad(HighRes_model[model_name[mod]+variable[0]+resolution[mod][0]].lat)) ### Used for step 2
    HighRes_model_tos_ds_Atl_dt = {}                                                                ### Used for step 3
    HighRes_model_tos_ds_Atl_dt_JJA = {}                                                            ### Used for step 4
    HighRes_model_tos_ds_Atl_dt_JJA_mean = {}                                                       ### Used for step 5
    for res in range(len(resolution[mod])):
        print(resolution[mod][res])
        if mod==0:   # 1. Deseasonalize - ERA5
            A=xr.DataArray(np.full((92), 1979)) ### Used for step 5
            HighRes_model_tos_ds [model_name[mod]+variable[0]+'_ds'] = HighRes_model[model_name[mod]+variable[0]].var34.groupby("time.dayofyear") - (HighRes_model[model_name[mod]+variable[0]].var34.groupby("time.dayofyear").mean("time") - HighRes_model[model_name[mod]+variable[0]].var34.mean('time'))
        elif mod==1: # 1. Deseasonalize - HadGEM
            A=xr.DataArray(np.full((90), 1950)) ### Used for step 5
            HighRes_model_tos_ds [model_name[mod]+variable[0]+resolution[mod][res]+'_ds'] = HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tos.groupby("time.dayofyear") - (HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tos.groupby("time.dayofyear").mean("time") - HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tos.mean('time'))
        elif mod>1 and mod<7:  # 1. Deseasonalize - remaining HighResMIP models
            A=xr.DataArray(np.full((92), 1950)) ### Used for step 5
            HighRes_model_tos_ds [model_name[mod]+variable[0]+resolution[mod][res]+'_ds'] = HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tos.groupby("time.dayofyear") - (HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tos.groupby("time.dayofyear").mean("time") - HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tos.mean('time'))
        elif mod==7: # 1. Deseasonalize FOCI
            A=xr.DataArray(np.full((92), 1979)) ### Used for step 5
            HighRes_model_tos_ds [model_name[mod]+variable[0]+resolution[mod][res]+'_ds'] = HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tsw.groupby("time.dayofyear") - (HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tsw.groupby("time.dayofyear").mean("time") - HighRes_model[model_name[mod]+variable[0]+resolution[mod][res]].tsw.mean('time'))

        # 2. Calc NA SST average
        HighRes_model_tos_ds_Atl [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl'] = HighRes_model_tos_ds[model_name[mod]+variable[0]+resolution[mod][res]+'_ds'].sel(lon=slice(320.0, 345.0), lat=slice(60.0, 45.0)).weighted(weights).mean('lon', skipna=True).mean('lat', skipna=True)
        # 3. Detrend
        HighRes_model_tos_ds_Atl_dt [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt'] = xr.DataArray(data=(signal.detrend(HighRes_model_tos_ds_Atl [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl'].squeeze(), axis=0)), dims=["time"], coords=dict(time=HighRes_model_tos_ds_Atl [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl'].time))
        # 4. Select JJA values
        HighRes_model_tos_ds_Atl_dt_JJA [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt_JJA'] = HighRes_model_tos_ds_Atl_dt [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt'].sel(time=HighRes_model_tos_ds_Atl_dt [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt']['time.season']==['JJA'])
        HighRes_model_tos_ds_Atl_dt_JJA_mean [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt_JJA_mean'] = HighRes_model_tos_ds_Atl_dt_JJA [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt_JJA'].resample(time="QS-DEC").mean().dropna(dim='time')
        # 5. Add coordinate 'Year' to the data set
        years=[]
        for i in range(len(HighRes_model_tos_ds_Atl_dt_JJA_mean [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt_JJA_mean'] )):
            years.append(A)
            A=A+1
        years_array = xr.concat(years, dim='dim_0')
        HighRes_model_tos_ds_Atl_dt_JJA [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt_JJA'] = HighRes_model_tos_ds_Atl_dt_JJA [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt_JJA'].assign_coords(year=('time',years_array.data))


    HighRes_model_null = {}
    HighRes_model_detrend = {}
    HighRes_model_dt_ds = {}
    HighRes_model_dt_ds_anom = {}
    for var in [3]:
        print(variable[var][:-1])
        for res in range(len(resolution[mod])):
            if var==0:
                HighRes_model_null[model_name[mod]+variable[var]+resolution[mod][res]+'_null'] = HighRes_model[model_name[mod]+variable[var]+resolution[mod][res]][var_name [mod][var]].fillna(HighRes_model[model_name[mod]+variable[1]+resolution[mod][res]][var_name [mod][1]])
            else:
                HighRes_model_null[model_name[mod]+variable[var]+resolution[mod][res]+'_null'] = HighRes_model[model_name[mod]+variable[var]+resolution[mod][res]][var_name [mod][var]]

        # 6 Bootstrap method
            sample_num= 1000
            sample_len= [[12],
                        [37,41,37],
                        [45,42],
                        [23,43],
                        [25,40,40],
                        [42,47],
                        [35,28,38,42],
                        [31,37]]
            pool = mp.Pool(processes=64)
            if mod==6 and res>1: # change time axis from float to datetime in MPI-ESM by overwriting time dimension in ER with ER tos time dim
                HighRes_model_null[model_name[mod]+variable[var]+resolution[mod][res]+'_null']['time'] = HighRes_model [model_name[mod]+variable[0]+resolution[mod][res]]['time']

            # Detrend and produce xarray
            HighRes_model_detrend [model_name[mod]+variable[var]+resolution[mod][res]+'_detrend'] =  xr.DataArray(data=signal.detrend(HighRes_model_null[model_name[mod]+variable[var]+resolution[mod][res]+'_null'], axis=0), dims=["time","lat","lon"], coords=dict(time=HighRes_model_null[model_name[mod]+variable[var]+resolution[mod][res]+'_null'].time, lat=HighRes_model_null[model_name[mod]+variable[var]+resolution[mod][res]+'_null'].lat, lon=HighRes_model_null[model_name[mod]+variable[var]+resolution[mod][res]+'_null'].lon))
            # Deseasonalize
            HighRes_model_dt_ds [model_name[mod]+variable[var]+resolution[mod][res]+'_dt_ds'] =  HighRes_model_detrend [model_name[mod]+variable[var]+resolution[mod][res]+'_detrend'].groupby("time.dayofyear") - (HighRes_model_detrend [model_name[mod]+variable[var]+resolution[mod][res]+'_detrend'].groupby("time.dayofyear").mean("time") -  HighRes_model_detrend [model_name[mod]+variable[var]+resolution[mod][res]+'_detrend'].mean('time'))
            # Subtract time mean
            HighRes_model_dt_ds_anom [model_name[mod]+variable[var]+resolution[mod][res]+'_dt_ds_anom'] = HighRes_model_dt_ds [model_name[mod]+variable[var]+resolution[mod][res]+'_dt_ds'] - HighRes_model_dt_ds [model_name[mod]+variable[var]+resolution[mod][res]+'_dt_ds'].mean('time')

            boots_map[model_name[mod]+variable[var]+resolution[mod][res]]= pool.starmap(bootstrap_composites_map, zip(repeat(HighRes_model_dt_ds_anom [model_name[mod]+variable[var]+resolution[mod][res]+'_dt_ds_anom']),
                                                                                                              repeat(HighRes_model_tos_ds_Atl_dt_JJA [model_name[mod]+variable[0]+resolution[mod][res]+'_ds_Atl_dt_JJA']),
                                                                                                              range(sample_num),
                                                                                                              repeat(sample_len[mod][res])))
            boots_mean_map[model_name[mod]+variable[var]+resolution[mod][res]]= xr.concat(boots_map[model_name[mod]+variable[var]+resolution[mod][res]], dim='sample_number')
            boots_mean_quant_025_map_lag0[model_name[mod]+variable[var]+resolution[mod][res]] = boots_mean_map[model_name[mod]+variable[var]+resolution[mod][res]].quantile(0.025,'sample_number')*unit_factor[mod][var]
            boots_mean_quant_975_map_lag0[model_name[mod]+variable[var]+resolution[mod][res]] = boots_mean_map[model_name[mod]+variable[var]+resolution[mod][res]].quantile(0.975,'sample_number')*unit_factor[mod][var]
            print(boots_mean_quant_025_map_lag0[model_name[mod]+variable[var]+resolution[mod][res]].shape)
            pool.close()    # Close Pool and let all the processes complete
            pool.join()     # postpones the execution of next line of code until all processes in the queue are done.
    print('It took', "{:.1f}".format((time.time()-start)/60), 'minutes.')

boots_mean_quant_025_map_lag0_ds = xr.Dataset(boots_mean_quant_025_map_lag0)
boots_mean_quant_975_map_lag0_ds = xr.Dataset(boots_mean_quant_975_map_lag0)

boots_mean_quant_025_map_lag0_ds.to_netcdf('/work/bm1398/m301111/Models/HighResMIP/boots_mean_quant_025_map_lag0_pr.nc')
boots_mean_quant_975_map_lag0_ds.to_netcdf('/work/bm1398/m301111/Models/HighResMIP/boots_mean_quant_975_map_lag0_pr.nc')
print('NetCDF Files of pr bootstrapping composites successfully saved')
