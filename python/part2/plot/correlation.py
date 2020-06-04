#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 21:05:53 2018

@author: Luke
"""


#==============================================================================
#SUMMARY
#==============================================================================


# 24 February 2020

# This script loops through annual fieldmeans (1981-2018) of each lake model's 
# icestart, iceend and icedur to compute correlation coefficients between 
# hist and pi chunks (first takes anomalies w.r.t 1981-1999)

# Plots histogram of correlation coefficients between hist and pi chunks
#   this version uses only mmmhist to correlate (instead of all hist)
#   this version also provides an option between pearson and spearman-rank cc
#   finally, this performs a statistical test

# This version (iceindex, v3) produces 3 correlation plots, one for each ice index variable
# Uses all pichunks available from CLM45, VIC-lake and SIMSTRAT (also takes their histmmm)


#==============================================================================
#IMPORT
#==============================================================================


import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm
import seaborn as sns


#==============================================================================
#FUNCTIONS
#==============================================================================


# reading in netCDF files

def reader(filename):
    ds = xr.open_dataset(filename, decode_times=False)
    
    ds = ds.squeeze(drop=True)
    
    # open based on ice index
    if 'icestart' in filename:
        da = ds.icestart
    elif 'iceend' in filename:
        da = ds.iceend
    elif 'icedur' in filename:
        da = ds.icedur
    elif 'watertemp' in filename and 'clm' in filename and not 'chunk' in filename:
        da = ds.watertemp.isel(time=slice(0,-1))
    elif 'watertemp' in filename and 'albm' in filename and not 'chunk' in filename:
        da = ds.watertemp.isel(time=slice(0,-1))
    elif 'watertemp' in filename and 'chunk' in filename:
        da = ds.watertemp.isel(time=slice(0,-2))
    elif 'watertemp' in filename and 'simstrat' in filename and not 'chunk' in filename:
        da = ds.watertemp    # temporary until new fldmeans computed
        da = da.where(da<310).interpolate_na(dim='time',method='linear').isel(time=slice(0,-1))
    elif 'lmlt' in filename:
        da = ds.lmlt.isel(time=slice(0,-1))
        
    # temporal centering
    base = da.mean(dim='time').values
    da = da - base
    
    return da


# =============================================================================
# def reader(filename):
#     ds = xr.open_dataset(filename, decode_times=False)
#     ds = ds.squeeze(drop=True)
#     # open based on ice index
#     if 'icestart' in filename:
#         da = ds.icestart
#     elif 'iceend' in filename:
#         da = ds.iceend
#     elif 'icedur' in filename:
#         da = ds.icedur
#     elif 'watertemp' in filename and 'clm' in filename and not 'chunk' in filename:
#         #da = ds.watertemp.isel(time=slice(0,-1))
#         da = ds.watertemp
#     elif 'watertemp' in filename and 'albm' in filename and not 'chunk' in filename:
#         da = ds.watertemp
#     elif 'watertemp' in filename and 'chunk' in filename:
#         da = ds.watertemp.isel(time=slice(0,-1))
#     elif 'watertemp' in filename and 'simstrat' in filename and not 'chunk' in filename:
#         da = ds.watertemp    # temporary until new fldmeans computed
#         da = da.where(da<310).interpolate_na(dim='time',method='linear')
#     elif 'lmlt' in filename:
#         #da = ds.lmlt.isel(time=slice(0,-1))
#         da = ds.lmlt
#     # take centered data
#     base = da.mean(dim='time').values
#     da = da - base
#     return da
# =============================================================================

# ensemble math
def ensembler(data):
    concat_dim = np.arange(len(data))
    aligned = xr.concat(data,dim=concat_dim)
    mean = aligned.mean(dim='concat_dim')
    return mean


def correlation_proc(inDIR,piDIR,obsDIR):
    
    #==============================================================================
    #INITIALIZE DIRECTORIES + GRAB OBS
    #==============================================================================
    
    era5_watertemp_file = 'era5-land_lmlt_fldmean_1981_2018.nc'
    era5_icestart_file = 'era5-land_icestart_fldmean_1981_2018.nc'
    era5_iceend_file = 'era5-land_iceend_fldmean_1981_2018.nc'
    era5_icedur_file = 'era5-land_icedur_fldmean_1981_2018.nc'
    
    # set directory
    os.chdir(inDIR)

    #==============================================================================
    #DEFINE SETTINGS
    #==============================================================================
    
    # settings
    flag_prod=0; # 0: fldmean
                 # 1: sig
                 # 2: scaled_sig
                 # 3: eval
    
    flag_test=1;# 0: pearson correlation coefficient
                # 1: spearman (rank) correlation coefficient
        
    # list of models (filename format)
    models = ['clm45','albm','simstrat-uog','vic-lake','lakemod'] # temporarily removed "lake"
    
    # list of end variables
    endvariables = ['watertemp','icestart', 'iceend', 'icedur']
    
    # list of products
    products = ['fldmean','sig','scaled_sig','eval']
    
    # list of statistical tests
    test_types = ['pearson', 'spearman']
    
    # settings
    prod = products[flag_prod]
    testtype = test_types[flag_test]
    
    #==============================================================================
    #DEFINE PLOT SETTINGS
    #==============================================================================    
    
    #========== FONT ==========#
    
    title_font = 14
    
    tick_font = 12
    
    #========== LEGEND ==========#
    
    # legend font size
    legend_font = 8
    
    # bbox
    le_x0 = 0.9
    le_y0 = 0.6
    le_xlen = 0.15
    le_ylen = 0.25
    
    # space between entries
    legend_entrypad = 0.5
    
    # length per entry
    legend_entrylen = 0.75
    
    #==============================================================================
    #READ IN DATA
    #==============================================================================
    
    hist_files = {}
    pi_files = {}
    
    hist_data = {}
    pichunks = {}
    
    # access all historical+rcp evaluation fldmean files
    for var in endvariables:
        hist_files[var] = []
        for file in [file for file in sorted(os.listdir(inDIR))\
                 if var in file and prod in file]:
            hist_files[var].append(file)
    
    # open historical data
    for var in endvariables:
        hist_data[var] = []
        for file in hist_files[var]:
            hist_data[var].append(reader(file))
    
    # filter model names iterable for existing model data for plotting
    available_models = []
    for mod in models:
        count = 0
        for file in hist_files['icestart']:
            if mod in file:
                count = count+1
                if count == 1:
                    available_models.append(mod)
          
    # grab relevant picontrol files for     
    os.chdir(piDIR) 
    for var in endvariables:
        pi_files[var] = []
        for file in [file for file in sorted(os.listdir(piDIR))\
                 if var in file and\
                         'chunk' in file and\
                             'picontrol' in file]:
            pi_files[var].append(file)  
            
    # mean of historical ensemble
    hist_mmm = {}
    for var in endvariables:
        hist_mmm[var] = ensembler(hist_data[var]).values 
        # read in pichunks
        pichunks[var] = []
        for file in pi_files[var]:
            pichunks[var].append(reader(file).values)
    
    # manual insert of era5-land reanalysis data
    os.chdir(obsDIR)
    era5_obs = {}
    era5_obs['watertemp'] = reader(era5_watertemp_file)
    era5_obs['icestart'] = reader(era5_icestart_file)
    era5_obs['iceend'] = reader(era5_iceend_file)
    era5_obs['icedur'] = reader(era5_icedur_file)
    
    #==============================================================================
    #CORRELATION SECTION
    #==============================================================================
    
    # using mean of historical ("pcc" -> pearson cor coeff, "spcc" -> spearman)
    pi_histmmm_pcc = {}
    pi_histmmm_spcc = {}
    
    histmmm_obs_pcc = {}
    histmmm_obs_spcc = {}
    
    samples = {}
    mean = {}
    std = {}
    n = {}
    
    cc_99 = {}
    cc_95 = {}
    cc_90 = {}
    
    for var in endvariables:
        
        pi_histmmm_pcc[var] = []
        pi_histmmm_spcc[var] = []
        
        for chunk in pichunks[var]:
            pi_histmmm_pcc[var].append(stats.pearsonr(chunk,hist_mmm[var])[0])
            pi_histmmm_spcc[var].append(stats.spearmanr(chunk,hist_mmm[var])[0])
    
        histmmm_obs_pcc[var] = stats.pearsonr(hist_mmm[var],era5_obs[var].values)[0]
        histmmm_obs_spcc[var] = stats.spearmanr(hist_mmm[var],era5_obs[var].values)[0]
    
        # 90, 95 % confidence level val for Ho: mean = 0 and  (one tailed)
        if testtype == 'pearson':
            samples[var] = pi_histmmm_pcc[var]
        elif testtype == 'spearman':
            samples[var] = pi_histmmm_spcc[var]
    
        mean[var] = np.mean(samples[var])
        std[var] = np.std(samples[var])
        n[var] = len(samples[var])
    
        cc_99[var] = stats.norm.ppf(q=0.99,loc=mean[var],scale=std[var])
        cc_95[var] = stats.norm.ppf(q=0.95,loc=mean[var],scale=std[var])
        cc_90[var] = stats.norm.ppf(q=0.90,loc=mean[var],scale=std[var])
        
    return samples,mean,n,std,\
           histmmm_obs_pcc,histmmm_obs_spcc,\
           pi_histmmm_pcc,pi_histmmm_spcc,\
           cc_99,cc_95,cc_90
         
    
 