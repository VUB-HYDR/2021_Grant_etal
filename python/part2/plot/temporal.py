#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 21:05:53 2018

@author: Luke
"""


#==============================================================================
#SUMMARY
#==============================================================================


# 13 March 2020

# functions to be called for ice index proc+plotting


#==============================================================================
#IMPORT
#==============================================================================


import os
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import math
from random import shuffle


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
    
    # take anomaly w.r.t. 1981-1999
    base = da.isel(time=slice(0,18)).mean(dim='time').values
    da = da - base
    
    return da.rolling(time=5, center=False).mean()


# =============================================================================
# def reader(filename,window):
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
#     # take anomaly w.r.t. 1981-1999
#     base = da.isel(time=slice(0,18)).mean(dim='time').values
#     da = da - base
#     
#     return da.rolling(time=window, center=False).mean().dropna("time")
# =============================================================================

# ensemble math
def ensembler(data):
    concat_dim = np.arange(len(data))
    aligned = xr.concat(data,dim=concat_dim)
    mean = aligned.mean(dim='concat_dim')
    std = aligned.std(dim='concat_dim')
    se = std/math.sqrt(len(concat_dim))
    return [mean,std,se]


def temporal_proc(inDIR,piDIR,obsDIR,shuff,window):
        
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
     
    flag_prod=0; # 0: fldmean
                 # 1: sig
                 # 2: scaled_sig
                 # 3: eval    
    
    # list of products
        # list of products
    products = ['fldmean','sig','scaled_sig','eval']
    endvariables = ['watertemp','icestart','iceend','icedur']
    models = ['clm45','albm','simstrat-uog','vic-lake']
    
    # settings
    prod = products[flag_prod]    
    
    #==============================================================================
    #READ IN DATA
    #==============================================================================
    
    hist_files = {}
    pi_files = {}
    
    hist_data = {}
    pi_data = {}
    
    hist_mmm = {}
    pi_mmm = {}
    
    # access all historical+rcp evaluation fldmean files
    for var in endvariables:
        hist_files[var] = []
        for file in [file for file in sorted(os.listdir(inDIR))\
                 if var in file and prod in file and 'rcp85' in file]:
            hist_files[var].append(file)
        hist_data[var] = []
        for file in hist_files[var]:
            hist_data[var].append(reader(file)) 
        hist_mmm[var] = ensembler(hist_data[var])
          
    ice_models = []
    for mod in models:
        count = 0
        for file in hist_files['icedur']:
            if mod in file:
                count = count+1
                if count == 1:
                    ice_models.append(mod)
                    
    watertemp_models = []
    for mod in models:
        count = 0
        for file in hist_files['watertemp']:
            if mod in file:
                count = count+1
                if count == 1:
                    watertemp_models.append(mod)
        
    # grab relevant picontrol files
    os.chdir(piDIR) 
    
    for var in endvariables[1:]:    # match ice index endvars with icemodels
        pi_files[var] = []
        for mod in ice_models:
            modlist = []
            for file in [file for file in sorted(os.listdir(piDIR))\
                     if var in file and\
                         mod in file and\
                             'chunk' in file and\
                                 'picontrol' in file]:
                modlist.append(file)
            if shuff == 'yes': # use 8 random pi files
                shuffle(modlist)
                i = 0
                while i < 8:   
                    pi_files[var].append(modlist[i]) 
                    i += 1
            elif shuff == 'no': # use all pi files
                for file in modlist: 
                    pi_files[var].append(file)
        
        pi_data[var] = []
        for file in pi_files[var]:
            pi_data[var].append(reader(file))
        pi_mmm[var] = ensembler(pi_data[var])  
        
    var = 'watertemp'
    pi_files[var] = []
    for mod in watertemp_models:
        modlist = []
        for file in [file for file in sorted(os.listdir(piDIR))\
                 if var in file and\
                     mod in file and\
                         'chunk' in file and\
                             'picontrol' in file]:
            modlist.append(file)
        if shuff == 'yes': # use 8 random pi files
            shuffle(modlist)
            i = 0
            while i < 8:   
                pi_files[var].append(modlist[i]) 
                i += 1
        elif shuff == 'no': # use all pi files
            for file in modlist: 
                pi_files[var].append(file)
    pi_data[var] = []
    for file in pi_files[var]:
        pi_data[var].append(reader(file)) 
    pi_mmm[var] = ensembler(pi_data[var])   
        
    # manual insert of era5-land reanalysis data
    os.chdir(obsDIR)
    era5_obs = {}
    era5_obs['watertemp'] = reader(era5_watertemp_file)
    era5_obs['icestart'] = reader(era5_icestart_file)
    era5_obs['iceend'] = reader(era5_iceend_file)
    era5_obs['icedur'] = reader(era5_icedur_file)
    

            
    return hist_mmm,pi_mmm,era5_obs

 