#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 21:05:53 2018

@author: Luke
"""


#==============================================================================
#SUMMARY
#==============================================================================


# 12 March 2020

# This script processes data for an icedur scaling plot


#==============================================================================
#IMPORT
#==============================================================================


import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D


#==============================================================================
#FUNCTIONS
#==============================================================================


# reading in netCDF files
def reader(filename):
    ds = xr.open_dataset(filename, decode_times=False)
    
    # open based on ice index
    if 'icestart' in filename:
        da = ds.icestart
    if 'iceend' in filename:
        da = ds.iceend
    if 'icedur' in filename:
        da = ds.icedur
    
    da = da.squeeze()
    return da

#data ensembler
def ensembler(data,roller_dim):
    concat_dim = np.arange(len(data))
    aligned = xr.concat(data,dim=concat_dim)    
    mean = aligned.mean(dim='concat_dim').rolling(dim={roller_dim:21},center=True).mean()
    scen_max = aligned.max(dim='concat_dim').rolling(dim={roller_dim:21},center=True).mean()
    scen_min = aligned.min(dim='concat_dim').rolling(dim={roller_dim:21},center=True).mean()
    return [mean,scen_max,scen_min]

# reader for air temperature GCM data
def stringer(folder,special,special85):
    airfiles = sorted([folder + '/' + file for file in os.listdir(folder) if not 'runmean' in file])
    if folder+'/'+'.DS_Store' in airfiles: airfiles.remove(folder+'/'+'.DS_Store')
    fullpi = np.loadtxt(airfiles[1])[:,1][:-200]
    base = np.mean(fullpi,axis=0)
    historical = np.loadtxt(airfiles[0])[:,1]
    pi = np.loadtxt(airfiles[1])[:,1][:200]
    
    #special is for files without extended projections until 2299
    if special is True:                     
        rcp26 = np.loadtxt(airfiles[2])[:,1]
    if special is False:
        rcp26 = np.loadtxt(airfiles[2])[:,1][:-200]
    rcp60 = np.loadtxt(airfiles[3])[:,1]
    if special85 is True:
        rcp85 = np.loadtxt(airfiles[4])[:,1]
    if special85 is False:
        rcp85 = np.loadtxt(airfiles[4])[:,1][:-200]
        
    rcp_26 = (np.concatenate((pi,historical,rcp26),axis=0)-base)
    rcp_60 = (np.concatenate((pi,historical,rcp60),axis=0)-base)
    rcp_85 = (np.concatenate((pi,historical,rcp85),axis=0)-base)
    
    return xr.DataArray(rcp_26),xr.DataArray(rcp_60),xr.DataArray(rcp_85)


# processing function
def iiscale_proc(inDIR,airDIR):
    
    # set directory
    os.chdir(inDIR)
    
    #==============================================================================
    #DEFINE SETTINGS
    #==============================================================================
                    
    flag_var=5;  # 0: watertemp
                 # 1: lakeicefrac
                 # 2: icethick
                 # 3: icestart
                 # 4: iceend
                 # 5: icedur
     
    flag_prod=0; # 0: fldmean
                 # 1: sig
                 # 2: scaled_sig
                 # 3: eval        
    
    # list of models (filename format)
    models = ['clm45','albm','simstrat-uog','vic-lake','lakemod'] # temporarily removed "lake"
    
    # list of variables (1st root variables, then processed)
    variables = ['watertemp','lakeicefrac','icethick','icestart','iceend','icedur']
    
    # list of forcings
    forcings = ['gfdl-esm2m','hadgem2-es','ipsl-cm5a-lr','miroc5']
    
    # list of products
    products = ['fldmean','sig','scaled_sig','eval']
    
    # list of rcps (no flagging)
    rcps = ['rcp26','rcp60','rcp85']
    
    # settings
    var = variables[flag_var]
    prod = products[flag_prod]
    
    #==============================================================================
    #READ IN DATA
    #==============================================================================
    
    files = []
    
    # access all icedur fldmean files
    for file in [file for file in sorted(os.listdir(inDIR))\
                 if var in file and prod in file]:
        files.append(file)
    
    filesets = {}
    datasets = {}
    #[rcp][gcm] to ensemble<stat> on gcm variants
    for rcp in rcps:
        filesets[rcp] = {}
        datasets[rcp] = []
        for gcm in forcings:
            filesets[rcp][gcm] = []
            for mod in models:
                for file in files:
                    if mod in file and gcm in file and rcp in file:
                        filesets[rcp][gcm].append(file)
                    if mod in file and gcm in file and 'picontrol' in file:
                        filesets[rcp][gcm].append(file)
                    if mod in file and gcm in file and 'historical' in file:
                        filesets[rcp][gcm].append(file)
                        
    existingmods=[]
    for mod in models:
        count = 0
        for file in files:
            if mod in file: 
                count = count+1
                if count == 1:
                    existingmods.append(mod)
                           
    for mod in existingmods:
        for rcp in rcps:
            for gcm in forcings:
                templist_a = [] # list for pi+hist+rcp under same mod and gcm
                templist_b = [] # list for pi_pi+pi_hist+pi_rcp under same mod and gcm
                for file in filesets[rcp][gcm]:
                    if mod in file and 'picontrol' in file:
                        templist_b.append(reader(file))
                    if mod in file and 'picontrol' in file and '1661_1860' in file:
                        templist_a.append(reader(file))
                    if mod in file and 'historical' in file and '1861_2005' in file:
                        templist_a.append(reader(file))
                    if mod in file and rcp in file and '2006_2099' in file:
                        templist_a.append(reader(file))
                if len(templist_b)>1 and len(templist_a)>1: # clm45 and simstrat
                    b = xr.concat(templist_b,dim='time').mean(dim='time').values
                    a = xr.concat(templist_a,dim='time')
                elif len(templist_b)==1 and len(templist_a)==1: # albm case
                    b = templist_b[0].mean(dim='time').values
                    a = templist_a[0]
                elif len(templist_b)==1 and len(templist_a)>1: # vic-lake case
                    b = templist_b[0].mean(dim='time').values
                    a = xr.concat(templist_a,dim='time')
                c = a-b
                datasets[rcp].append(c)
            
    ensembles = {}
    for rcp in rcps:
        ensembles[rcp] = ensembler(data=datasets[rcp],roller_dim='time')
    
    #==============================================================================
    #READ IN AIR TEMPERATURES
    #==============================================================================
    
    # gfdl-esm2m
    gfdl_26,gfdl_60,gfdl_85 = stringer(folder=airDIR+'/GFDL-ESM2M',\
                                       special=True,special85=True)
    
    # hadgem2-es
    hadgem_26,hadgem_60,hadgem_85 = stringer(folder=airDIR+'/HadGEM2-ES',\
                                             special=False,special85=True)
    
    # ipsl-cm5a-lr
    ipsl_26,ipsl_60,ipsl_85 = stringer(folder=airDIR+'/IPSL-CM5A-LR/',\
                                       special=False,special85=False)
    
    # miroc5
    miroc_26,miroc_60,miroc_85 = stringer(folder=airDIR+'/MIROC5/',\
                                          special=False,special85=True)
    
    # load air arrays for ensembling
    rcp26_air = [gfdl_26,hadgem_26,ipsl_26,miroc_26]
    rcp60_air = [gfdl_60,hadgem_60,ipsl_60,miroc_60]
    rcp85_air = [gfdl_85,hadgem_85,ipsl_85,miroc_85]
    
    # ensemble stats of air data
    air_ensembles = {}
    air_ensembles['rcp26'] = ensembler(data=rcp26_air,roller_dim='dim_0')
    air_ensembles['rcp60'] = ensembler(data=rcp60_air,roller_dim='dim_0')
    air_ensembles['rcp85'] = ensembler(data=rcp85_air,roller_dim='dim_0')
    
    #==============================================================================
    #IMPACT YEARS
    #==============================================================================
    
    #gather impact year plotting points
    impact_yrs = {}

    impact_yrs['rcp26_x2030'] = air_ensembles['rcp26'][0][368].values
    impact_yrs['rcp26_y2030'] = ensembles['rcp26'][0][368].values
    impact_yrs['rcp26_x2050'] = air_ensembles['rcp26'][0][388].values
    impact_yrs['rcp26_y2050'] = ensembles['rcp26'][0][388].values
    impact_yrs['rcp26_x2100'] = air_ensembles['rcp26'][0][425].values
    impact_yrs['rcp26_y2100'] = ensembles['rcp26'][0][425].values
    
    impact_yrs['rcp60_x2030'] = air_ensembles['rcp60'][0][368].values
    impact_yrs['rcp60_y2030'] = ensembles['rcp60'][0][368].values
    impact_yrs['rcp60_x2050'] = air_ensembles['rcp60'][0][388].values
    impact_yrs['rcp60_y2050'] = ensembles['rcp60'][0][388].values
    impact_yrs['rcp60_x2100'] = air_ensembles['rcp60'][0][425].values
    impact_yrs['rcp60_y2100'] = ensembles['rcp60'][0][425].values
    
    impact_yrs['rcp85_x1980'] = air_ensembles['rcp85'][0][318].values
    impact_yrs['rcp85_y1980'] = ensembles['rcp85'][0][318].values
    impact_yrs['rcp85_x2000'] = air_ensembles['rcp85'][0][338].values
    impact_yrs['rcp85_y2000'] = ensembles['rcp85'][0][338].values
    impact_yrs['rcp85_x2020'] = air_ensembles['rcp85'][0][358].values
    impact_yrs['rcp85_y2020'] = ensembles['rcp85'][0][358].values
    impact_yrs['rcp85_x2030'] = air_ensembles['rcp85'][0][368].values
    impact_yrs['rcp85_y2030'] = ensembles['rcp85'][0][368].values
    impact_yrs['rcp85_x2050'] = air_ensembles['rcp85'][0][388].values
    impact_yrs['rcp85_y2050'] = ensembles['rcp85'][0][388].values
    impact_yrs['rcp85_x2100'] = air_ensembles['rcp85'][0][425].values
    impact_yrs['rcp85_y2100'] = ensembles['rcp85'][0][425].values
    
    return air_ensembles,ensembles,impact_yrs
   