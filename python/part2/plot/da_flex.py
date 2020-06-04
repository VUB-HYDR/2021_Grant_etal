#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 09:51:03 2020

@author: luke
"""


#%%==============================================================================
# IMPORT
#==============================================================================

import numpy as np
import xarray as xr
from scipy import signal, stats
import os
import sys
import subprocess
import PyDnA as pda
import scipy.linalg as spla
import scipy.stats as sps
from random import shuffle

#==============================================================================
# FUNCTION
#==============================================================================

def data_loader(inDIR,piDIR,obsDIR,var,fp_scen,window,block):
    
    
# =============================================================================
#     inDIR -> directory with RCP fp data
#     piDIR -> directory with ctl data (also to be used as fp if chosen)
#     obsDIR -> directory with observations (y)
#     var -> watertemp, icestart, iceend or icedur choice for scaling factors
#     fp_scen -> rcp85,rcp60,ctl options as strings "rcp85", "rcp60" and "ctl" for scaling factor
# =============================================================================
    
    
    fp_files = []
    pi_files = []
    fp_data = []
    pi_data = []
    

    # access all historical+,8.5 files separate
    os.chdir(inDIR)
    for file in [file for file in sorted(os.listdir(inDIR))\
             if var in file and fp_scen in file]:
            fp_files.append(file)
            
    for file in fp_files:
        fp_data.append(pda.reader(file,window,block)) 
        
    fp_mmm = pda.ensembler(fp_data)
    
    # access all ctrl files
    os.chdir(piDIR)
    for file in [file for file in sorted(os.listdir(piDIR))\
             if var in file]:
        pi_files.append(file)
    shuffle(pi_files)
        
    for file in pi_files:
        pi_data.append(pda.reader(file,window,block)) 
    
    
    # access obs
    os.chdir(obsDIR)
    era5_watertemp_file = 'era5-land_lmlt_fldmean_1981_2018.nc'
    era5_icestart_file = 'era5-land_icestart_fldmean_1981_2018.nc'
    era5_iceend_file = 'era5-land_iceend_fldmean_1981_2018.nc'
    era5_icedur_file = 'era5-land_icedur_fldmean_1981_2018.nc'
    
    if var == 'watertemp':
        obs = pda.reader(obsDIR+'/'+era5_watertemp_file,window,block)
    elif var == 'icestart':
        obs = pda.reader(obsDIR+'/'+era5_icestart_file,window,block)
    elif var == 'iceend':
        obs = pda.reader(obsDIR+'/'+era5_iceend_file,window,block)
    elif var == 'icedur':
        obs = pda.reader(obsDIR+'/'+era5_icedur_file,window,block)
    
    
    # data management for da
    obs = obs.values
    if fp_scen == 'rcp60' or fp_scen == 'rcp85':
        nx = np.array(([len(fp_data)]))
    elif fp_scen == 'ctl':
        nx = np.array(([len(pi_data)]))
    ctl_list = []
    for array in pi_data:
        ctl_list.append(array.values)
    ctl = np.stack(ctl_list,axis=0)
    fp = np.stack([fp_mmm.values],axis=0)
    
    
    return obs,fp,ctl,nx

    






























      

