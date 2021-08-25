#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 08:15:24 2020

@author: Luke
"""

#==============================================================================
# SUMMARY
#==============================================================================


# main script for plotting global mean anomalies (f4)


#==============================================================================
# IMPORT
#==============================================================================


import sys
import os
import pickle as pk


#==============================================================================
# PATH
#==============================================================================


# change current working directory to pathway of this file
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# object for current dir
curDIR = os.getcwd()

# data input directory
inDIR = os.path.join(curDIR, 'data/lakemodels/')
airDIR = os.path.join(curDIR, 'data/airtemps/')

# script directory
plotDIR = os.path.join(curDIR, 'plot/')

# figure output directory
outDIR = os.path.join(curDIR, 'figures/final')


#==============================================================================
# OPTIONS
#==============================================================================


# adjust these settings for refreshing plots based on sig vs scaledsig, variable  
# and main text versus SI

flag_svplt=0;     # 0: do not save plot
                  # 1: save plot in outDIR
             
flag_main_SI=0;   # 0: process main script figures
                  # 1: process SI figures
               
flag_SI_var=0;    # 0: watertemp
                  # 1: icethick
                  # 2: icedur
                  
flag_prod=0;      # 0: timeseries
                  # 1: scaling

# res of SI pics (dpi for dots per inch)
flag_dpi=2;   # 0: 100
              # 1: 200
              # 2: 300
              # 3: 400
              # 4: 500
              # 5: 600

# variables
variables = ['watertemp','icethick','iceindex']

# products
products = ['timeseries','scaling']
dpis = [100,200,300,400,500,600]

# assert settings
var = variables[flag_SI_var]
prod = products[flag_prod]
dpi = dpis[flag_dpi]

#==============================================================================
# PROC + PLOT
#==============================================================================


sys.path.append(plotDIR)

# conditions for updating main text/SI plots        
if flag_main_SI == 0:

    os.chdir(inDIR)
    
    if not os.path.isfile('pk_data_p4.pkl'):
    
        print('no pickled data')
        
        # watertemp timeseries extraction
        from watertemp_timeseries import *
        wt_mmm,keys_test = wt_proc(inDIR,flag_main_SI)
        
        # icethick timeseries extraction
        from icethick_timeseries import *
        it_mmm = it_proc(inDIR,flag_main_SI)
        
        # iceindex timeseries extraction
        from icedur_timeseries import *
        ii_mmm,arrays_test = ii_proc(inDIR,flag_main_SI)
        
        # watertemp scaling extraction
        from watertemp_scaling import *
        wt_air_ens,wt_ens,wt_iy,wt_test = wtscale_proc(inDIR,airDIR)
        
        # icethick scaling extraction
        from icethick_scaling import *
        it_air_ens,it_ens,it_iy = itscale_proc(inDIR,airDIR)
        
        # iceindex scaling extraction
        from icedur_scaling import *
        ii_air_ens,ii_ens,ii_iy = iiscale_proc(inDIR,airDIR)
        
        pickle_data = {}
        
        pickle_data['wt_mmm'] = wt_mmm
        pickle_data['it_mmm'] = it_mmm
        pickle_data['ii_mmm'] = ii_mmm
        
        pickle_data['wt_air_ens'] = wt_air_ens
        pickle_data['wt_ens'] = wt_ens
        pickle_data['wt_iy'] = wt_iy
        
        pickle_data['it_air_ens'] = it_air_ens
        pickle_data['it_ens'] = it_ens
        pickle_data['it_iy'] = it_iy        

        pickle_data['ii_air_ens'] = ii_air_ens
        pickle_data['ii_ens'] = ii_ens
        pickle_data['ii_iy'] = ii_iy        
        
        output = open('pk_data_p4.pkl','wb')
        pk.dump(pickle_data,output)
        output.close()
        
    elif os.path.isfile('pk_data_p4.pkl'):
        
        print('use pickled data')
        
        pkl_file = open('pk_data_p4.pkl','rb')
        pickle_data = pk.load(pkl_file)
        pkl_file.close()
        
        wt_mmm = pickle_data['wt_mmm']
        it_mmm = pickle_data['it_mmm']
        ii_mmm = pickle_data['ii_mmm']
        
        wt_air_ens = pickle_data['wt_air_ens']
        wt_ens = pickle_data['wt_ens']
        wt_iy = pickle_data['wt_iy']
        
        it_air_ens = pickle_data['it_air_ens']
        it_ens = pickle_data['it_ens']
        it_iy = pickle_data['it_iy']

        ii_air_ens = pickle_data['ii_air_ens']
        ii_ens = pickle_data['ii_ens']
        ii_iy = pickle_data['ii_iy']


    # run plot
    from plot_p4 import *
    plot_p4(wt_mmm,it_mmm,ii_mmm,\
            wt_air_ens,wt_ens,wt_iy,\
            it_air_ens,it_ens,it_iy,\
            ii_air_ens,ii_ens,ii_iy,\
            outDIR,flag_svplt,dpi)
    

elif flag_main_SI == 1:
    
    # water temperature
    if flag_SI_var == 0:
        
        if flag_prod == 0:
            
            from watertemp_timeseries import *
            datasets = wt_proc(inDIR,flag_main_SI)
            wt_plot(datasets,outDIR,var,prod,flag_svplt,dpi)
    
    # ice thickness    
    elif flag_SI_var == 1:
        
        if flag_prod == 0:
            
            from icethick_timeseries import *
            datasets,filesets,datasets_test=it_proc(inDIR,flag_main_SI)
            it_plot(datasets,outDIR,var,prod,flag_svplt,dpi)
    
    # ice index    
    elif flag_SI_var == 2:
        
        if flag_prod == 0:
        
            from icedur_timeseries import *
            datasets=ii_proc(inDIR,flag_main_SI)
            ii_plot(datasets,outDIR,var,prod,flag_svplt,dpi)
            
        

