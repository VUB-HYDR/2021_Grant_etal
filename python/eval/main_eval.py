#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 08:15:24 2020

@author: Luke
"""

#==============================================================================
# SUMMARY
#==============================================================================


# 17 March 2020

# this serves as the wrapper script for processing + plotting SI evaluation figures
# for ice index and watertemp


#==============================================================================
# IMPORT
#==============================================================================


import sys
import os


#==============================================================================
# PATH
#==============================================================================


# change current working directory to pathway of this file
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# object for current dir
curDIR = os.getcwd()

# data input directory
inDIR = os.path.join(curDIR, 'data/lakemodels/')
obsDIR = os.path.join(curDIR, 'data/obs/')
valDIR = os.path.join(curDIR, 'data/validation/')

# script directory
plotDIR = os.path.join(curDIR, 'plot/')

# figure output directory
outDIR = os.path.join(curDIR, 'figures/')


#==============================================================================
# OPTIONS
#==============================================================================


# adjust these settings for refreshing plots based on sig vs scaledsig, variable  
# and main text versus SI

flag_svplt=1;     # 0: do not save plot
                  # 1: save plot in outDIR

flag_eval=0;      # 0: evaluation of isimip w/ era5-land
                  # 1: validation era5-land with globo/amsr

flag_var=1;    # 0: watertemp
               # 1: ice indexes                 


#==============================================================================
# PROC + PLOT
#==============================================================================


sys.path.append(plotDIR)

# eval of isimip    
if flag_eval == 0:     
    
    if flag_var == 0:
        
        from watertemp_eval import *
        available_models,era5,isimip,lat,lon = wteval_proc(inDIR,obsDIR)
                
        wteval_plot(outDIR,flag_svplt,\
                    available_models,\
                    era5,isimip,lat,lon)
    
    elif flag_var == 1:
        
        from iceindex_eval import *
        available_models,era5, multimodel,lat,lon = iieval_proc(inDIR,obsDIR)
        
        iieval_plot(outDIR,flag_svplt,\
                    available_models,\
                    era5,multimodel,lat,lon)
 
# compare era5-land w/ other prods       
elif flag_eval == 1:
    
    if flag_var == 0:
        
        from globolakes_bias import *
        globolakes_plot(valDIR,outDIR,flag_svplt)
        
    elif flag_var == 1:
        
        from amsr_bias import *
        amsr_plot(valDIR,outDIR,flag_svplt)
            
        
        
        


















