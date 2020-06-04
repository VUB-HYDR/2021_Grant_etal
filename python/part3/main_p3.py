#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 08:15:24 2020

@author: Luke
"""

#==============================================================================
# SUMMARY
#==============================================================================


# 18 May 2020

# main script for calling maps of climate change impacts (f3)


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
inDIR = os.path.join(curDIR, 'data/')

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
             
flag_main_SI=1;   # 0: process main script figures
                  # 1: process SI figures
               
flag_SI_var=0;    # 0: watertemp
                  # 1: icethick
                  # 2: iceindex

# this option only applies to the choice of rcp for the main text fig    
  
flag_rcp=2;       # 0: rcp26
                  # 1: rcp60
                  # 2: rcp85
                  
flag_prod=0;      # 0: sig
                  # 1: scaled_sig

# res of SI pics (dpi for dots per inch)

flag_dpi=1;   # 0: 100
              # 1: 200
              # 2: 300
              # 3: 400
              # 4: 500
              # 5: 600


rcps = ['rcp26','rcp60','rcp85']
variables = ['watertemp','icethick','iceindex']
products = ['sig','scaled_sig']
dpis = [100,200,300,400,500,600]


rcp = rcps[flag_rcp]
var = variables[flag_SI_var]
prod = products[flag_prod]
dpi = dpis[flag_dpi]


#==============================================================================
# PROC + PLOT
#==============================================================================


sys.path.append(plotDIR)

# conditions for updating main text/SI plots        
if flag_main_SI == 0:
    
    # watertemp annual extraction
    from watertemp_sig_annual import *
    annual_data,lat,lon,ann_mods = wtsig_proc(inDIR)
    wt_annual = annual_data[rcp][0] 
    
    # watertemp jja scaled extraction
    from watertemp_scaled_sig import *
    scale_data,lat,lon,sca_mods = wtscalesig_proc(inDIR)
    wt_jja_scaled = scale_data[rcp][2] 

    # iceindex signal extraction here
    from iceindex_sig import *
    iisig_data,lat,lon = iisig_proc(inDIR)
    iisig_start = iisig_data[rcp][0]
    iisig_end = iisig_data[rcp][1]
    iisig_dur = iisig_data[rcp][2]
    iisig_data = [iisig_start,iisig_end,iisig_dur]
    
    # run plot
    from plot_p3 import *
    plot_p3(wt_annual,wt_jja_scaled,iisig_data,lat,lon,outDIR,flag_svplt)
    
        
elif flag_main_SI == 1:
    
    # water temperature
    if flag_SI_var == 0:
        
        if flag_prod == 0:
            
            from watertemp_sig import *
            multimodel,lat,lon,seas_mods=wtsig_proc(inDIR)
            wtsig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi)
            del wtsig_plot # new function for annual below
            del wtsig_proc
        
            from watertemp_sig_annual import *
            multimodel,lat,lon,ann_mods=wtsig_proc(inDIR)
            wtsig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi)
        
        if flag_prod == 1:
            
            from watertemp_scaled_sig import *
            multimodel,lat,lon,sca_mods=wtscalesig_proc(inDIR)
            wtscalesig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi)
    
    # ice thickness    
    elif flag_SI_var == 1:
        
        if flag_prod == 0:
            
            from icethick_sig import *
            multimodel,lat,lon,arrays=itsig_proc(inDIR)
            itsig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi)
        
        if flag_prod == 1:
            
            from icethick_scaled_sig import *
            multimodel,lat,lon=itscalesig_proc(inDIR)
            itscalesig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi)
    
    # ice index    
    elif flag_SI_var == 2:
        
        if flag_prod == 0:
        
            from iceindex_sig import *
            multimodel,lat,lon,lakemod_85=iisig_proc(inDIR)
            iisig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi)
            del reader # new reader function for scaled iceindex
        
        if flag_prod == 1:
            
            from iceindex_scaled_sig import *
            multimodel,lat,lon=iiscalesig_proc(inDIR)
            iiscalesig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi)
        
    

    




