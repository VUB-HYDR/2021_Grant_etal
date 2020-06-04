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

# main script for calling plotting of reconstructed signals (f1)


#==============================================================================
# IMPORT
#==============================================================================


import sys
import os
import numpy as np


#==============================================================================
# PATH
#==============================================================================


# change current working directory to pathway of this file
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# object for current dir
curDIR = os.getcwd()

# data input directory
iceDIR = os.path.join(curDIR, 'data/icecover/')
mltDIR = os.path.join(curDIR, 'data/mlt/')

# script directory
plotDIR = os.path.join(curDIR, 'plot/')

# figure output directory
outDIR = os.path.join(curDIR, 'figures/')


#==============================================================================
# OPTIONS
#==============================================================================


# adjust these settings for either testing OF output or producing main text fig

flag_svplt=1;     # 0: do not save plot
                  # 1: save plot in outDIR

             
flag_var=1;   # 0: ice cover signal
              # 1: mlt signal

# res of SI pics (dpi for dots per inch)

flag_dpi=1;   # 0: 100
              # 1: 200
              # 2: 300
              # 3: 400
              # 4: 500
              # 5: 600


variables = ['icecover','mlt']
dpis = [100,200,300,400,500,600]


var = variables[flag_var]
dpi = dpis[flag_dpi]

#==============================================================================
# PROC MAIN PLOT
#==============================================================================


sys.path.append(plotDIR)


if flag_var == 0:
    
    from plot_recon_ice_signals import *
    ice_plot(iceDIR,outDIR,flag_svplt)
            
    
elif flag_var == 1:
    
    from plot_recon_mlt_signals import *
    mlt_plot(mltDIR,outDIR,flag_svplt,dpi)
    
            
    

    
    
