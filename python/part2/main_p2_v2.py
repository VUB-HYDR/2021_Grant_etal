#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 08:15:24 2020

@author: Luke
"""

#==============================================================================
# SUMMARY
#==============================================================================


# main script for calling detection and attribution analysis (f2)


#==============================================================================
# IMPORT
#==============================================================================


import sys
import os
import numpy as np
import pickle as pk


#==============================================================================
# PATH
#==============================================================================


# change current working directory to pathway of this file
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# object for current dir
curDIR = os.getcwd()

# data input directory
inDIR = os.path.join(curDIR, 'data/all_v2/')
piDIR = os.path.join(curDIR, 'data/pichunks_v2/')
obsDIR = os.path.join(curDIR, 'data/obs_v2/')

# script directory
plotDIR = os.path.join(curDIR, 'plot/')

# figure output directory
outDIR = os.path.join(curDIR, 'figures/final/')


#==============================================================================
# OPTIONS
#==============================================================================


# adjust these settings for either testing OF output or producing main text fig

flag_svplt=0;     # 0: do not save plot
                  # 1: save plot in outDIR

             
flag_main_OF=0;   # 0: process main text fig
                  # 1: process OF output

# choices for individual OF outputs    
           
flag_fp_var=0;    # 0: watertemp
                  # 1: icestart
                  # 2: iceend
                  # 3: icedur

flag_fp_scen=2;   # 0: rcp85
                  # 1: rcp60
                  # 2: ctl

flag_reg=1;       # 0: OLS
                  # 1: TLS

flag_constest=3;  # 0: OLS_AT99 
                  # 1: OLS_Corr
                  # 2: AS03 (TLS only)
                  # 3: MC (TLS only)

# confidence internval calculation in case that TLS regression chosen 

flag_ci_tls=1;    # 0: AS03
                  # 1: ODP

# len of rolling mean window

flag_rolllen=3;   # 0: no roll
                  # 1: 3-year
                  # 2: 4-year
                  # 3: 5-year

flag_block=1;     # 0: no block means in ROF (if not block means, rolling window taken)
                  # 1: yes block means in ROF (take rolllen as block size)

# shuffle option for reading in t-series (was necessary with temporal detection)

flag_shuffle=0    # 0: No shuffling/limitation of num pi series in tseries
                  # 1: shuffle (temporal approach)

# bootstrap repetitions

flag_bs=4         # 0: No bootstrapping of covariance matrix build
                  # 1: 50 (e.g. 50 reps of ROF, each with shuffled pichunks for Cov_matrices)
                  # 2: 100
                  # 3: 500
                  # 4: 1000

# confidence intervals on scaling factors

flag_ci_bnds=2    # 0: 10-90
                  # 1: 5-95
                  # 2: 2.5-97.5
                  # 3: 0.5-99.5


endvariables = ['watertemp','icestart','iceend','icedur']
scenarios = ['rcp85','rcp60','ctl']
regressions = ['OLS','TLS']
consistency_tests = ['OLS_AT99','OLS_Corr','AS03','MC']
tls_cis = ['AS03','ODP']
shuffle_opts = ['no', 'yes']
blockmeans = ['no', 'yes']
rolls = [0,3,4,5]
bootstrap_reps = [0,50,100,500,1000]
confidence_intervals = [0.8,0.9,0.95,0.99]


fp_var = endvariables[flag_fp_var]
fp_scen = scenarios[flag_fp_scen]
reg = regressions[flag_reg]
cons_test = consistency_tests[flag_constest]
formule_ic_tls = tls_cis[flag_ci_tls]
window = rolls[flag_rolllen]
block = blockmeans[flag_block]
shuff = shuffle_opts[flag_shuffle]
bs_reps = bootstrap_reps[flag_bs]
ci_bnds = confidence_intervals[flag_ci_bnds]


#==============================================================================
# PROC MAIN PLOT
#==============================================================================


sys.path.append(plotDIR)


if flag_main_OF == 0:
    
    os.chdir(obsDIR)
    
    if not os.path.isfile('pk_data_p2.pkl'):
        
        print('no pickled data')
    
        from temporal_v2 import *
        hist_mmm,pi_mmm,era5_obs = temporal_proc(inDIR,piDIR,obsDIR,shuff,window)
                
        from da_flex_v2 import data_loader   
        import PyDnA_v2 as pda
        var_sfs = {}
        scen = 'rcp85'
        bhi = {}
        b = {}
        blow = {}
        pval = {}
        var_fin = {}
        var_xruns = {}
        var_ctlruns = {}
        
        for var in endvariables:
            bhi[var] = []
            b[var] = []
            blow[var] = []
            pval[var] = []
            if bs_reps == 0: # for no bs, run ROF once
                bs_reps += 1
            for i in np.arange(0,bs_reps):
                obs,fp,ctl,nx = data_loader(inDIR,piDIR,obsDIR,var,scen,window,block)
                y = obs
                X = fp
                nb_runs_x= nx
                trunc = 0 # constant for our global means
                var_sfs[var],var_ctlruns[var] = pda.da_main(y,X,ctl,nb_runs_x,reg,cons_test,formule_ic_tls,trunc,ci_bnds)
                bhi[var].append(var_sfs[var][0])
                b[var].append(var_sfs[var][1])
                blow[var].append(var_sfs[var][2])
                pval[var].append(var_sfs[var][3])
                print(str(var)+str(i))
            bhi_med = np.median(bhi[var])
            b_med = np.median(b[var])
            blow_med = np.median(blow[var])
            pval_med = np.median(pval[var])
            var_fin[var] = [bhi_med,b_med,blow_med,pval_med]
            var_xruns[var] = nb_runs_x[0]
        print('Number of ens members for EXT for watertemp, icestart, iceend and icedur are: '\
              + str(var_xruns['watertemp']) + ', '\
              + str(var_xruns['icestart']) + ', '\
              + str(var_xruns['iceend']) + ', '\
              + str(var_xruns['icedur']) + ', ')
    
        print('Number of PIC chunks for watertemp, icestart, iceend and icedur are: '\
              + str(var_ctlruns['watertemp']) + ', '\
              + str(var_ctlruns['icestart']) + ', '\
              + str(var_ctlruns['iceend']) + ', '\
              + str(var_ctlruns['icedur']) + ', ')
    
                
                
        from correlation_v2 import *
        samples,mean,n,std,\
                   histmmm_obs_pcc,histmmm_obs_spcc,\
                   pi_histmmm_pcc,pi_histmmm_spcc,\
                   cc_99,cc_95,cc_90 = correlation_proc(inDIR,piDIR,obsDIR)
                   
        pickle_data = {}
        
        pickle_data['hist_mmm'] = hist_mmm
        pickle_data['pi_mmm'] = pi_mmm
        pickle_data['era5_obs'] = era5_obs
        
        pickle_data['samples'] = samples
        pickle_data['mean'] = mean
        pickle_data['n'] = n
        pickle_data['std'] = std
        
        pickle_data['histmmm_obs_pcc'] = histmmm_obs_pcc
        pickle_data['histmmm_obs_spcc'] = histmmm_obs_spcc
        
        pickle_data['pi_histmmm_pcc'] = pi_histmmm_pcc
        pickle_data['pi_histmmm_spcc'] = pi_histmmm_spcc
        
        pickle_data['cc_99'] = cc_99
        pickle_data['cc_95'] = cc_95
        pickle_data['cc_90'] = cc_90
        
        pickle_data['var_fin'] = var_fin
        
        output = open('pk_data_p2.pkl','wb')
        pk.dump(pickle_data,output)
        output.close()
        
    elif os.path.isfile('pk_data_p2.pkl'):
        
        print('use pickled data')
        
        pkl_file = open('pk_data_p2.pkl','rb')
        pickle_data = pk.load(pkl_file)
        pkl_file.close()
        
        hist_mmm = pickle_data['hist_mmm']
        pi_mmm = pickle_data['pi_mmm'] 
        era5_obs = pickle_data['era5_obs'] 
        
        samples = pickle_data['samples']
        mean = pickle_data['mean']
        n = pickle_data['n']
        std = pickle_data['std']
        
        histmmm_obs_pcc = pickle_data['histmmm_obs_pcc']
        histmmm_obs_spcc = pickle_data['histmmm_obs_spcc']
        
        pi_histmmm_pcc = pickle_data['pi_histmmm_pcc']
        pi_histmmm_spcc = pickle_data['pi_histmmm_spcc']
        
        cc_99 = pickle_data['cc_99']
        cc_95 = pickle_data['cc_95']
        cc_90 = pickle_data['cc_90']
        
        var_fin = pickle_data['var_fin']
    
        
    from plot_p2_v2 import *
    plot_p2(outDIR,flag_svplt,endvariables,\
                 hist_mmm,pi_mmm,era5_obs,\
                 samples,mean,n,std,\
                 histmmm_obs_pcc,histmmm_obs_spcc,\
                 pi_histmmm_pcc,pi_histmmm_spcc,\
                 cc_99,cc_95,cc_90,\
                 var_fin)
        
        
#==============================================================================
# PROC OF INDIVIDUALLY
#==============================================================================


elif flag_main_OF == 1:
    
    from da_flex_v2 import data_loader    
    obs,fp,ctl,nx = data_loader(inDIR,piDIR,obsDIR,fp_var,fp_scen,window)
    
    
    import PyDnA_v2 as pda
    y = obs
    X = fp
    nb_runs_x= nx
    trunc = 0 # constant for our global means
    
    beta = pda.da_main(y,X,ctl,nb_runs_x,reg,cons_test,formule_ic_tls,trunc)
    
    
