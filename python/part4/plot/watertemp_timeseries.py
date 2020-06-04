#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 21:05:53 2018

@author: Luke
"""


#==============================================================================
#SUMMARY
#==============================================================================


# 11 March 2020

#This script loops through annual fieldmeans of each lake model's watertemp data for
#all periods for plotting future projections (x axis time, y axis lake temp)


#==============================================================================
#IMPORT
#==============================================================================


import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


#==============================================================================
#FUNCTIONS
#==============================================================================


# reading netcdfs
def reader(filename):
    ds = xr.open_dataset(filename, decode_times=False)
    da = ds.watertemp
    da = da.where(310>da).interpolate_na(dim='time',method='linear')
    da = da.where(da>260).interpolate_na(dim='time',method='linear')
    dims = list(da.dims)
    discard = []
    if 'lat' in dims:
        discard.append('lat')
    if 'lon' in dims:
        discard.append('lon')
    if 'levlak' in dims:
        discard.append('levlak')
    return da.squeeze(dim=discard,drop=True)

# ensemble math
def ensembler(data):
    concat_dim = np.arange(len(data))
    aligned = xr.concat(data,dim=concat_dim)
    mean = aligned.mean(dim='concat_dim')
    std = aligned.std(dim='concat_dim')
    return [mean,std]

# plotting function for subplots (will copy this into main text section plotting function)
def plotter(datasets,model,ax,time_pre,time_his,time_fut,\
            lw_mean,col_pimean,col_pifill,col_histmean,col_histfill,\
            col_rcp26mean,col_rcp26fill,col_rcp60mean,col_rcp60fill,\
            col_rcp85mean,col_rcp85fill,ub_alpha):
    
    if 'picontrol_1661_1860' in datasets[model]: 
        h = ax.plot(time_pre, datasets[model]['picontrol_1661_1860'][0], lw=lw_mean, color=col_pimean, label='PI control', zorder=1)
        ax.fill_between(time_pre, ((datasets[model]['picontrol_1661_1860'][0]) + datasets[model]['picontrol_1661_1860'][1]),\
                        ((datasets[model]['picontrol_1661_1860'][0]) - datasets[model]['picontrol_1661_1860'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_pifill, zorder=1)
    if 'picontrol_1861_2005' in datasets[model]:
        h = ax.plot(time_his, datasets[model]['picontrol_1861_2005'][0], lw=lw_mean, color=col_pimean, label='_nolegend_', zorder=1)
        ax.fill_between(time_his, ((datasets[model]['picontrol_1861_2005'][0]) + datasets[model]['picontrol_1861_2005'][1]),\
                        ((datasets[model]['picontrol_1861_2005'][0]) - datasets[model]['picontrol_1861_2005'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_pifill, zorder=1)
    if 'picontrol_2006_2099' in datasets[model]:
        h = ax.plot(time_fut, datasets[model]['picontrol_2006_2099'][0], lw=lw_mean, color=col_pimean, label='_nolegend_', zorder=1)
        ax.fill_between(time_fut, ((datasets[model]['picontrol_2006_2099'][0]) + datasets[model]['picontrol_2006_2099'][1]),\
                        ((datasets[model]['picontrol_2006_2099'][0]) - datasets[model]['picontrol_2006_2099'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_pifill, zorder=1)
    if 'historical_1861_2005' in datasets[model]:
        h = ax.plot(time_his, datasets[model]['historical_1861_2005'][0], lw=lw_mean, color=col_histmean, label='Historical', zorder=4)
        ax.fill_between(time_his, ((datasets[model]['historical_1861_2005'][0]) + datasets[model]['historical_1861_2005'][1]),\
                        ((datasets[model]['historical_1861_2005'][0]) - datasets[model]['historical_1861_2005'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_histfill, zorder=4)
    if 'rcp26_2006_2099' in datasets[model]:
        h = ax.plot(time_fut, datasets[model]['rcp26_2006_2099'][0], lw=lw_mean, color=col_rcp26mean, label='RCP 2.6', zorder=3)
        ax.fill_between(time_fut, ((datasets[model]['rcp26_2006_2099'][0]) + datasets[model]['rcp26_2006_2099'][1]),\
                        ((datasets[model]['rcp26_2006_2099'][0]) - datasets[model]['rcp26_2006_2099'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_rcp26fill, zorder=3) 
    if 'rcp60_2006_2099' in datasets[model]:
        h = ax.plot(time_fut, datasets[model]['rcp60_2006_2099'][0], lw=lw_mean, color=col_rcp60mean, label='RCP 6.0', zorder=2)
        ax.fill_between(time_fut, ((datasets[model]['rcp60_2006_2099'][0]) + datasets[model]['rcp60_2006_2099'][1]),\
                        ((datasets[model]['rcp60_2006_2099'][0]) - datasets[model]['rcp60_2006_2099'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_rcp60fill, zorder=2)  
    if 'rcp85_2006_2099' in datasets[model]:
        h = ax.plot(time_fut, datasets[model]['rcp85_2006_2099'][0], lw=lw_mean, color=col_rcp85mean, label='RCP 8.5', zorder=1)
        ax.fill_between(time_fut, ((datasets[model]['rcp85_2006_2099'][0]) + datasets[model]['rcp85_2006_2099'][1]),\
                    ((datasets[model]['rcp85_2006_2099'][0]) - datasets[model]['rcp85_2006_2099'][1]),
                    lw=0.1, alpha=ub_alpha, color=col_rcp85fill, zorder=1) 
    return h,ax

# processing function
def wt_proc(inDIR,flag_main_SI):
    
    # inDIR and flag_man_SI provided by wrapper (main_p2):
        # flag_main_SI = 0; perform calculations for mmm
        # flag_main_SI = 1; perform calculations for SI subplots

    #==============================================================================
    #INITIALIZE DIRECTORIES
    #==============================================================================
    
    # set directory
    os.chdir(inDIR)
    
    #==============================================================================
    #DEFINE SETTINGS
    #==============================================================================
    
    # settings                
    flag_var=0;  # 0: watertemp
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
    models = ['clm45','albm','lake','simstrat-uog','vic-lake']
    
    # list of variables (1st root variables, then processed)
    variables = ['watertemp','lakeicefrac','icethick','icestart','iceend','icedur']
    
    # list of scenarios
    scenarios = ['picontrol','historical','rcp26','rcp60','rcp85']
    
    # list of forcings
    forcings = ['gfdl-esm2m','hadgem2-es','ipsl-cm5a-lr','miroc5']
    
    # list of products
    products = ['fldmean','sig','scaled_sig','eval']
    
    # list of end periods
    endperiods = ['1661_1860','1861_2005','2006_2099']
    
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
        
    available_models = []
    for mod in models:
        count = 0
        for file in files:
            if mod in file:
                count = count+1
                if count == 1:
                    available_models.append(mod)
        
    # proc for mmm    
    if flag_main_SI == 0:
        
        # first get pi_refs per mod per gcm
        pi_refs = {}
        for mod in available_models:
            pi_refs[mod] = {}
            for gcm in forcings:
                temp_list = []
                for file in files:
                    if mod in file and gcm in file and 'picontrol' in file:
                        temp_list.append(reader(file))
                pi_refs[mod][gcm] = xr.concat(temp_list,dim='time').mean(dim='time').values
        
        # read in each file as an anomaly
        arrays = {}
        for mod in available_models:
            arrays[mod] = {}
            for gcm in forcings:
                arrays[mod][gcm] = {}
                for scen in scenarios:
                    for endper in endperiods:
                        for file in files:
                            if mod in file\
                            and gcm in file\
                            and scen in file\
                            and endper in file: # take anomaly with pi_ref
                                arrays[mod][gcm][scen+'_'+endper] = reader(file) - pi_refs[mod][gcm]
        
        # gather ensembles per scenario + endperiod
        datasets = {}
        keys = list(arrays['clm45']['gfdl-esm2m'].keys()) # scenario + endperiod keys
        for key in keys:
            datasets[key] = []
            for mod in available_models:
                for gcm in forcings:
                    if key in arrays[mod][gcm]:
                        datasets[key].append(arrays[mod][gcm][key])
                    
        # calculate mmm
        mmm = {}
        for key in keys:
            mmm[key] = ensembler(datasets[key])
# =============================================================================
#             if '2006_2099' in key:
#                 new_mean = mmm[key][0].shift(time=-1)
#                 new_std = mmm[key][1].shift(time=-1)
#                 mmm[key] = [new_mean,new_std]
# =============================================================================

        return mmm,keys

    
    # proc for subplots
    if flag_main_SI == 1:
        
        # arrange + read files
        filesets = {}
        datasets = {}
        for mod in models:
            filesets[mod] = {}
            datasets[mod] = {}
            for scen in scenarios:
                for endper in endperiods:
                    filesets[mod][scen+'_'+endper] = []
                    for file in files:
                        if mod in file and scen in file and endper in file:
                            filesets[mod][scen+'_'+endper].append(reader(file))
                    if not filesets[mod][scen+'_'+endper]:  #check if dict is empty -> remove
                        del filesets[mod][scen+'_'+endper]
        
        # ensemble statistics for plot
        for mod in filesets:
            for scen in scenarios:
                for endper in endperiods:
                    if scen+'_'+endper in filesets[mod]:
                        datasets[mod][scen+'_'+endper] = ensembler(filesets[mod][scen+'_'+endper])
                        
        # picontrol ref mean calculation              
        pi_refs = {}
        for mod in models: 
            if datasets[mod]:
                temp_list = []
                if 'picontrol_1661_1860' in datasets[mod]:
                    temp_list.append(datasets[mod]['picontrol_1661_1860'][0])
                if 'picontrol_1861_2005' in datasets[mod]:
                    temp_list.append(datasets[mod]['picontrol_1861_2005'][0])            
                if 'picontrol_2006_2099' in datasets[mod]:
                    temp_list.append(datasets[mod]['picontrol_2006_2099'][0])   
                pi_refs[mod] = xr.concat(temp_list,dim='time').mean(dim='time').values
            if not datasets[mod]:
                del datasets[mod]
                
        # picontrol ref mean subtraction
        anoms = {}
        for mod in datasets:
            anoms[mod] = {}
            for key in datasets[mod]:
                anoms[mod][key] = []
                anoms[mod][key].append(datasets[mod][key][0] - pi_refs[mod]) # insert mean anomaly
                anoms[mod][key].append(datasets[mod][key][1]) # insert std
                
        return anoms
    
        
# plotting function for subplots
def wt_plot(datasets,outDIR,var,prod,flag_svplt,dpi):
    
    # list of model titles
    model_titles = ['CLM4.5', 'ALBM', 'SIMSTRAT-UoG', 'VIC-LAKE']
    
    #==============================================================================
    #DEFINE PLOT SETTINGS
    #==============================================================================
    
    # list of figure panel ids
    letters = ['a', 'b', 'c',\
               'd', 'e', 'f',\
               'g', 'h', 'i',\
               'j', 'k', 'l']
    
    #========== FONTS ==========#
    
    title_font = 9
    tick_font = 8
    axis_font = 9
    legend_font = 8
    
    #========== LINE THICKNESS ==========#
    
    # mean line thickness
    lw_mean = 1
    
    #========== PLOT COLORS ==========#
    
    col_pimean = 'blue'         # picontrol mean color
    col_pifill = '#a6bddb'      # picontrol fill color
    col_histmean = '0.3'       # historical mean color
    col_histfill = '0.75'       # historical fill color
    col_rcp26mean = 'darkgreen'       # rcp26 mean color
    col_rcp26fill = '#adebad'     # rcp26 fill color
    col_rcp60mean = 'darkgoldenrod'   # rcp60 mean color
    col_rcp60fill = '#ffec80'     # rcp60 fill color
    col_rcp85mean = 'darkred'       # rcp85 mean color
    col_rcp85fill = '#F08080'     # rcp85 fill color
    ub_alpha = 0.5
    
    #========== AXII ==========#
    
    # figsize = (x,y)
    x = 5
    y = 8
    
    # subplots_adjust
    hspace = 0.3
    top = 0.9
    
    ymin = -1   # ymin
    ymax = 5    # ymax
    xmin = 1900 # xmin
    xmax = 2100 # xmax
    
    # x ticks/labels
    xticks = np.arange(1900,2125,25)
    xtick_labels = [None,None,1950,None,2000,None,2050,None,2100]
    
    # x axis label
    xlabel = 'Years'
    xlabel_xpos = 0.5
    xlabel_ypos = 0.05
    
    # y axis label
    ylabel = 'Lake temperature anomaly (°C)'
    ylabel_xpos = 0.025
    ylabel_ypos = 0.535
    
    # xaxis tick label sharing
    axis_share = False
    
    #========== LEGEND ==========#
    
    # labels
    lab_pimean = 'Pre-industrial control'
    lab_histmean = 'Historical'
    lab_26mean = 'RCP 2.6'
    lab_60mean = 'RCP 6.0'
    lab_85mean = 'RCP 8.5'
    
    # bbox
    x0 = -0.1
    y0 = 1.15
    xlen = 1.15
    ylen = 0.9
    
    # space between entries
    legend_entrypad = 0.5
    
    # length per entry
    legend_entrylen = 0.75
    
    #==============================================================================
    #TIME
    #==============================================================================
    
    # define time variables
    time_pre = np.arange(1661,1861,1)
    time_his = np.arange(1861,2006,1)
    time_fut = np.arange(2006,2100,1)
    
    #==============================================================================
    #PLOTTING
    #==============================================================================
      
    # initiate plots
    f, axes = plt.subplots(len(datasets),1,figsize=(x,y),sharex=axis_share,)
    
    # load data 
    
    count = 0
    
    for model,ax in zip(datasets,axes.flatten()):
        
        count += 1
        
        h,ax = plotter(datasets,model,ax,time_pre,time_his,time_fut,\
            lw_mean,col_pimean,col_pifill,col_histmean,col_histfill,\
            col_rcp26mean,col_rcp26fill,col_rcp60mean,col_rcp60fill,\
            col_rcp85mean,col_rcp85fill,ub_alpha)
            
    # figure adjustments
        ax.set_title(model_titles[count-1],loc='right',fontsize=title_font)
        ax.set_title(letters[count-1],loc='left',fontsize=title_font,fontweight='bold')
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels(xtick_labels)
        ax.tick_params(labelsize=tick_font,axis="x",direction="in", left="off",labelleft="on")
        ax.tick_params(labelsize=tick_font,axis="y",direction="in")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.grid(color='0.8', linestyle='dashed', linewidth=0.5)
        ax.xaxis.grid(color='0.8', linestyle='dashed', linewidth=0.5)
        ax.set_axisbelow(True)
    
    # legend
    legendcols = [col_pimean,col_histmean,col_rcp26mean,col_rcp60mean,col_rcp85mean]
    handles = [Line2D([0],[0],linestyle='-',lw=2,color=legendcols[0]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[1]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[2]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[3]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[4])]
    labels= [lab_pimean,lab_histmean,lab_26mean,lab_60mean,lab_85mean]
    axes[0].legend(handles, labels, bbox_to_anchor=(x0, y0, xlen, ylen), loc=3,   #bbox: (x, y, width, height)
               ncol=6,fontsize=legend_font, mode="expand", borderaxespad=0.,\
               frameon=False, columnspacing=0.05, handlelength=legend_entrylen, handletextpad=legend_entrypad)
    
    # labels
    f.text(xlabel_xpos, xlabel_ypos, xlabel, ha='center', fontsize=axis_font)
    f.text(ylabel_xpos, ylabel_ypos, ylabel, va='center', rotation='vertical', fontsize=axis_font)
    
    f.subplots_adjust(hspace=hspace,top=top)
    
    plt.show()
    
    # save figure
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        f.savefig(outDIR+'/si_f23.png',bbox_inches='tight',dpi=dpi)
    