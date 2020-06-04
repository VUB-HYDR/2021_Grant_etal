#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:50:59 2020

@author: Luke
"""

#==============================================================================
# SUMMARY
#==============================================================================


# 18 May 2020

# plots detection and attribution analysis


#==============================================================================
# IMPORT
#==============================================================================


import os
import sys
import xarray as xr
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats
import seaborn as sns
from scipy.stats import norm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#==============================================================================
# FUNCTION
#==============================================================================

def scale_take(array): #must take diff between b and sup/inf, store in separate lists
    b = array[1]
    b_inf = b - array[0]
    b_sup = array[2] - b
    p = array[3]
    return b,b_inf,b_sup,p

def plot_p2(outDIR,flag_svplt,endvariables,\
             hist_mmm,pi_mmm,era5_obs,\
             samples,mean,n,std,\
             histmmm_obs_pcc,histmmm_obs_spcc,\
             pi_histmmm_pcc,pi_histmmm_spcc,\
             cc_99,cc_95,cc_90,\
             var_fin):
    
    #==============================================================================
    # INITIALIZE
    #==============================================================================
    
    # figure size
    f = plt.figure(figsize=(10,14))
    
    #========== TEMPORAL ==========#
    
    # opt fing rect, rect=[left, bottom, right, top]
    t_left = 0.475
    t_bottom = 0.025
    t_right = 1
    t_top = 1.0
    t_rect = [t_left, t_bottom, t_right, t_top]
    
    # temporal - gs1; ax1,ax2,ax3, ax4
    gs1 = gridspec.GridSpec(4,1)
    ax1 = f.add_subplot(gs1[0])    
    ax2 = f.add_subplot(gs1[1])    
    ax3 = f.add_subplot(gs1[2])    
    ax4 = f.add_subplot(gs1[3])    
    gs1.tight_layout(figure=f, rect=t_rect, h_pad=5)
    
    temp_axes = [ax1,ax2,ax3,ax4]

    #========== CORRELATION ==========#     
    
    gs5 = gridspec.GridSpec(4,1)
    
    # corr detection rect, rect=[left, bottom, right, top]
    c_left = 0
    c_bottom = 0.025
    c_right = 0.425
    c_top = 1.0
    c_rect = [c_left, c_bottom, c_right, c_top]
    ax13 = f.add_subplot(gs5[0])
    ax14 = f.add_subplot(gs5[1])
    ax15 = f.add_subplot(gs5[2])
    ax16 = f.add_subplot(gs5[3])
    corr_axes = [ax13,ax14,ax15,ax16]
    gs5.tight_layout(figure=f, rect=c_rect, h_pad=5)
    
    #==============================================================================
    # GENERAL TIMESERIES + SCALING SETTINGS
    #==============================================================================
    
    # list of figure panel ids
    letters = ['a', 'b', 'c',\
               'd', 'e', 'f',\
               'g', 'h', 'i',\
               'j', 'k', 'l']
    
    #========== LINE THICKNESS ==========#
    
    # mean line thickness
    lw_mean = 1.0
    
    # era5-land line thickness
    lw_era5 = 1.0
    
    #========== PLOT COLORS ==========#
    
    col_pimean = 'dodgerblue'         # picontrol mean color
    col_pifill = '#a6bddb'      # picontrol fill color
    col_PICmean = 'mediumblue'    # PIC block mean color
    col_histmean = '0.3'       # historical mean color
    col_histfill = '0.75'       # historical fill color
    col_rcp26mean = 'darkgreen'       # rcp26 mean color
    col_rcp26fill = '#adebad'     # rcp26 fill color
    col_rcp60mean = 'darkgoldenrod'   # rcp60 mean color
    col_rcp60fill = '#ffec80'     # rcp60 fill color
    col_rcp85mean = 'darkred'       # rcp85 mean color
    col_rcp85fill = '#F08080'     # rcp85 fill color
    col_ALLmean = 'red'     # ALL block mean color   
    ub_alpha = 0.5
    col_era5 = 'k'
    col_OBSmean = 'k'     # OBS block mean color
    
    #========== AXII ==========#
    
    # ymin
    ymin_ice = -10   # ymin ice vars
    ymax_ice = 10   # ymax ice vars
    ymin_wat = -0.25   # ymin watertemp 
    ymax_wat = 1   # ymax watertemp
    xmin = 1985   # xmin
    xmax = 2018   # xmax
    
    # y ticks; ice
    yticks_ice = np.arange(ymin_ice,ymax_ice+2.5,2.5)
    ytick_labels_ice = [-10,None,-5,None,0,None,5,None,10]
    
    # y ticks; watertemp
    yticks_wat = np.arange(ymin_wat,ymax_wat+0.25,0.25)
    ytick_labels_wat = [None,0,None,0.5,None,1.0]
    
    # x ticks; timeseries
    xticks_ts = np.arange(1990,2020,5)
    xtick_labels_ts = [1990,None,2000,None,2010,None]
    
    # x ticks temporal OF insets
    xticks_OF = [0.5]
    xtick_labels_OF = ['EXT']
    
    # y ticks temporal OF insets
    yticks_OF = np.arange(-0.5,2.5,0.5)
    ytick_labels_OF = [None, '0', None, '1', None, '2']
        
    #========== FONTS ==========#
    
    title_font = 14
    tick_font = 10
    axis_font = 11
    legend_font = 11
    inset_font = 9
    
    #==============================================================================
    # CORRELATION
    #==============================================================================
    
    # main plot x axis label
    xlabel = 'Spearman (rank) correlation coefficient'
    xlabel_xpos = 0.225
    xlabel_ypos = 0.0
    
    ylabel = 'Density [-]'
    
    #========== LEGEND ==========#
    
    # bbox
    le_x0 = 0.365
    le_y0 = 0.48
    le_xlen = 0.15
    le_ylen = 0.25
    
    # space between entries
    legend_entrypad = 0.5
    
    # length per entry
    legend_entrylen = 0.75
    
    # space between entries
    legend_spacing = 1.5
    
    yticks_corr = np.arange(0,3,0.5)
    ytick_labels_corr = [0,None,1.0,None,2.0,None,]
    
    xticks_corr = np.arange(-0.8,1.0,0.2)
    xtick_labels_corr = [-0.8,None,-0.4,None,0,None,0.4,None,0.8]
    
    count = 0
    
    for endvar,ax in zip(endvariables,corr_axes):
        
        count += 1
        
        sns.distplot(samples[endvar], bins=20, hist_kws={"color":"0.4"},fit=norm, fit_kws={"color":"0.4"},\
                     norm_hist=True, kde=False, label='PIC, EXT', ax=ax)
        ax.axvline(x=histmmm_obs_spcc[endvar], color='red', linewidth=2, label='EXT, ERA5-land')
        ax.vlines(x=cc_99[endvar], ymin=0, ymax=0.125, colors='blue', linewidth=2.5, linestyle='-', label='99%', zorder=0)
        ax.vlines(x=cc_95[endvar], ymin=0, ymax=0.125, colors='blue', linewidth=1, linestyle='-', label='95%', zorder=0)
        ax.tick_params(labelsize=tick_font)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(letters[count-1],loc='left',fontsize=title_font,fontweight='bold')
        ax.yaxis.set_ticks(yticks_corr)
        ax.yaxis.set_ticklabels(ytick_labels_corr)
        ax.set_ylabel(ylabel)
        ax.set_ylim(yticks_corr[0],yticks_corr[-1])
        ax.xaxis.set_ticks(xticks_corr)
        ax.xaxis.set_ticklabels(xtick_labels_corr)
        
    ax13.legend(frameon=False,bbox_to_anchor=(le_x0, le_y0, le_xlen, le_ylen),\
                fontsize=legend_font,labelspacing=legend_spacing)
    
    f.text(xlabel_xpos, xlabel_ypos, xlabel, ha='center', fontsize=axis_font)
    
    #==============================================================================
    # OF
    #==============================================================================
    
    # main plot x axis label
    xlabel = 'Years'
    xlabel_xpos = 0.735
    xlabel_ypos = 0.0
    
    #========== LEGEND ==========#
    
    # bbox
    le_x0 = 0.75
    le_y0 = 0.975
    le_xlen = 0.2
    le_ylen = 0.1
    
    # space between marker and label
    legend_entrypad = 0.5
    
    # length per entry
    legend_entrylen = 0.75
    
    # space between entries
    legend_spacing = 2.2
    
    time_ice = np.arange(1981,2018)
    time_wt = np.arange(1981,2018)
    
    for endvar,ax in zip(endvariables,temp_axes):
        
        if endvar == 'watertemp':
            time = time_wt
        else:
            time = time_ice
        
        count += 1
        
        # OF data prep
        b85,b_inf85,b_sup85,p85 = scale_take(var_fin[endvar])
        
            
        infers = [b_inf85]
        supers = [b_sup85]
        err = np.stack([infers,supers],axis=0)
        x = [0.5]
        y = [b85]

        # timeseries
        ax.plot(time, pi_mmm[endvar][0], lw=lw_mean, color=col_pimean, label='PIC', zorder=1)
        ax.fill_between(time, (pi_mmm[endvar][0] + pi_mmm[endvar][1]),\
                        ((pi_mmm[endvar][0]) - pi_mmm[endvar][1]),\
                        lw=0.1, color=col_pifill, zorder=1)
        ax.plot(time, hist_mmm[endvar][0], lw=lw_mean, color=col_rcp85mean, label='EXT', zorder=1)
        ax.fill_between(time, (hist_mmm[endvar][0] + hist_mmm[endvar][1]),\
                    (hist_mmm[endvar][0] - hist_mmm[endvar][1]),
                    lw=0.1, color=col_rcp85fill, zorder=1 ,alpha=ub_alpha)
        ax.plot(time, era5_obs[endvar], lw=lw_era5, color=col_era5, label='OBS', zorder=4)
        
        if count == 5:
            ax_ins = inset_axes(ax, width="20%", height="35%", loc=2, borderpad=3)
        elif count == 6:
            ax_ins = inset_axes(ax, width="20%", height="35%", loc=2, borderpad=3)
        elif count == 7:
            ax_ins = inset_axes(ax, width="20%", height="35%", loc=3, borderpad=3)
        elif count == 8:
            ax_ins = inset_axes(ax, width="20%", height="35%", loc=3, borderpad=3)
            
        ax_ins.errorbar(x=x,y=y,yerr=err,
                        fmt='o',
                        markersize=3,
                        ecolor=col_rcp85mean,
                        markerfacecolor=col_rcp85mean,
                        mec=col_rcp85mean,
                        capsize=5,
                        elinewidth=2,
                        markeredgewidth=1)
        
        ax_ins.set_ylim(-0.5,2)
        ax_ins.set_xlim(0,1)
        
        ax_ins.hlines(y=1,xmin=0,xmax=3,colors='k',linestyle='dashed',linewidth=1)
        ax_ins.hlines(y=0,xmin=0,xmax=3,colors='k',linestyle='solid',linewidth=0.25)
        
        
        ax_ins.xaxis.set_ticks(xticks_OF)
        ax_ins.xaxis.set_ticklabels(xtick_labels_OF,fontsize=inset_font)
        
        ax_ins.yaxis.set_ticks(yticks_OF)
        ax_ins.yaxis.set_ticklabels(ytick_labels_OF,fontsize=inset_font)
    
        # settings for timeseries plots
        if 'ice' in endvar:
            if count == 2:
                ax.set_ylim(ymin_ice+5,ymax_ice)
                ax.yaxis.set_ticks(yticks_ice[2:])
                ax.yaxis.set_ticklabels(ytick_labels_ice[2:])
            elif count == 3:
                ax.set_ylim(ymin_ice,ymax_ice-5)
                ax.yaxis.set_ticks(yticks_ice[:-2])
                ax.yaxis.set_ticklabels(ytick_labels_ice[:-2])
            elif count ==4:
                ax.set_ylim(ymin_ice,ymax_ice-5)
                ax.yaxis.set_ticks(yticks_ice[:-2])
                ax.yaxis.set_ticklabels(ytick_labels_ice[:-2])
        elif 'water' in endvar:
            ax.set_ylim(ymin_wat,ymax_wat)
            ax.yaxis.set_ticks(yticks_wat)
            ax.yaxis.set_ticklabels(ytick_labels_wat)
            
        ax.set_xlim(xmin,xmax)   
        ax.xaxis.set_ticks(xticks_ts)
        ax.xaxis.set_ticklabels(xtick_labels_ts)
        ax.set_title(letters[count-1],loc='left',fontsize=title_font,fontweight='bold')
        ax.tick_params(labelsize=tick_font,axis="x",direction="in",labelleft="on",left=True)
        ax.tick_params(labelsize=tick_font,axis="y",direction="in",labelbottom="on",bottom=True,top=False,)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.grid(color='0.8', linestyle='dashed', linewidth=0.5)
        ax.xaxis.grid(color='0.8', linestyle='dashed', linewidth=0.5)
        ax.set_axisbelow(True)
        
        # labels
        if endvar == 'icedur':
            ylabel = 'Ice duration anomaly (days)'
        elif endvar == 'icestart':
            ylabel = 'Ice onset anomaly (days)'
        elif endvar == 'iceend':
            ylabel = 'Ice break-up anomaly (days)'
        elif endvar == 'watertemp':
            ylabel = 'Water temperature anomaly (°C)'
        
        ax.set_ylabel(ylabel, va='center', rotation='vertical', fontsize=axis_font, labelpad=10)

        
        labels = ['EXT', 'OBS', 'PIC']
        handles = [Rectangle((0,0),1,1,color=col_rcp85fill),\
           Line2D([0],[0],linestyle='-',lw=2,color=col_era5),\
           Rectangle((0,0),1,1,color=col_pifill)]
            
    f.legend(handles, labels, bbox_to_anchor=(le_x0, le_y0, le_xlen, le_ylen), loc=3,
                   ncol=3, mode="expand", borderaxespad=0.,\
                   frameon=False, columnspacing=0.1, handlelength=legend_entrylen,\
                   handletextpad=legend_entrypad,\
                   fontsize=legend_font,labelspacing=legend_spacing)
           
    f.text(xlabel_xpos, xlabel_ypos, xlabel, ha='center', fontsize=axis_font)
        
    plt.show()
    
    # save figure
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        f.savefig(outDIR+'/f2.png',bbox_inches='tight',dpi=500)
