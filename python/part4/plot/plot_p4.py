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

# plots global mean anomalies


#==============================================================================
# IMPORT
#==============================================================================


import xarray as xr
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


#==============================================================================
# FUNCTION
#==============================================================================


def ts_plotter(mmm,ax,time_pre,time_his,time_fut,\
            lw_mean,col_pimean,col_pifill,col_histmean,col_histfill,\
            col_rcp26mean,col_rcp26fill,col_rcp60mean,col_rcp60fill,\
            col_rcp85mean,col_rcp85fill,ub_alpha):
    
    if 'picontrol_1661_1860' in mmm: 
        h = ax.plot(time_pre, mmm['picontrol_1661_1860'][0], lw=lw_mean, color=col_pimean, label='PI control', zorder=1)
        ax.fill_between(time_pre, ((mmm['picontrol_1661_1860'][0]) + mmm['picontrol_1661_1860'][1]),\
                        ((mmm['picontrol_1661_1860'][0]) - mmm['picontrol_1661_1860'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_pifill, zorder=1)
    if 'picontrol_1861_2005' in mmm:
        h = ax.plot(time_his, mmm['picontrol_1861_2005'][0], lw=lw_mean, color=col_pimean, label='_nolegend_', zorder=1)
        ax.fill_between(time_his, ((mmm['picontrol_1861_2005'][0]) + mmm['picontrol_1861_2005'][1]),\
                        ((mmm['picontrol_1861_2005'][0]) - mmm['picontrol_1861_2005'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_pifill, zorder=1)
    if 'picontrol_2006_2099' in mmm:
        h = ax.plot(time_fut, mmm['picontrol_2006_2099'][0], lw=lw_mean, color=col_pimean, label='_nolegend_', zorder=1)
        ax.fill_between(time_fut, ((mmm['picontrol_2006_2099'][0]) + mmm['picontrol_2006_2099'][1]),\
                        ((mmm['picontrol_2006_2099'][0]) - mmm['picontrol_2006_2099'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_pifill, zorder=1)
    if 'historical_1861_2005' in mmm:
        h = ax.plot(time_his, mmm['historical_1861_2005'][0], lw=lw_mean, color=col_histmean, label='Historical', zorder=4)
        ax.fill_between(time_his, ((mmm['historical_1861_2005'][0]) + mmm['historical_1861_2005'][1]),\
                        ((mmm['historical_1861_2005'][0]) - mmm['historical_1861_2005'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_histfill, zorder=4)
    if 'rcp26_2006_2099' in mmm:
        h = ax.plot(time_fut, mmm['rcp26_2006_2099'][0], lw=lw_mean, color=col_rcp26mean, label='RCP 2.6', zorder=3)
        ax.fill_between(time_fut, ((mmm['rcp26_2006_2099'][0]) + mmm['rcp26_2006_2099'][1]),\
                        ((mmm['rcp26_2006_2099'][0]) - mmm['rcp26_2006_2099'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_rcp26fill, zorder=3) 
    if 'rcp60_2006_2099' in mmm:
        h = ax.plot(time_fut, mmm['rcp60_2006_2099'][0], lw=lw_mean, color=col_rcp60mean, label='RCP 6.0', zorder=2)
        ax.fill_between(time_fut, ((mmm['rcp60_2006_2099'][0]) + mmm['rcp60_2006_2099'][1]),\
                        ((mmm['rcp60_2006_2099'][0]) - mmm['rcp60_2006_2099'][1]),\
                        lw=0.1, alpha=ub_alpha, color=col_rcp60fill, zorder=2)  
    if 'rcp85_2006_2099' in mmm:
        h = ax.plot(time_fut, mmm['rcp85_2006_2099'][0], lw=lw_mean, color=col_rcp85mean, label='RCP 8.5', zorder=1)
        ax.fill_between(time_fut, ((mmm['rcp85_2006_2099'][0]) + mmm['rcp85_2006_2099'][1]),\
                    ((mmm['rcp85_2006_2099'][0]) - mmm['rcp85_2006_2099'][1]),
                    lw=0.1, alpha=ub_alpha, color=col_rcp85fill, zorder=1) 
    return h,ax


def plot_p4(wt_mmm,it_mmm,ii_mmm,\
            wt_air_ens,wt_ens,wt_iy,\
            it_air_ens,it_ens,it_iy,\
            ii_air_ens,ii_ens,ii_iy,\
            outDIR,flag_svplt):
    
    #==============================================================================
    # INITIALIZE
    #==============================================================================
    
    # final size will be (10,12); temporarily change to 5,6
    f = plt.figure(figsize=(15,18))
    
    gs1 = gridspec.GridSpec(3,1)

    ax1 = f.add_subplot(gs1[0])    
    ax2 = f.add_subplot(gs1[1])
    ax3 = f.add_subplot(gs1[2]) 
    
    # rect=[left, bottom, right, top]
    gs1.tight_layout(figure=f, rect=[0, 0.1, 0.5, 0.9], h_pad=5)
    
    gs2 = gridspec.GridSpec(3,1)
    
    ax4 = f.add_subplot(gs2[0])
    ax5 = f.add_subplot(gs2[1])
    ax6 = f.add_subplot(gs2[2])
    
    
    gs2.tight_layout(figure=f, rect=[0.55, 0.1, 1, 0.9], h_pad=5)
    
    
    #==============================================================================
    # GENERAL TIMESERIES + SCALING SETTINGS
    #==============================================================================
    
    # list of figure panel ids
    letters = ['a', 'b', 'c',\
               'd', 'e', 'f',\
               'g', 'h', 'i',\
               'j', 'k', 'l']
    
    lettercount = 0
    
    #========== FONTS ==========#
    
    title_font = 14
    tick_font = 12
    axis_font = 14
    legend_font = 14
    impactyr_font =  11
    
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
    col_fill = 'peachpuff'     # fill color
    col_curyear = 'k'     # color of current impact year
    col_2019marker = 'k'        # marker for 2019
    ub_alpha = 0.5
    
    #========== IMPACT YEARS ==========#
    
    # marker size
    markersize = 4
    
    # marker style
    markerstyle = 'o'
    
    #========== AXII ==========#
    
    col_grid = '0.8'     # color background grid
    style_grid = 'dashed'     # style background grid
    lw_grid = 0.5     # lineweight background grid
    
    col_bis = 'black'     # color bisector
    style_bis = '--'     # style bisector
    lw_bis = 1     # lineweight bisector
    
    # x extents for timeseries
    xmin_ts = 1900 # xmin
    xmax_ts = 2100 # xmax
    
    # x ticks/labels for timeseries
    xticks_ts = np.arange(1900,2125,25)
    xtick_labels_ts = [None,None,1950,None,2000,None,2050,None,2100]
    
    # x axis label for timeseries
    xlabel_ts = 'Years'
    xlabel_xpos_ts = 0.33
    xlabel_ypos_ts = 0.05
    
    # x extents for scaling
    xmin_sc = -0.25 # xmin
    xmax_sc = 5 # xmax
    
    # x ticks/labels for scaling
    xticks_sc = np.arange(0,6)
    xtick_labels_sc = np.arange(0,6)
    
    # x axis label for scaling
    xlabel_sc = 'Air temperature (surface,°C)'
    xlabel_xpos_sc = 0.66
    xlabel_ypos_sc = 0.05
    
    #========== LEGEND ==========#
    
    # labels
    lab_pimean = 'Pre-industrial control'
    lab_histmean = 'Historical'
    lab_26mean = 'RCP 2.6'
    lab_60mean = 'RCP 6.0'
    lab_85mean = 'RCP 8.5'
    lab_fill = 'Range'
    
    # bbox
    x0 = 0.4275
    y0 = 1.15
    xlen = 1.3
    ylen = 0.9
    
    # space between entries
    legend_entrypad = 0.5
    
    # length per entry
    legend_entrylen = 0.75
    
    # height/width between panels
    yspace = 0.3
    xspace = 0.3
    
    # figure adjustments for timeseries
    for ax in [ax1,ax2,ax3]:
        lettercount += 1
        ax.set_title(letters[lettercount-1],loc='left',fontsize=title_font,fontweight='bold')
        ax.set_xlim(xmin_ts,xmax_ts)
        ax.xaxis.set_ticks(xticks_ts)
        ax.xaxis.set_ticklabels(xtick_labels_ts)
        ax.tick_params(labelsize=tick_font,axis="x",direction="in", left="off",labelleft="on")
        ax.tick_params(labelsize=tick_font,axis="y",direction="in")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.grid(color=col_grid, linestyle=style_grid, linewidth=lw_grid)
        ax.xaxis.grid(color=col_grid, linestyle=style_grid, linewidth=lw_grid)
        ax.set_axisbelow(True) 
    
    # figure adjustments for scaling
    for ax in [ax4,ax5,ax6]:
        lettercount += 1
        ax.set_title(letters[lettercount-1],loc='left',fontsize=title_font,fontweight='bold')
        ax.set_xlim(xmin_sc,xmax_sc)
        ax.xaxis.set_ticks(xticks_sc)
        ax.xaxis.set_ticklabels(xtick_labels_sc)
        ax.tick_params(labelsize=tick_font,axis="x",direction="in", left="off",labelleft="on")
        ax.tick_params(labelsize=tick_font,axis="y",direction="in")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.grid(color=col_grid, linestyle=style_grid, linewidth=lw_grid)
        ax.xaxis.grid(color=col_grid, linestyle=style_grid, linewidth=lw_grid)
        ax.set_axisbelow(True) 
        
    #==============================================================================
    # WATERTEMP TIMESERIES PLOTTING - ax1
    #==============================================================================
    
    # define time variables
    time_pre = np.arange(1661,1861,1)
    time_his = np.arange(1861,2006,1)
    time_fut = np.arange(2005,2099,1)
    
    ymin = -1   # ymin
    ymax = 5    # ymax
    
    # y axis label
    ylabel = 'Lake temperature anomaly (°C)'
    
    # load data 
    h,ax1 = ts_plotter(wt_mmm,ax1,time_pre,time_his,time_fut,\
        lw_mean,col_pimean,col_pifill,col_histmean,col_histfill,\
        col_rcp26mean,col_rcp26fill,col_rcp60mean,col_rcp60fill,\
        col_rcp85mean,col_rcp85fill,ub_alpha)
    
    ax1.set_ylim(ymin,ymax)
    
    ax1.set_ylabel(ylabel, va='center', rotation='vertical', fontsize=axis_font, labelpad=10)
    
    #==============================================================================
    # ICETHICK TIMESERIES PLOTTING - ax2
    #==============================================================================
    
    ymin = -0.4   # ymin
    ymax = 0.1    # ymax
    
    # y axis label
    ylabel = 'Ice thickness anomaly (m)'
    
    # load data 
    h,ax2 = ts_plotter(it_mmm,ax2,time_pre,time_his,time_fut,\
        lw_mean,col_pimean,col_pifill,col_histmean,col_histfill,\
        col_rcp26mean,col_rcp26fill,col_rcp60mean,col_rcp60fill,\
        col_rcp85mean,col_rcp85fill,ub_alpha)

    ax2.set_ylim(ymin,ymax)
    
    ax2.set_ylabel(ylabel, va='center', rotation='vertical', fontsize=axis_font, labelpad=10)
    
    #==============================================================================
    # ICEDUR TIMESERIES PLOTTING - ax3
    #==============================================================================
    
    # define time variables
    time_pre = np.arange(1661,1860,1)
    time_his = np.arange(1861,2005,1)
    time_fut = np.arange(2005,2098,1)
    
    ymin = -60   # ymin
    ymax = 10    # ymax
    
    # axis labels
    xlabel = 'Time'
    ylabel = 'Ice duration anomaly (days)'
    
    # load data 
    h,ax3 = ts_plotter(ii_mmm,ax3,time_pre,time_his,time_fut,\
        lw_mean,col_pimean,col_pifill,col_histmean,col_histfill,\
        col_rcp26mean,col_rcp26fill,col_rcp60mean,col_rcp60fill,\
        col_rcp85mean,col_rcp85fill,ub_alpha)

    ax3.set_ylim(ymin,ymax)
    
    ax3.set_xlabel(xlabel, va='center', rotation='horizontal', fontsize=axis_font, labelpad=20)
    ax3.set_ylabel(ylabel, va='center', rotation='vertical', fontsize=axis_font, labelpad=20)
    
    #==============================================================================
    # WATERTEMP SCALING PLOTTING - ax4
    #==============================================================================
    
    ymin = -1   # ymin
    ymax = 5    # ymax
    
    # historical
    ax4.plot(wt_air_ens['rcp26'][0].values[201:345],wt_ens['rcp26'][0].values[201:345], lw=lw_mean, color='0.5', zorder=4, label='Historical')
    
    
    # rcp26
    ax4.plot(wt_air_ens['rcp26'][0].values[201:],wt_ens['rcp26'][0].values[201:], lw=lw_mean, color=col_rcp26mean, zorder=3, label='RCP 2.6')
    ax4.fill_between(wt_air_ens['rcp26'][0].values[201:], (wt_ens['rcp26'][1].values[201:]),\
                            (wt_ens['rcp26'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0) 
    
    ax4.plot(wt_iy['rcp26_x2030'], wt_iy['rcp26_y2030'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax4.annotate("2030", (wt_iy['rcp26_x2030'], wt_iy['rcp26_y2030']),xytext=(wt_iy['rcp26_x2030'], wt_iy['rcp26_y2030']-0.3),\
                 color=col_rcp26mean,fontsize=impactyr_font,zorder=5)
    ax4.plot(wt_iy['rcp26_x2050'], wt_iy['rcp26_y2050'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)   
    ax4.annotate("2050, 2090", (wt_iy['rcp26_x2050'], wt_iy['rcp26_y2050']),xytext=(1.1, 1.75),va = "bottom", ha="center",\
                 arrowprops=dict(color=col_rcp26mean,connectionstyle='angle,angleA=180,angleB=90',arrowstyle=']-'),\
        color=col_rcp26mean,fontsize=impactyr_font,zorder=5)                
    ax4.plot(wt_iy['rcp26_x2100'], wt_iy['rcp26_y2100'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)

    
    # rcp60
    ax4.plot(wt_air_ens['rcp60'][0].values[201:],wt_ens['rcp60'][0].values[201:], lw=lw_mean, color=col_rcp60mean, zorder=2, label='RCP 6.0') 
    ax4.fill_between(wt_air_ens['rcp60'][0].values[201:], (wt_ens['rcp60'][1].values[201:]),\
                            (wt_ens['rcp60'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0)
    
    ax4.plot(wt_iy['rcp60_x2030'], wt_iy['rcp60_y2030'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax4.annotate("2030", (wt_iy['rcp60_x2030'], wt_iy['rcp60_y2030']),xytext=(0.8, 1.3),va = "bottom", ha="center",\
                 arrowprops=dict(color=col_rcp60mean,connectionstyle='angle,angleA=180,angleB=90',arrowstyle=']-'),\
        color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                           
    ax4.plot(wt_iy['rcp60_x2050'], wt_iy['rcp60_y2050'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)  
    ax4.annotate("2050", (wt_iy['rcp60_x2050'], wt_iy['rcp60_y2050']),xytext=(wt_iy['rcp60_x2050']+0.075, wt_iy['rcp60_y2050']-0.18),\
                 color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                           
    ax4.plot(wt_iy['rcp60_x2100'], wt_iy['rcp60_y2100'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)  
    ax4.annotate("2090", (wt_iy['rcp60_x2100'], wt_iy['rcp60_y2100']),xytext=(wt_iy['rcp60_x2100'], wt_iy['rcp60_y2100']-0.25),\
                 color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                                          
    
    
    # rcp85                        
    ax4.plot(wt_air_ens['rcp85'][0].values[201:],wt_ens['rcp85'][0].values[201:], lw=lw_mean, color=col_rcp85mean, zorder=1, label='RCP 8.5')
    ax4.fill_between(wt_air_ens['rcp85'][0].values[201:], (wt_ens['rcp85'][1].values[201:]),\
                            (wt_ens['rcp85'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0)    
    
    ax4.plot(wt_iy['rcp85_x2030'], wt_iy['rcp85_y2030'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax4.annotate("2030", (wt_iy['rcp85_x2030'], wt_iy['rcp85_y2030']),xytext=(wt_iy['rcp85_x2030']+0.07, wt_iy['rcp85_y2030']-0.2),\
                 color=col_rcp85mean,fontsize=impactyr_font,zorder=5)
    ax4.plot(wt_iy['rcp85_x2050'], wt_iy['rcp85_y2050'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize)                   
    ax4.annotate("2050", (wt_iy['rcp85_x2050'], wt_iy['rcp85_y2050']),xytext=(wt_iy['rcp85_x2050'], wt_iy['rcp85_y2050']-0.25),\
                 color=col_rcp85mean,fontsize=impactyr_font,zorder=5)                    
    ax4.plot(wt_iy['rcp85_x2100'], wt_iy['rcp85_y2100'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax4.annotate("2090", (wt_iy['rcp85_x2100'], wt_iy['rcp85_y2100']),xytext=(wt_iy['rcp85_x2100'], wt_iy['rcp85_y2100']-0.35),\
                 color=col_rcp85mean,fontsize=impactyr_font,zorder=5)                    
    
    # plot impact years for today and historical sample years
    ax4.plot(wt_iy['rcp85_x2020'], wt_iy['rcp85_y2020'], color=col_2019marker,zorder=4,marker=markerstyle,markersize=markersize) 
    ax4.annotate("2020", (wt_iy['rcp85_x2020'], wt_iy['rcp85_y2020']),xytext=(wt_iy['rcp85_x2020']-0.1, wt_iy['rcp85_y2020']-0.35),\
                 color=col_curyear,fontsize=impactyr_font,zorder=5)                    
    
    ax4.plot(wt_iy['rcp85_x2000'], wt_iy['rcp85_y2000'], color=col_histmean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax4.annotate("2000", (wt_iy['rcp85_x2000'], wt_iy['rcp85_y2000']),xytext=(wt_iy['rcp85_x2000'], wt_iy['rcp85_y2000']-0.3),\
                 color=col_histmean,fontsize=impactyr_font,zorder=5)                    
    
    ax4.plot(wt_iy['rcp85_x1980'], wt_iy['rcp85_y1980'], color=col_histmean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax4.annotate("1980", (wt_iy['rcp85_x1980'], wt_iy['rcp85_y1980']),xytext=(wt_iy['rcp85_x1980'], wt_iy['rcp85_y1980']-0.3),\
                 color=col_histmean,fontsize=impactyr_font,zorder=5)                    
    
    # bisector
    ax4.plot([0,1,2,3,4,4.5],[0,1,2,3,4,4.5], color=col_bis,zorder=4,linestyle=style_bis,linewidth=lw_bis)
    
    ax4.set_ylim(ymin,ymax)
    
    
    #==============================================================================
    # ICETHICK SCALING PLOTTING - ax5
    #==============================================================================

    ymin = -0.4   # ymin
    ymax = 0.1    # ymax
    
    # historical
    ax5.plot(it_air_ens['rcp26'][0].values[201:345],it_ens['rcp26'][0].values[201:345], lw=lw_mean, color='0.5', zorder=4, label='Historical')
    
    
    # rcp26
    ax5.plot(it_air_ens['rcp26'][0].values[201:],it_ens['rcp26'][0].values[201:], lw=lw_mean, color=col_rcp26mean, zorder=3, label='RCP 2.6')
    ax5.fill_between(it_air_ens['rcp26'][0].values[201:], (it_ens['rcp26'][1].values[201:]),\
                            (it_ens['rcp26'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0) 
    ax5.plot(it_iy['rcp26_x2030'], it_iy['rcp26_y2030'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax5.annotate("2030", (it_iy['rcp26_x2030'], it_iy['rcp26_y2030']),xytext=(it_iy['rcp26_x2030'], it_iy['rcp26_y2030']-0.03),va="bottom", ha="center",\
                     color=col_rcp26mean,fontsize=impactyr_font,zorder=5)
    ax5.plot(it_iy['rcp26_x2050'], it_iy['rcp26_y2050'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)   
    ax5.annotate("2050,\n 2090", (it_iy['rcp26_x2050'], it_iy['rcp26_y2050']),xytext=(it_iy['rcp26_x2050']+0.4,it_iy['rcp26_y2050']-0.05),va="bottom", ha="center",\
                arrowprops=dict(color=col_rcp26mean,connectionstyle='angle,angleA=180,angleB=90',arrowstyle=']-'),color=col_rcp26mean,fontsize=impactyr_font,zorder=5)                
    ax5.plot(it_iy['rcp26_x2100'], it_iy['rcp26_y2100'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)

    
    # rcp60
    ax5.plot(it_air_ens['rcp60'][0].values[201:],it_ens['rcp60'][0].values[201:], lw=lw_mean, color=col_rcp60mean, zorder=2, label='RCP 6.0') 
    ax5.fill_between(it_air_ens['rcp60'][0].values[201:], (it_ens['rcp60'][1].values[201:]),\
                            (it_ens['rcp60'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0)
    ax5.plot(it_iy['rcp60_x2030'], it_iy['rcp60_y2030'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax5.annotate("2030", (it_iy['rcp60_x2030'], it_iy['rcp60_y2030']),xytext=(it_iy['rcp60_x2030']-0.05,it_iy['rcp60_y2030']+0.025),va="bottom", ha="center",\
                color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                           
                       
    ax5.plot(it_iy['rcp60_x2050'], it_iy['rcp60_y2050'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)  
    ax5.annotate("2050", (it_iy['rcp60_x2050'], it_iy['rcp60_y2050']),xytext=(it_iy['rcp60_x2050']+0.2, it_iy['rcp60_y2050']+0.02),va="bottom", ha="center",\
                color=col_rcp60mean,fontsize=impactyr_font,zorder=5)
    ax5.plot(it_iy['rcp60_x2100'], it_iy['rcp60_y2100'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)  
    ax5.annotate("2090", (it_iy['rcp60_x2100'], it_iy['rcp60_y2100']),xytext=(it_iy['rcp60_x2100'], it_iy['rcp60_y2100']+.02),\
                color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                                          
    
    
    # rcp85                        
    ax5.plot(it_air_ens['rcp85'][0].values[201:],it_ens['rcp85'][0].values[201:], lw=lw_mean, color=col_rcp85mean, zorder=1, label='RCP 8.5')
    ax5.fill_between(it_air_ens['rcp85'][0].values[201:], (it_ens['rcp85'][1].values[201:]),\
                            (it_ens['rcp85'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0)
    ax5.plot(it_iy['rcp85_x2030'], it_iy['rcp85_y2030'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax5.annotate("2030", (it_iy['rcp85_x2030'], it_iy['rcp85_y2030']),xytext=(it_iy['rcp85_x2030']-0.06, it_iy['rcp85_y2030']+0.02),\
                color=col_rcp85mean,fontsize=impactyr_font,zorder=5)
    ax5.plot(it_iy['rcp85_x2050'], it_iy['rcp85_y2050'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize)                   
    ax5.annotate("2050", (it_iy['rcp85_x2050'], it_iy['rcp85_y2050']),xytext=(it_iy['rcp85_x2050'], it_iy['rcp85_y2050']-0.03),\
                color=col_rcp85mean,fontsize=impactyr_font,zorder=5)                    
    ax5.plot(it_iy['rcp85_x2100'], it_iy['rcp85_y2100'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax5.annotate("2090", (it_iy['rcp85_x2100'], it_iy['rcp85_y2100']),xytext=(it_iy['rcp85_x2100'], it_iy['rcp85_y2100']+0.02),\
                color=col_rcp85mean,fontsize=impactyr_font,zorder=5)                    
    
    
    # plot impact years for today and historical sample years
    ax5.plot(it_iy['rcp85_x2020'], it_iy['rcp85_y2020'], color=col_2019marker,zorder=4,marker=markerstyle,markersize=markersize) 
    ax5.annotate("2020", (it_iy['rcp85_x2020'], it_iy['rcp85_y2020']),xytext=(it_iy['rcp85_x2020']-0.275, it_iy['rcp85_y2020']-0.02),\
                color=col_curyear,fontsize=impactyr_font,zorder=5)                    
    
    ax5.plot(it_iy['rcp85_x2000'], it_iy['rcp85_y2000'], color=col_histmean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax5.annotate("2000", (it_iy['rcp85_x2000'], it_iy['rcp85_y2000']),xytext=(it_iy['rcp85_x2000'], it_iy['rcp85_y2000']+0.02),\
                color=col_histmean,fontsize=impactyr_font,zorder=5)                    
    
    ax5.plot(it_iy['rcp85_x1980'], it_iy['rcp85_y1980'], color=col_histmean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax5.annotate("1980", (it_iy['rcp85_x1980'], it_iy['rcp85_y1980']),xytext=(it_iy['rcp85_x1980'], it_iy['rcp85_y1980']+0.02),\
                color=col_histmean,fontsize=impactyr_font,zorder=5)                    
    
    ax5.set_ylim(ymin,ymax)
    
    #==============================================================================
    # ICEDUR SCALING PLOTTING - ax6
    #==============================================================================
    
    ymin = -60   # ymin
    ymax = 10    # ymax
    
    # x-axis label
    xlabel = 'Air temperature anomaly (°C)'
    
    # historical
    ax6.plot(ii_air_ens['rcp26'][0].values[201:345],ii_ens['rcp26'][0].values[201:345], lw=lw_mean, color='0.5', zorder=4, label='Historical')
    
    
    # rcp26
    ax6.plot(ii_air_ens['rcp26'][0].isel(dim_0=slice(0,-3)).values[201:],ii_ens['rcp26'][0].values[201:], lw=lw_mean, color=col_rcp26mean, zorder=3, label='RCP 2.6')
    ax6.fill_between(ii_air_ens['rcp26'][0].isel(dim_0=slice(0,-3)).values[201:], (ii_ens['rcp26'][1].values[201:]),\
                            (ii_ens['rcp26'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0) 
    ax6.plot(ii_iy['rcp26_x2030'], ii_iy['rcp26_y2030'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax6.annotate("2030", (ii_iy['rcp26_x2030'], ii_iy['rcp26_y2030']),xytext=(ii_iy['rcp26_x2030']-0.2, ii_iy['rcp26_y2030']-3),color=col_rcp26mean,fontsize=impactyr_font,zorder=5)
    ax6.plot(ii_iy['rcp26_x2050'], ii_iy['rcp26_y2050'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)   
    ax6.annotate("2050, 2090", (ii_iy['rcp26_x2050'], ii_iy['rcp26_y2050']),xytext=(ii_iy['rcp26_x2050'],ii_iy['rcp26_y2050']-6),va = "bottom", ha="center",arrowprops=dict(color=col_rcp26mean,connectionstyle='angle,angleA=180,angleB=90',arrowstyle=']-'),color=col_rcp26mean,fontsize=impactyr_font,zorder=5)                
    ax6.plot(ii_iy['rcp26_x2100'], ii_iy['rcp26_y2100'], color=col_rcp26mean,zorder=4,marker=markerstyle,markersize=markersize)
    
    
    # rcp60
    ax6.plot(ii_air_ens['rcp60'][0].isel(dim_0=slice(0,-3)).values[201:],ii_ens['rcp60'][0].values[201:], lw=lw_mean, color=col_rcp60mean, zorder=2, label='RCP 6.0') 
    ax6.fill_between(ii_air_ens['rcp60'][0].isel(dim_0=slice(0,-3)).values[201:], (ii_ens['rcp60'][1].values[201:]),\
                            (ii_ens['rcp60'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0)
    ax6.plot(ii_iy['rcp60_x2030'], ii_iy['rcp60_y2030'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax6.annotate("2030", (ii_iy['rcp60_x2030'], ii_iy['rcp60_y2030']),xytext=(ii_iy['rcp60_x2030']+0.025,ii_iy['rcp60_y2030']+2.5),va = "bottom", ha="center",color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                           
    ax6.plot(ii_iy['rcp60_x2050'], ii_iy['rcp60_y2050'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)  
    ax6.annotate("2050", (ii_iy['rcp60_x2050'], ii_iy['rcp60_y2050']),xytext=(ii_iy['rcp60_x2050']+0.05, ii_iy['rcp60_y2050']+1),color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                           
    ax6.plot(ii_iy['rcp60_x2100'], ii_iy['rcp60_y2100'], color=col_rcp60mean,zorder=4,marker=markerstyle,markersize=markersize)  
    ax6.annotate("2090", (ii_iy['rcp60_x2100'], ii_iy['rcp60_y2100']),xytext=(ii_iy['rcp60_x2100'], ii_iy['rcp60_y2100']+3),color=col_rcp60mean,fontsize=impactyr_font,zorder=5)                                          
    
    
    # rcp85                        
    ax6.plot(ii_air_ens['rcp85'][0].isel(dim_0=slice(0,-3)).values[201:],ii_ens['rcp85'][0].values[201:], lw=lw_mean, color=col_rcp85mean, zorder=1, label='RCP 8.5')
    ax6.fill_between(ii_air_ens['rcp85'][0].isel(dim_0=slice(0,-3)).values[201:], (ii_ens['rcp85'][1].values[201:]),\
                            (ii_ens['rcp85'][2].values[201:]),\
                            lw=0.1, color=col_fill, zorder=0)
    ax6.plot(ii_iy['rcp85_x2030'], ii_iy['rcp85_y2030'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize)
    ax6.annotate("2030", (ii_iy['rcp85_x2030'], ii_iy['rcp85_y2030']),xytext=(ii_iy['rcp85_x2030'], ii_iy['rcp85_y2030']+1.5),color=col_rcp85mean,fontsize=impactyr_font,zorder=5)
    ax6.plot(ii_iy['rcp85_x2050'], ii_iy['rcp85_y2050'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize)                   
    ax6.annotate("2050", (ii_iy['rcp85_x2050'], ii_iy['rcp85_y2050']),xytext=(ii_iy['rcp85_x2050'], ii_iy['rcp85_y2050']+3),color=col_rcp85mean,fontsize=impactyr_font,zorder=5)                    
    ax6.plot(ii_iy['rcp85_x2100'], ii_iy['rcp85_y2100'], color=col_rcp85mean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax6.annotate("2090", (ii_iy['rcp85_x2100'], ii_iy['rcp85_y2100']),xytext=(ii_iy['rcp85_x2100'], ii_iy['rcp85_y2100']+3),color=col_rcp85mean,fontsize=impactyr_font,zorder=5)                    
    
    
    # plot impact years for today and historical sample years
    ax6.plot(ii_iy['rcp85_x2020'], ii_iy['rcp85_y2020'], color=col_2019marker,zorder=4,marker=markerstyle,markersize=markersize) 
    ax6.annotate("2020", (ii_iy['rcp85_x2020'], ii_iy['rcp85_y2020']),xytext=(ii_iy['rcp85_x2020']-0.3, ii_iy['rcp85_y2020']-3),color=col_curyear,fontsize=impactyr_font,zorder=5)                    
    
    ax6.plot(ii_iy['rcp85_x2000'], ii_iy['rcp85_y2000'], color=col_histmean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax6.annotate("2000", (ii_iy['rcp85_x2000'], ii_iy['rcp85_y2000']),xytext=(ii_iy['rcp85_x2000'], ii_iy['rcp85_y2000']+1),color=col_histmean,fontsize=impactyr_font,zorder=5)                    
    
    ax6.plot(ii_iy['rcp85_x1980'], ii_iy['rcp85_y1980'], color=col_histmean,zorder=4,marker=markerstyle,markersize=markersize) 
    ax6.annotate("1980", (ii_iy['rcp85_x1980'], ii_iy['rcp85_y1980']),xytext=(ii_iy['rcp85_x1980'], ii_iy['rcp85_y1980']+1),color=col_histmean,fontsize=impactyr_font,zorder=5)                    
    
    ax6.set_ylim(ymin,ymax)
    
    ax6.set_xlabel(xlabel, va='center', rotation='horizontal', fontsize=axis_font, labelpad=20)
    
    # legend
    legendcols = [col_pimean,col_histmean,col_rcp26mean,col_rcp60mean,col_rcp85mean,col_fill]
    handles = [Line2D([0],[0],linestyle='-',lw=2,color=legendcols[0]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[1]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[2]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[3]),\
               Line2D([0],[0],linestyle='-',lw=2,color=legendcols[4]),\
               Rectangle((0,0),1,1,color=legendcols[5])]
    labels= [lab_pimean,lab_histmean,lab_26mean,lab_60mean,lab_85mean,lab_fill]
    ax1.legend(handles, labels, bbox_to_anchor=(x0, y0, xlen, ylen), loc=3,   #bbox: (x, y, width, height)
               ncol=6,fontsize=legend_font, mode="expand", borderaxespad=0.,\
               frameon=False, columnspacing=0.05, handlelength=legend_entrylen, handletextpad=legend_entrypad)
    
    plt.show()
    
    # save figure
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        f.savefig(outDIR+'/f4.png',bbox_inches='tight',dpi=500)

