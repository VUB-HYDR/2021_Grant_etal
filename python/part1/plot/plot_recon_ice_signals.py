#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:18:37 2019

@author: Luke
"""

#==============================================================================
#SUMMARY
#==============================================================================

# 18 May 2020

# plots era5-land ice cover sigs

#==============================================================================
#IMPORT
#==============================================================================

import xarray as xr
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

#==============================================================================
#FUNCTIONS
#==============================================================================

def reader(file):
    if 'start' in file:
        da = xr.open_dataset(file,decode_times=False).icestart.squeeze(dim='time')
    if 'end' in file:
        da = xr.open_dataset(file,decode_times=False).iceend.squeeze(dim='time')
    if 'dur' in file:
        da = xr.open_dataset(file,decode_times=False).iceduration.squeeze(dim='time')
    return da

def ice_plot(iceDIR,outDIR,flag_svplt):
    
    #==============================================================================
    #SETTINGS
    #==============================================================================
    
    # continent fill
    col_cont='white'
    
    # ocean fill
    col_ocean='whitesmoke'
    
    # zero change color
    col_zero='gray'
    
    # list of ice index titles
    ice_titles = ['Ice onset', 'Ice break-up', 'Ice duration']
    
    letters = ['a', 'b', 'c']
    
    # font settings
    title_font = 15
    cbtitle_font = 15
    tick_font = 12
    arrow_font = 14
    
    # bbox (arrow plot relative to this axis)
    cb_x0 = 0.225
    cb_y0 = 0.26
    cb_xlen = 0.55
    cb_ylen = 0.015
    
    # colorbar label
    cblabel = 'Δ reconstructed ice index (days)'
    
    #========== ARROWS ==========#
    
    # blue arrow label
    bluelabel = 'Later date (panels a,b) or longer duration (panel c)'
    x0_bluelab = 0.85
    y0_bluelab = -2.7
    
    # blue arrow
    x0_bluearr = 0.505
    y0_bluearr = -3.3
    xlen_bluearr = 0.8
    ylen_bluearr = 0
    
    # red arrow label
    redlabel = 'Earlier date (panels a,b) or shorter duration (panel c)'
    x0_redlab = 0.14
    y0_redlab = -2.7
    
    # red arrow
    x0_redarr = 0.495
    y0_redarr = -3.3
    xlen_redarr = -0.8
    ylen_redarr = 0
    
    # general
    arrow_width = 0.25
    arrow_linew = 0.1
    arrow_headwidth = 0.5
    arrow_headlength = 0.06
    
    #==============================================================================
    #INITIALIZE DIRECTORIES
    #==============================================================================
    
    os.chdir(iceDIR)
    
    files = []
    for file in sorted(os.listdir(iceDIR)):
        if '.nc' in file:
            files.append(file)
    
    #==============================================================================
    #DATA 
    #==============================================================================
    
    start_plottable = reader(files[2])
    end_plottable = reader(files[1])
    dur_plottable = reader(files[0])
    
    signals = [start_plottable,end_plottable,dur_plottable]
    
    lat = start_plottable.latitude.values
    lon = start_plottable.longitude.values

    
    #==============================================================================
    #PLOT MAIN FIG
    #==============================================================================
    
    # identify colors
    cmap_whole = plt.cm.get_cmap('RdBu_r')
    cmap55 = cmap_whole(0.01)
    cmap50 = cmap_whole(0.05)   #blue
    cmap45 = cmap_whole(0.1)
    cmap40 = cmap_whole(0.15)
    cmap35 = cmap_whole(0.2)
    cmap30 = cmap_whole(0.25)
    cmap25 = cmap_whole(0.3)
    cmap20 = cmap_whole(0.325)
    cmap10 = cmap_whole(0.4)
    cmap5 = cmap_whole(0.475)
    cmap0 = col_zero
    cmap_5 = cmap_whole(0.525)
    cmap_10 = cmap_whole(0.6)
    cmap_20 = cmap_whole(0.625)
    cmap_25 = cmap_whole(0.7)
    cmap_30 = cmap_whole(0.75)
    cmap_35 = cmap_whole(0.8)
    cmap_40 = cmap_whole(0.85)
    cmap_45 = cmap_whole(0.9)
    cmap_50 = cmap_whole(0.95)  #red
    cmap_55 = cmap_whole(0.99)
    
    colorlist = [cmap_45,cmap_40,cmap_35,cmap_30,cmap_25,cmap_10,cmap0,\
                 cmap5,cmap10,cmap20,cmap30,cmap35,cmap40]
    
    cmap = mpl.colors.ListedColormap(colorlist,N=len(colorlist))
    cmap.set_over(cmap55)
    cmap.set_under(cmap_55)
    
    values = [-30,-25,-20,-15,-10,-5,-1,1,5,10,15,20,25,30]
    
    tick_locs = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]
    
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    #=============================================================================
    #SET PLOTS
    #=============================================================================
    
    f, axes = plt.subplots(1,3,figsize=(15,12));
    
    lon, lat = np.meshgrid(lon, lat)
        
    count = 0
    for rcpmap,ax in zip(signals,axes.flatten()):
        count = count+1
        m = Basemap(projection='npaeqd',round=True,boundinglat=20,\
                    lat_0=80,lon_0=0,resolution='l');
        m.ax = ax
        m.drawcoastlines(linewidth=0.05);
        m.drawmapboundary(linewidth=0.15,fill_color=col_ocean);
        m.fillcontinents(color=col_cont);
        m.pcolormesh(lon,lat,rcpmap,latlon=True,cmap=cmap,norm=norm,vmax=30,vmin=-30,zorder=3)
        if count<=3:
            ax.set_title(ice_titles[count-1],loc='center',fontsize=title_font)
            ax.set_title(letters[count-1],loc='left',fontsize=title_font,fontweight='bold')
    
        
    #==============================================================================
    #COLORBAR
    #==============================================================================
            
    cbax = f.add_axes([cb_x0, cb_y0, cb_xlen, cb_ylen])
    cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                               norm=norm,
                               spacing='proportional',
                               orientation='horizontal',
                               extend='both',
                               ticks=tick_locs)
    cb.set_label(cblabel,size=cbtitle_font)
    cb.ax.xaxis.set_label_position('top');
    cb.ax.tick_params(labelcolor='0.2',labelsize=tick_font,color='0.2',\
                      length=2.5,width=0.35,direction='out'); #change color of ticks?
    cb.ax.set_xticklabels([r'$\leq$-30','-25','-20','-15','-10',\
                           '-5','0','5','10','15','20','25',r'30$\leq$'])
    cb.outline.set_edgecolor('0.2')
    cb.outline.set_linewidth(0.4)
    
    #plot arrows
    bluelabel = 'Later date (panels a,b) or longer duration (panel c)'
    redlabel = 'Earlier date (panels a,b) or shorter duration (panel c)'
    
    plt.text(0.85, -2.8, bluelabel, size=arrow_font, ha='center', va='center')
    plt.text(0.15, -2.8, redlabel, size=arrow_font, ha='center', va='center')
    
    plt.arrow(0.505, -3.5, 0.75, 0, width=0.25, linewidth=0.1, label=bluelabel,\
              shape='right', head_width=0.5, head_length=0.06,\
              facecolor=cmap40, edgecolor='k', clip_on=False)
    plt.arrow(0.495, -3.5, -0.75, 0, width=0.25, linewidth=0.1, label=redlabel,\
              shape='left', head_width=0.5, head_length=0.06,\
              facecolor=cmap_40, edgecolor='k', clip_on=False)
    
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0.15, hspace=0.05)
    
    plt.show()
    
    
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        #save figure
        f.savefig(outDIR+'f1',bbox_inches='tight',dpi=900 )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
