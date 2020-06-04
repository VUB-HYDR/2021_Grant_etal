#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:18:37 2019

@author: Luke
"""

#==============================================================================
#SUMMARY
#==============================================================================

#This script generates ERA5 icedepth (communicated as ice cover start/end/duration)
#signal plots

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
        da = xr.open_dataset(file,decode_times=False).icestart.squeeze(dim='time')#.squeeze(dim='time')
    if 'end' in file:
        da = xr.open_dataset(file,decode_times=False).iceend.squeeze(dim='time')#.squeeze(dim='time')
    if 'dur' in file:
        da = xr.open_dataset(file,decode_times=False).icedur.squeeze(dim='time')#.squeeze(dim='time')
    return da


def amsr_plot(valDIR,outDIR,flag_svplt):
    
    #==============================================================================
    #SETTINGS
    #==============================================================================
    
    # font settings
    title_font = 13
    cbtitle_font = 13
    tick_font = 11
    arrow_font = 13
    
    # list of figure panel ids
    letters = ['a', 'b', 'c',\
               'd', 'e', 'f',\
               'g', 'h', 'i',\
               'j', 'k', 'l']
    
    #==============================================================================
    #INITIALIZE DIRECTORIES
    #==============================================================================
    
    os.chdir(valDIR)
    
    files = []
    for file in sorted(os.listdir(valDIR)):
        if '.nc' in file and 'amsr' in file:
            files.append(file)
    
    #==============================================================================
    #DATA 
    #==============================================================================
    
    start_plottable = reader(files[2])
    end_plottable = reader(files[1])
    dur_plottable = reader(files[0])
    
    data = [start_plottable,end_plottable,dur_plottable]
    
    lat = start_plottable.lat.values
    lon = start_plottable.lon.values
    
    ice_titles = ['Ice onset', 'Ice break-up', 'Ice duration']
    letters = ['a)', 'b)', 'c)']
    
    #==============================================================================
    #PLOT TEST HIST
    #==============================================================================
    
    #hist to see spread of values dictating colorbar range for signal plots
    f, ax1 = plt.subplots(1,3,figsize=(15,5),sharey=True)
    
    count=0
    for ice,ax in zip(data,ax1.flat):
        count=count+1
        ax.hist(ice.values.flatten(),bins=20,range=(-80,80),rwidth=0.9)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(letters[count-1],loc='left',fontsize=title_font, fontweight='bold')
        ax.set_title(ice_titles[count-1],loc='center',fontsize=title_font)
        if count == 2:
            ax.set_xlabel('Bias in ice index (days)',fontsize=title_font,labelpad=5)
    
    f.text(0.0, 0.5, 'Grid cells [-]', va='center', rotation='vertical',fontsize=title_font)
    
    plt.tight_layout()
    plt.show()  

    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        #save figure    
        f.savefig(outDIR+'/'+'era5-land_amsr_bias_histogram.png',bbox_inches='tight',dpi=500)
    
    #==============================================================================
    #PLOT MAIN FIG
    #==============================================================================
        
    # continent fill
    col_cont='white'
    
    # ocean fill
    col_ocean='whitesmoke'
    
    # zero change color
    col_zero='gray'
    
    # list of ice index titles
    ice_titles = ['Ice onset', 'Ice break-up', 'Ice duration']
    
    #========== COLORBAR ==========#
        
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
    
    
    colorlist = [cmap_50,cmap_40,cmap_35,cmap_30,cmap_20,cmap_5,cmap0,\
                 cmap5,cmap20,cmap30,cmap35,cmap40,cmap50]
    
    # declare list of colors for discrete colormap of colorbar
    cmap = mpl.colors.ListedColormap(colorlist,N=len(colorlist))
    
    # set color of over/under arrows in colorbar
    cmap.set_over(cmap55)
    cmap.set_under(cmap_55)
    
    # colorbar args
    values = [-30,-25,-20,-15,-10,-5,5,10,15,20,25,30]
    tick_locs = [-30,-25,-20,-15,-10,0,10,15,20,25,30]
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    #=============================================================================
    #SET PLOTS
    #=============================================================================
    
    f, axes = plt.subplots(1,3,figsize=(15,15));
    
    lon, lat = np.meshgrid(lon, lat)
        
    count = 0
    for ice,ax in zip(data,axes.flatten()):
        count += 1
        m = Basemap(projection='npaeqd',round=True,boundinglat=20,\
                    lat_0=80,lon_0=0,resolution='l');
        m.ax = ax
        m.drawcoastlines(linewidth=0.2);
        m.drawmapboundary(linewidth=0.15);
        m.fillcontinents(color='whitesmoke');
        m.pcolormesh(lon,lat,ice,latlon=True,cmap=cmap,norm=norm,vmax=values[-1],vmin=values[0],zorder=3)
        ax.set_title(letters[count-1],loc='left',fontsize=title_font, fontweight='bold')
        ax.set_title(ice_titles[count-1],loc='center',fontsize=title_font)
    
    #==============================================================================
    #COLORBAR
    #==============================================================================
            
    cbaxes = f.add_axes([0.215, 0.375, 0.6, 0.015])
    cb = mpl.colorbar.ColorbarBase(ax=cbaxes,cmap=cmap,
                                   norm=norm,
                                   spacing='uniform',
                                   orientation='horizontal',
                                   extend='both',
                                   ticks=tick_locs)
    cb.set_label('Bias in ice onset, break-up or duration (days)',size=title_font)
    cb.ax.xaxis.set_label_position('top');
    cb.ax.tick_params(labelcolor='0.2',labelsize=tick_font,color='0.2',\
                      length=2.5,width=0.35,direction='out'); 
    cb.ax.set_xticklabels([r'$\leq$-30','-25','-20','-15','-10',\
                           '0','10','15','20','25',r'30$\leq$'])
    cb.outline.set_edgecolor('0.2')
    cb.outline.set_linewidth(0.4)
    
    #plot arrows
    bluelabel = 'Later date (panels a,b) or longer duration (panel c)'
    redlabel = 'Earlier date (panels a,b) or shorter duration (panel c)'
    
    plt.text(0.8, -2.8, bluelabel, size=arrow_font, ha='center', va='center')
    plt.text(0.2, -2.8, redlabel, size=arrow_font, ha='center', va='center')
    
    plt.arrow(0.510, -3.5, 0.6, 0, width=0.25, linewidth=0.1, label=bluelabel,\
              shape='right', head_width=0.5, head_length=0.06,\
              facecolor=cmap40, edgecolor='k', clip_on=False)
    plt.arrow(0.490, -3.5, -0.6, 0, width=0.25, linewidth=0.1, label=redlabel,\
              shape='left', head_width=0.5, head_length=0.06,\
              facecolor=cmap_40, edgecolor='k', clip_on=False)
    
    plt.subplots_adjust(left=0.175, right=0.85, bottom=0.2, top=0.875, wspace=0.03, hspace=0.05)
    
    plt.show()
    
    
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        #save figure
        f.savefig(outDIR+'/'+'era5-land_amsr_bias.png',bbox_inches='tight',dpi=500)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
