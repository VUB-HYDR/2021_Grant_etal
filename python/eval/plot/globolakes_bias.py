#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:18:37 2019

@author: Luke
"""
#==============================================================================
#SUMMARY
#==============================================================================


# 14 January 2020


#==============================================================================
#IMPORT
#==============================================================================


import xarray as xr
import os
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec


#==============================================================================
#FUNCTIONS
#==============================================================================


def reader(file):
    ds = xr.open_dataset(file, decode_times=False)
    da = ds.lmlt.mean(dim='time')#.squeeze(dim='time',drop=True)
    return da

def globolakes_plot(valDIR,outDIR,flag_svplt):
    
    #==============================================================================
    #INITIALIZE
    #==============================================================================
    

    os.chdir(valDIR)
    # see evernote re meeting with wim on 13 for definition per file below
    file = 'era5-land_globolakes_lwst_timmean_bias.nc'
    
    data = reader(file)
    lat = data.latitude.values
    lon = data.longitude.values
    data = data.values
    
    #==============================================================================
    #DEFINE PLOT SETTINGS
    #==============================================================================
        
    # font settings
    title_font = 13
    cbtitle_font = 13
    tick_font = 11
    arrow_font = 13
    
    # continent fill color
    col_cont='white'
    
    # ocean fill color
    col_ocean='whitesmoke'
    
    # zero change color
    col_zero='gray'
    
    #========== COLORBAR ==========#
    
    # colorbar colormap setup
    cmap_whole = plt.cm.get_cmap('RdBu')
    cmap55 = cmap_whole(0.01)   
    cmap50 = cmap_whole(0.05)   #red
    cmap45 = cmap_whole(0.1)
    cmap40 = cmap_whole(0.15)
    cmap35 = cmap_whole(0.2)
    cmap30 = cmap_whole(0.25)
    cmap25 = cmap_whole(0.3)
    cmap20 = cmap_whole(0.35)
    cmap10 = cmap_whole(0.4)
    cmap0 = col_zero
    cmap_5 = cmap_whole(0.55)
    cmap_10 = cmap_whole(0.6)
    cmap_20 = cmap_whole(0.65)
    cmap_25 = cmap_whole(0.7)
    cmap_30 = cmap_whole(0.75)
    cmap_35 = cmap_whole(0.8)
    cmap_40 = cmap_whole(0.85)
    cmap_45 = cmap_whole(0.9)    
    cmap_50 = cmap_whole(0.95)  #blue
    cmap_55 = cmap_whole(0.99)
     
    colorlist = [cmap_50,cmap_45,cmap_40,cmap_35,cmap_30,cmap_25,cmap_20,cmap0,\
                 cmap20,cmap25,cmap30,cmap35,cmap40,cmap45,cmap50]
    
    cmap = mpl.colors.ListedColormap(colorlist,N=len(colorlist))
    
    # set color of over/under arrows on colorbar
    cmap.set_over(cmap55)
    cmap.set_under(cmap_55)
    
    # colorbar args
    values = [-7,-6,-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5,6,7]
    tick_locs = [-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7]
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    # bbox
    cb_x0 = 0.0395
    cb_y0 = 0.3
    cb_xlen = 0.5
    cb_ylen = 0.01
    
    # colorbar label
    cblabel = 'Bias in mixed layer temperature (°C)'
     
    #========== ARROWS ==========#
    
    # blue arrow label
    bluelabel = 'Colder'
    x0_bluelab = 0.25
    y0_bluelab = -2.9
    
    # blue arrow
    x0_bluearr = 0.495
    y0_bluearr = -3.5
    xlen_bluearr = -0.44
    ylen_bluearr = 0
    
    # red arrow label
    redlabel = 'Warmer'
    x0_redlab = 0.75
    y0_redlab = -2.9
    
    # red arrow
    x0_redarr = 0.505
    y0_redarr = -3.5
    xlen_redarr = 0.44
    ylen_redarr = 0
    
    # general
    arrow_width = 0.25
    arrow_linew = 0.1
    arrow_headwidth = 0.5
    arrow_headlength = 0.06

    #==============================================================================
    #INITIALIZE PLOTTING
    #==============================================================================
    
    lon, lat = np.meshgrid(lon, lat)
    
    f = plt.figure(figsize=(15,15))
    
    # maps rect for left column, rect=[left, bottom, right, top]
    m_left = 0
    m_bottom = 0.25
    m_right = 0.55
    m_top = 0.65
    m_rect = [m_left, m_bottom, m_right, m_top]
    
    gs1 = gridspec.GridSpec(1,1)
    ax1 = f.add_subplot(gs1[0])    

    gs1.tight_layout(figure=f, rect=m_rect, h_pad=1)
    
    # histogram rect for right column, rect=[left, bottom, right, top]
    h_left = 0.6
    h_bottom = 0.275
    h_right = 1.0
    h_top = 0.575
    h_rect = [h_left, h_bottom, h_right, h_top]
    
    gs2 = gridspec.GridSpec(1,1)
    ax2 = f.add_subplot(gs2[0])    

    gs2.tight_layout(figure=f, rect=h_rect, h_pad=3.225)
    
    # map
    m = Basemap(llcrnrlon=-170, llcrnrlat=-60, urcrnrlon=180, urcrnrlat=90, suppress_ticks=False);
    m.ax = ax1
    m.drawcoastlines(linewidth=0.05);
    m.drawmapboundary(linewidth=0.15,fill_color=col_ocean);
    m.fillcontinents(color=col_cont);
    m.pcolormesh(lon,lat,data,latlon=True,\
                 cmap=cmap,norm=norm,vmax=values[-1],vmin=values[0],zorder=3)
    
    ax1.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=False, top=False, left=False, right=False, color='0.2',\
                     labelcolor='0.2',width=0.4,direction="in",length=2.5)
    ax1.set_title('a',loc='left',fontsize=title_font,fontweight='bold')
    ax1.spines['bottom'].set_color('0.2')
    ax1.spines['bottom'].set_linewidth(0.4)
    ax1.spines['top'].set_color('0.2')
    ax1.spines['top'].set_linewidth(0.4)
    ax1.xaxis.label.set_color('0.2')
    ax1.spines['left'].set_color('0.2')
    ax1.spines['left'].set_linewidth(0.4)
    ax1.spines['right'].set_color('0.2')
    ax1.spines['right'].set_linewidth(0.4)
    ax1.yaxis.label.set_color('0.2')
        
    # colorbar setup
    cbaxes = f.add_axes([cb_x0, cb_y0, cb_xlen, cb_ylen])
    cb = mpl.colorbar.ColorbarBase(ax=cbaxes, cmap=cmap,
                                   norm=norm,
                                   spacing='proportional',
                                   orientation='horizontal',
                                   extend='both',
                                   ticks=tick_locs)
    cb.set_label(cblabel,size=cbtitle_font)
    cb.ax.xaxis.set_label_position('top');
    cb.ax.tick_params(labelcolor='0.2',labelsize=tick_font,color='0.2',\
                      length=2.5,width=0.35,direction='out'); 
    cb.ax.set_xticklabels([r'$\leq$-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6',r'7$\leq$'])
    cb.outline.set_edgecolor('0.2')
    cb.outline.set_linewidth(0.4)
    
    # arrows
    plt.text(x0_redlab, y0_redlab, redlabel, size=arrow_font, ha='center', va='center')
    plt.text(x0_bluelab, y0_bluelab, bluelabel, size=arrow_font, ha='center', va='center')
    
    plt.arrow(x0_redarr, y0_redarr, xlen_redarr, ylen_redarr, width=arrow_width, linewidth=arrow_linew,\
              shape='right', head_width=arrow_headwidth, head_length=arrow_headlength,\
              facecolor=cmap40, edgecolor='k', clip_on=False)
    plt.arrow(x0_bluearr, y0_bluearr, xlen_bluearr, ylen_bluearr, width=arrow_width, linewidth=arrow_linew,\
              shape='left', head_width=arrow_headwidth, head_length=arrow_headlength,\
              facecolor=cmap_40, edgecolor='k', clip_on=False)
    
    # histogram
    ax2.hist(data.flatten(),bins=20,rwidth=0.9)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_ylabel('Grid cells [-]',fontsize=title_font)
    ax2.set_title('b',loc='left',fontsize=title_font,fontweight='bold')
    ax2.set_xlabel('Bias (°C)',fontsize=title_font)

    plt.show()
    
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        #save figure
        f.savefig(outDIR+'/'+'era5-land_globolakes_bias.png',bbox_inches='tight',dpi=500)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
