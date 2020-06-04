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

# plots maps of climate change impacts


#==============================================================================
# IMPORT
#==============================================================================


import xarray as xr
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#==============================================================================
# FUNCTION
#==============================================================================


def plot_p3(wt_annual,wt_jja_scaled,iisig_data,lat,lon,outDIR,flag_svplt):
    
    #==============================================================================
    # general
    #==============================================================================
    
    # font settings
    title_font = 15
    cbtitle_font = 15
    tick_font = 12
    arrow_font = 14
    
    # list of figure panel ids
    letters = ['a', 'b', 'c',\
               'd', 'e', 'f',\
               'g', 'h', 'i',\
               'j', 'k', 'l']
    
    lettercount = 0
    
    lon, lat = np.meshgrid(lon, lat)
    
    #==============================================================================
    # initialize
    #==============================================================================
    
    f = plt.figure(figsize=(15,14))
    
    gs1 = gridspec.GridSpec(1,2)
    
    ax1 = f.add_subplot(gs1[0])
    ax2 = f.add_subplot(gs1[1])
    
    # rect=[left, bottom, right, top]
    gs1.tight_layout(figure=f, rect=[0, 0.65, 1, 1])
    
    gs2 = gridspec.GridSpec(1,3)
    
    ax3 = f.add_subplot(gs2[0])
    ax4 = f.add_subplot(gs2[1])
    ax5 = f.add_subplot(gs2[2])
    
    gs2.tight_layout(figure=f, rect=[0, 0.1, 1, 0.65])
    
    #==============================================================================
    # watertemp annual plotting
    #==============================================================================
    
    # continent fill color
    col_cont='white'
    
    # ocean fill color
    col_ocean='whitesmoke'
    
    # zero change color
    col_zero='gray'
    
    #========== COLORBAR ==========#
        
    
    # identify colors
    cmap_whole = plt.cm.get_cmap('RdBu')
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
    
    # declare list of colors for discrete colormap of colorbar
    cmap = mpl.colors.ListedColormap([cmap_45,cmap_35,cmap_25,cmap_10,cmap_5,cmap0,
                                      cmap5,cmap10,cmap20,cmap30,cmap40],N=11)
    
    # set color of over/under arrows in colorbar
    cmap.set_over(cmap50)
    cmap.set_under(cmap_50)
    
    # colorbar args
    values = [-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5]
    tick_locs = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    # colorbar label
    cblabel = 'Δ lake temperature (°C)'
    
    # bbox (arrow plot relative to this axis)
    cb_x0 = 0.063
    cb_y0 = 0.675
    cb_xlen = 0.4
    cb_ylen = 0.015
    
    #========== ARROWS ==========#
    
    # blue arrow label
    bluelabel = 'Colder'
    x0_bluelab = 0.25
    y0_bluelab = -2.8
    
    # blue arrow
    x0_bluearr = 0.495
    y0_bluearr = -3.3
    xlen_bluearr = -0.4
    ylen_bluearr = 0
    
    # red arrow label
    redlabel = 'Warmer'
    x0_redlab = 0.75
    y0_redlab = -2.8
    
    # red arrow
    x0_redarr = 0.505
    y0_redarr = -3.3
    xlen_redarr = 0.4
    ylen_redarr = 0
    
    # general
    arrow_width = 0.25
    arrow_linew = 0.1
    arrow_headwidth = 0.5
    arrow_headlength = 0.06
    
    #========== PLOTTING ==========#
    
    lettercount += 1
    
    m = Basemap(llcrnrlon=-170, llcrnrlat=-60, urcrnrlon=180, urcrnrlat=90, suppress_ticks=False);
    m.ax = ax1
    m.drawcoastlines(linewidth=0.1);
    m.drawmapboundary(fill_color=col_ocean)
    m.fillcontinents(color=col_cont);
    m.pcolormesh(lon,lat,wt_annual,latlon=True,cmap=cmap,norm=norm,vmax=5,vmin=-5,zorder=3)
    ax1.set_title(letters[lettercount-1],loc='left',pad=10,fontsize=title_font,fontweight='bold')
    ax1.set_title('Annual',loc='center',pad=10,fontsize=title_font)
    
    ax1.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=False, top=False, left=False, right=False, color='0.2',\
                     labelcolor='0.2',width=0.4,direction="in",length=2.5)
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
    cbax = f.add_axes([cb_x0, cb_y0, cb_xlen, cb_ylen])
    cb = mpl.colorbar.ColorbarBase(ax=cbax, cmap=cmap,
                                   norm=norm,
                                   spacing='proportional',
                                   orientation='horizontal',
                                   extend='both',
                                   ticks=tick_locs)
    cb.set_label(cblabel,size=cbtitle_font)
    cb.ax.xaxis.set_label_position('top');
    cb.ax.tick_params(labelcolor='0.2', labelsize=tick_font, color='0.2',length=2.5, width=0.4, direction='out'); #change color of ticks?
    cb.ax.set_xticklabels([r'$\leq$-5','-4','-3','-2','-1','0','1','2','3','4',r'5$\leq$'])
    cb.outline.set_edgecolor('0.2')
    cb.outline.set_linewidth(0.4)
    
    #arrows
    plt.text(x0_redlab, y0_redlab, redlabel, size=arrow_font, ha='center', va='center')
    plt.text(x0_bluelab, y0_bluelab, bluelabel, size=arrow_font, ha='center', va='center')
    
    plt.arrow(x0_bluearr, y0_bluearr, xlen_bluearr, ylen_bluearr, width=arrow_width, linewidth=arrow_linew,\
               shape='left', head_width=arrow_headwidth, head_length=arrow_headlength,\
               facecolor=cmap_40, edgecolor='k', clip_on=False)
    plt.arrow(x0_redarr, y0_redarr, xlen_redarr, ylen_redarr, width=arrow_width, linewidth=arrow_linew,\
               shape='right', head_width=arrow_headwidth, head_length=arrow_headlength,\
               facecolor=cmap40, edgecolor='k', clip_on=False)
    
    #==============================================================================
    # watertemp jja scaled plotting
    #==============================================================================
    
    # continent fill color
    col_cont='white'
    
    # ocean fill color
    col_ocean='whitesmoke'
    
    # zero change color
    col_zero='gray'
    
    #========== COLORBAR ==========#
    
    # colorbar colormap setup
    cmap_whole = plt.cm.get_cmap('PRGn')
    cmap55 = cmap_whole(0.01)   
    cmap50 = cmap_whole(0.05)   #purple
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
    cmap_50 = cmap_whole(0.95)  #green
    cmap_55 = cmap_whole(0.99)
    
    # list of discrete colors as colormap for colorbar
    cmap = mpl.colors.ListedColormap([cmap_50,cmap_40,cmap_30,cmap_20,
                                      cmap20,cmap30,cmap40,cmap50], N=8)  
    
    # set color of over/under arrows on colorbar
    cmap.set_over(cmap55)
    cmap.set_under(cmap_55)
    
    # colorbar args
    values = [0,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
    tick_locs = np.arange(0,2.25,0.25)
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    # colorbar label
    cblabel = 'Δ lake temperature / Δ global mean air temperature (°C/°C)'
                 
    # bbox (arrow plot relative to this axis)
    cb_x0 = 0.553
    cb_y0 = 0.675
    cb_xlen = 0.4
    cb_ylen = 0.015
    
    #========== ARROWS ==========#
    
    # green arrow label
    greenlabel = 'Low sensitivity'
    x0_greenlab = 0.25
    y0_greenlab = -2.8
    
    # green arrow
    x0_greenarr = 0.495
    y0_greenarr = -3.3
    xlen_greenarr = -0.5
    ylen_greenarr = 0
    
    # purple arrow label
    purplelabel = 'High sensitivity'
    x0_purplelab = 0.75
    y0_purplelab = -2.8
    
    # purple arrow
    x0_purplearr = 0.505
    y0_purplearr = -3.3
    xlen_purplearr = 0.5
    ylen_purplearr = 0
    
    # general
    arrow_width = 0.25
    arrow_linew = 0.1
    arrow_headwidth = 0.5
    arrow_headlength = 0.06
    
    #========== PLOTTING ==========#
    
    lettercount += 1
    
    m = Basemap(llcrnrlon=-170, llcrnrlat=-60, urcrnrlon=180, urcrnrlat=90, suppress_ticks=False);
    m.ax = ax2
    m.drawcoastlines(linewidth=0.1);
    m.drawmapboundary(fill_color=col_ocean)
    m.fillcontinents(color=col_cont);
    m.pcolormesh(lon,lat,wt_jja_scaled,latlon=True,cmap=cmap,norm=norm,vmax=2,vmin=0,zorder=3)
    ax2.set_title(letters[lettercount-1],loc='left',pad=10,fontsize=title_font,fontweight='bold')
    ax2.set_title('JJA',loc='center',pad=10,fontsize=title_font)
    
    ax2.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                     bottom=False, top=False, left=False, right=False, color='0.2',\
                     labelcolor='0.2',width=0.4,direction="in",length=2.5)
    ax2.spines['bottom'].set_color('0.2')
    ax2.spines['bottom'].set_linewidth(0.4)
    ax2.spines['top'].set_color('0.2')
    ax2.spines['top'].set_linewidth(0.4)
    ax2.xaxis.label.set_color('0.2')
    ax2.spines['left'].set_color('0.2')
    ax2.spines['left'].set_linewidth(0.4)
    ax2.spines['right'].set_color('0.2')
    ax2.spines['right'].set_linewidth(0.4)
    ax2.yaxis.label.set_color('0.2')
    
    # colorbar setup
    cbax2 = f.add_axes([cb_x0, cb_y0, cb_xlen, cb_ylen])
    cb = mpl.colorbar.ColorbarBase(ax=cbax2, cmap=cmap,
                                   norm=norm,
                                   spacing='proportional',
                                   orientation='horizontal',
                                   extend='both',
                                   ticks=tick_locs)
    cb.set_label(cblabel,size=cbtitle_font)
    cb.ax.xaxis.set_label_position('top');
    cb.ax.tick_params(labelcolor='0.2', labelsize=tick_font, color='0.2',length=2.5, width=0.4, direction='out'); #change color of ticks?
    cb.ax.set_xticklabels([r'$\leq$0','0.25','0.5','0.75','1','1.25','1.5','1.75',r'2$\leq$'])
    cb.outline.set_edgecolor('0.2')
    cb.outline.set_linewidth(0.4)
    
    #arrows
    plt.text(x0_purplelab, y0_purplelab, purplelabel, size=arrow_font, ha='center', va='center')
    plt.text(x0_greenlab, y0_greenlab, greenlabel, size=arrow_font, ha='center', va='center')
    
    plt.arrow(x0_greenarr, y0_greenarr, xlen_greenarr, ylen_greenarr, width=arrow_width, linewidth=arrow_linew,\
               shape='left', head_width=arrow_headwidth, head_length=arrow_headlength,\
               facecolor=cmap_40, edgecolor='k', clip_on=False)
    plt.arrow(x0_purplearr, y0_purplearr, xlen_purplearr, ylen_purplearr, width=arrow_width, linewidth=arrow_linew,\
               shape='right', head_width=arrow_headwidth, head_length=arrow_headlength,\
               facecolor=cmap40, edgecolor='k', clip_on=False)
    
    #==============================================================================
    # ice index plotting
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
    
    # declare list of colors for discrete colormap of colorbar
    cmap = mpl.colors.ListedColormap([cmap_45,cmap_40,cmap_35,cmap_30,cmap_25,cmap_10,cmap0,
                                      cmap5,cmap10,cmap20,cmap30,cmap35,cmap40],N=13)
    
    # set color of over/under arrows in colorbar
    cmap.set_over(cmap55)
    cmap.set_under(cmap_55)
    
    # colorbar args
    values = [-90,-75,-60,-45,-30,-15,-5,5,15,30,45,60,75,90]
    tick_locs = [-90,-75,-60,-45,-30,-15,0,15,30,45,60,75,90]
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    # bbox (arrow plot relative to this axis)
    cb_x0 = 0.235
    cb_y0 = 0.175
    cb_xlen = 0.55
    cb_ylen = 0.015
    
    # colorbar label
    cblabel = 'Δ ice index (days)'
    
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
    
    #========== PLOTTING ==========#
    
    ax_list = [ax3,ax4,ax5]
    
    count = 0
    
    for array,ax in zip(iisig_data,ax_list):
        
        count += 1
        
        lettercount += 1
    
        m = Basemap(projection='npaeqd',round=True,boundinglat=20,\
                    lat_0=80,lon_0=0,resolution='l');
        m.ax = ax
        m.drawcoastlines(linewidth=0.05);
        m.drawmapboundary(linewidth=0.15,fill_color=col_ocean);
        m.fillcontinents(color=col_cont);
        m.pcolormesh(lon,lat,array,latlon=True,cmap=cmap,norm=norm,vmax=90,vmin=-90,zorder=3)
        ax.set_title(letters[lettercount-1],loc='left',fontsize=title_font,fontweight='bold')
        if count<=3:
            ax.set_title(ice_titles[count-1],loc='center',pad=10,fontsize=title_font)
        
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
                  length=2.5,width=0.35,direction='out'); 
    cb.ax.set_xticklabels([r'$\leq$-90','-75','-60','-45','-30','-15',\
                   '0','15','30','45','60','75',r'90$\leq$'])
    cb.outline.set_edgecolor('0.2')
    cb.outline.set_linewidth(0.4)
    
    # arrow setup
    plt.text(x0_bluelab, y0_bluelab, bluelabel, size=arrow_font, ha='center', va='center')
    plt.text(x0_redlab, y0_redlab, redlabel, size=arrow_font, ha='center', va='center')
    plt.arrow(x0_bluearr, y0_bluearr, xlen_bluearr, ylen_bluearr, width=arrow_width, linewidth=arrow_linew,\
              shape='right', head_width=arrow_headwidth, head_length=arrow_headlength,\
              facecolor=cmap40, edgecolor='k', clip_on=False)
    plt.arrow(x0_redarr, y0_redarr, xlen_redarr, ylen_redarr, width=arrow_width, linewidth=arrow_linew,\
              shape='left', head_width=arrow_headwidth, head_length=arrow_headlength,\
              facecolor=cmap_40, edgecolor='k', clip_on=False)
    
    plt.show()
    
    # save figure
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        f.savefig(outDIR+'/f3.png',bbox_inches='tight',dpi=500)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
