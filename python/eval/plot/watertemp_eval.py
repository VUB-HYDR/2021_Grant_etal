#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:18:37 2019

@author: Luke
"""

#==============================================================================
#SUMMARY
#==============================================================================


# This script plots sim-obs bias for icedur in isimip lake models vs era5-land 
# icedur


#==============================================================================
#IMPORT
#==============================================================================


import xarray as xr
import os
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#==============================================================================
#FUNCTIONS
#==============================================================================


# reading netcdfs
def reader(filename):
    ds = xr.open_dataset(filename, decode_times=False)
    if 'era5-land' in filename:
        da = ds.lmlt
    if not 'era5-land' in filename:
        da = ds.watertemp
    dims = list(da.dims)
    discard = []
    if 'time' in dims:
        discard.append('time')
    if 'levlak' in dims:
        discard.append('levlak')
    return da.squeeze(dim=discard,drop=True)

# ensemble math
def ensembler(arrays):
    concat_dim = np.arange(len(arrays))
    aligned = xr.concat(arrays,dim=concat_dim)
    mean = aligned.mean(dim='concat_dim')
    return mean

# processing function
def wteval_proc(inDIR,obsDIR):
    
    #==============================================================================
    #INITIALIZE DIRECTORIES
    #==============================================================================
    
    era5_file = 'era5-land_lakes_mixlayertemp_annual_timmean_1981_2018.nc'
    
    # set directory
    os.chdir(inDIR)
    
    #==============================================================================
    #DEFINE LOAD SETTINGS
    #==============================================================================
                    
    flag_var=0;  # 0: watertemp
                 # 1: lakeicefrac
                 # 2: icethick
                 # 3: icestart
                 # 4: iceend
                 # 5: icedur
     
    flag_prod=4; # 0: fldmean
                 # 1: sig
                 # 2: scaled_sig
                 # 3: eval
                 # 4: timmean
    
    # list of models (filename format)
    models = ['clm45','albm','lake','simstrat-uog','vic-lake']
    
    # list of variables (1st root variables, then processed)
    variables = ['watertemp','lakeicefrac','icethick','icestart','iceend','icedur']
    
    # list of products
    products = ['fldmean','sig','scaled_sig','eval','timmean']
    
    # list of rcps (no flagging)
    rcps = ['rcp60','rcp85']
    
    # assert settings
    var = variables[flag_var]
    prod = products[flag_prod]
    
    #==============================================================================
    #ACCESS FILES
    #==============================================================================
    
    files = []
    for rcp in rcps:
        for mod in models:
            for file in [file for file in sorted(os.listdir(inDIR))\
                     if mod in file and rcp in file and var in file and prod in file]:
                files.append(file)
                
    # filter model names iterable for existing model data for evaluation
    available_models = []
    for mod in models:
        count = 0
        for file in files:
            if mod in file:
                count = count+1
                if count == 1:
                    available_models.append(mod)
    
    #==============================================================================
    #OPEN MODEL TIMMEANS
    #==============================================================================
    
    # arrange + read files  
    arrays = {}
    for mod in available_models:
        arrays[mod] = []
        for file in files:
            if mod in file:
                arrays[mod].append(reader(file))
            
    os.chdir(obsDIR)
    era5 = reader(era5_file)
    
    lat = arrays['clm45'][0].lsmlat.values
    lon = arrays['clm45'][0].lsmlon.values
        
    isimip = {}
    for mod in available_models:
        ensemble = ensembler(arrays[mod])
        array = ensemble.values - era5.values
        isimip[mod] = array
        
    return available_models,\
           era5,isimip,lat,lon
    
# processing function
def wteval_plot(outDIR,flag_svplt,\
                available_models,\
                era5,isimip,lat,lon):
    
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
    
    # list of figure panel IDs
    letters = ['a','b','c','d','e',\
               'f','g','h','i','j']
    
    # list of model titles
    model_titles = ['CLM4.5', 'SIMSTRAT-UoG', 'ALBM']
    
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
    
    # list of discrete colors as colormap for colorbar
    cmap = mpl.colors.ListedColormap([cmap_50,cmap_45,cmap_40,cmap_35,cmap_30,cmap_25,cmap_20,cmap0,
                                      cmap20,cmap25,cmap30,cmap35,cmap40,cmap45,cmap50], N=15)  
    
    # set color of over/under arrows on colorbar
    cmap.set_over(cmap55)
    cmap.set_under(cmap_55)
    
    # colorbar args
    values = [-7,-6,-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5,6,7]
    tick_locs = [-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7]
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    # bbox
    cb_x0 = 0.0395
    cb_y0 = 0.225
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
    
    f = plt.figure(figsize=(10,15))
    
    # maps rect for left column, rect=[left, bottom, right, top]
    m_left = 0
    m_bottom = 0.25
    m_right = 0.55
    m_top = 0.65
    m_rect = [m_left, m_bottom, m_right, m_top]
    
    gs1 = gridspec.GridSpec(len(available_models),1)
    ax1 = f.add_subplot(gs1[0])    
    ax2 = f.add_subplot(gs1[1])
    ax3 = f.add_subplot(gs1[2])

    gs1.tight_layout(figure=f, rect=m_rect, h_pad=1)
    
    # maps axes
    maps_axs = [ax1,ax2,ax3]
    
    # histogram rect for right column, rect=[left, bottom, right, top]
    h_left = 0.6
    h_bottom = 0.2655
    h_right = 1.0
    h_top = 0.65
    h_rect = [h_left, h_bottom, h_right, h_top]
    
    gs2 = gridspec.GridSpec(len(available_models),1)
    ax5 = f.add_subplot(gs2[0])    
    ax6 = f.add_subplot(gs2[1])
    ax7 = f.add_subplot(gs2[2])

    gs2.tight_layout(figure=f, rect=h_rect, h_pad=3.225)
    
    # histogram axes
    hist_axs = [ax5,ax6,ax7]
    
    
    # maps
    
    count = 0
    
    for mod,ax in zip(isimip,maps_axs):
        
        count += 1
        
        m = Basemap(llcrnrlon=-170, llcrnrlat=-60, urcrnrlon=180, urcrnrlat=90, suppress_ticks=False);
        m.ax = ax
        m.drawcoastlines(linewidth=0.05);
        m.drawmapboundary(linewidth=0.15,fill_color=col_ocean);
        m.fillcontinents(color=col_cont);
        m.pcolormesh(lon,lat,isimip[mod],latlon=True,cmap=cmap,norm=norm,vmax=7,vmin=-7,zorder=3)
        
        ax.set_ylabel(model_titles[count-1],fontsize=title_font,labelpad=15)
        ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                         bottom=False, top=False, left=False, right=False, color='0.2',\
                         labelcolor='0.2',width=0.4,direction="in",length=2.5)
        ax.set_title(letters[count-1],loc='left',fontsize=title_font,fontweight='bold')
        ax.spines['bottom'].set_color('0.2')
        ax.spines['bottom'].set_linewidth(0.4)
        ax.spines['top'].set_color('0.2')
        ax.spines['top'].set_linewidth(0.4)
        ax.xaxis.label.set_color('0.2')
        ax.spines['left'].set_color('0.2')
        ax.spines['left'].set_linewidth(0.4)
        ax.spines['right'].set_color('0.2')
        ax.spines['right'].set_linewidth(0.4)
        ax.yaxis.label.set_color('0.2')
        
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
    cb.ax.tick_params(labelcolor='0.2', labelsize=tick_font, color='0.2',length=2.5, width=0.4, direction='out'); #change color of ticks?
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
    
    # histograms    
    for mod,ax in zip(isimip,hist_axs):
        
        count += 1
        
        ax.hist(isimip[mod].flatten(),bins=40,range=(-3,8),rwidth=0.9)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylabel('Grid cells [-]',fontsize=title_font)
        ax.set_title(letters[count-1],loc='left',fontsize=title_font,fontweight='bold')
        
    ax7.set_xlabel('Bias (°C)',fontsize=title_font)
    
    plt.show()
    
    # save figure
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        f.savefig(outDIR+'/si_f26.png',bbox_inches='tight',dpi=500)