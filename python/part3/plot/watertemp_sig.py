#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:18:37 2019

@author: Luke
"""


#==============================================================================
#SUMMARY
#==============================================================================


# 09 March 2020

# This script is used to plot historical-future signals in watertemp for
# 3 seasons across all available lake models

# converted whole plotting scheme to function for main_p1.py
# have commented out (because main script handles these):
#   directory calls
#   flag for svplt
#   flag for mmm or subplot


#==============================================================================
#IMPORT
#==============================================================================


import xarray as xr
import os
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl


#==============================================================================
#FUNCTIONS
#==============================================================================


# reading netcdfs
def reader(filename):
    ds = xr.open_dataset(filename, decode_times=False)
    da = ds.watertemp
    dims = list(da.dims)
    discard = []
    if 'time' in dims:
        discard.append('time')
    if 'levlak' in dims:
        discard.append('levlak')
    return da.squeeze(dim=discard,drop=True)

# reading netcdfs (adjusted for LAKE)
def LAKEreader(filename,clm_lon,clm_lat):
    ds = xr.open_dataset(filename, decode_times=False)
    ds['lon'] = clm_lon
    ds['lat'] = clm_lat
    da = ds.watertemp
    dims = list(da.dims)
    discard = []
    if 'time' in dims:
        discard.append('time')
    if 'levlak' in dims:
        discard.append('levlak')
    da.squeeze(dim=discard,drop=True)
    return da

# ensemble math
def ensembler(arrays):
    concat_dim = np.arange(len(arrays))
    aligned = xr.concat(arrays,dim=concat_dim)
    mean = aligned.mean(dim='concat_dim')
    return mean

# processing function
def wtsig_proc(inDIR):
    
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
     
    flag_prod=1; # 0: fldmean
                 # 1: sig
                 # 2: scaled_sig
                 # 3: eval
        
    # list of models (filename format)
    models = ['clm45','albm','simstrat-uog','vic-lake'] # temporarily removed "lake"
    
    # list of variables (1st root variables, then processed)
    variables = ['watertemp','lakeicefrac','icethick','icestart','iceend','icedur']
    
    # list of products
    products = ['fldmean','sig','scaled_sig','eval']
    
    # list of rcps (no flagging)
    rcps = ['rcp26','rcp60','rcp85']
    
    # list of seasons 
    seasons = ['DJF','MAM','JJA','SON']
    
    # assert settings
    var = variables[flag_var]
    prod = products[flag_prod]
    
    #==============================================================================
    #ACCESS FILES
    #==============================================================================
        
    files = []
    for mod in models:
        for seas in seasons:
            for rcp in rcps:
                for file in [file for file in sorted(os.listdir(inDIR))\
                             if mod+'_'+rcp+'_'+var+'_'+seas+'_'+prod in file]:
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
    #OPEN MODEL SIGNALS
    #==============================================================================
    
    # arrange + read files  
    arrays = {}
    for rcp in rcps:
        arrays[rcp] = {}
        for mod in available_models:
            arrays[rcp][mod] = {}
            for seas in seasons:
                for file in files:
                    count = 0
                    if mod != 'lakemod':
                        if rcp in file and mod in file and seas in file:
                            arrays[rcp][mod][seas] = reader(file)
                            count += 1
                            if count == 1:
                                if mod == 'clm45': # grabbing new coord dims for LAKE data-arrays
                                    clm_lat = arrays[rcp][mod][seas].lat
                                    clm_lon = arrays[rcp][mod][seas].lon
                    elif mod == 'lakemod':
                        if rcp in file and mod in file and seas in file:
                            arrays[rcp][mod][seas] = LAKEreader(file,clm_lon,clm_lat)
    
    # add key for multi-model mean
    for rcp in rcps:    
        arrays[rcp]['mmm'] = {}
        for seas in seasons:
            temp_list = []
            for mod in available_models:
                temp_list.append(arrays[rcp][mod][seas])
            arrays[rcp]['mmm'][seas] = ensembler(temp_list)
    
    multimodel = {}
    for rcp in rcps:    
        multimodel[rcp] = [arrays[rcp]['mmm']['DJF'],arrays[rcp]['mmm']['MAM'],arrays[rcp]['mmm']['JJA'],arrays[rcp]['mmm']['SON'],\
                           arrays[rcp]['clm45']['DJF'],arrays[rcp]['clm45']['MAM'],arrays[rcp]['clm45']['JJA'],arrays[rcp]['clm45']['SON'],\
                           arrays[rcp]['simstrat-uog']['DJF'],arrays[rcp]['simstrat-uog']['MAM'],arrays[rcp]['simstrat-uog']['JJA'],arrays[rcp]['simstrat-uog']['SON'],\
                           arrays[rcp]['albm']['DJF'],arrays[rcp]['albm']['MAM'],arrays[rcp]['albm']['JJA'],arrays[rcp]['albm']['SON'],\
                           arrays[rcp]['vic-lake']['DJF'],arrays[rcp]['vic-lake']['MAM'],arrays[rcp]['vic-lake']['JJA'],arrays[rcp]['vic-lake']['SON']]
        
    
    lat = arrays['rcp26']['clm45']['DJF'].lat.values
    lon = arrays['rcp26']['clm45']['DJF'].lon.values
    
    return multimodel,lat,lon,available_models

# plotting function    
def wtsig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi):
    
    # list of rcps (no flagging)
    rcps = ['rcp26','rcp60','rcp85']
    
    # list of seasons 
    seasons = ['DJF','MAM','JJA','SON']
    
    #==============================================================================
    #DEFINE PLOT SETTINGS
    #==============================================================================
    
    # font settings
    title_font = 15
    cbtitle_font = 15
    tick_font = 12
    arrow_font = 14
    
    # list of model titles
    model_titles = ['Multi-model mean', 'CLM4.5', 'SIMSTRAT-UoG', 'ALBM', 'VIC-LAKE']
    
    # continent fill color
    col_cont='white'
    
    # ocean fill color
    col_ocean='whitesmoke'
    
    # zero change color
    col_zero='gray'
    
    # list of figure panel IDs
    letters = ['a', 'b', 'c', 'd',\
               'e', 'f', 'g', 'h',\
               'i', 'j', 'k', 'l',\
               'm', 'n', 'o', 'p',\
               'q', 'r', 's', 't']
    
    #========== COLORBAR ==========#
    
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
    cblabel = 'Change in lake temperature (at 2m depth) in °C'
    
    # bbox (arrow plot relative to this axis)
    cb_x0 = 0.275
    cb_y0 = 0.175
    cb_xlen = 0.45
    cb_ylen = 0.015
    
    #========== ARROWS ==========#
    
    # blue arrow label
    bluelabel = 'Colder'
    x0_bluelab = 0.25
    y0_bluelab = -2.8
    
    # blue arrow
    x0_bluearr = 0.495
    y0_bluearr = -3.4
    xlen_bluearr = -0.4
    ylen_bluearr = 0
    
    # red arrow label
    redlabel = 'Warmer'
    x0_redlab = 0.75
    y0_redlab = -2.8
    
    # red arrow
    x0_redarr = 0.505
    y0_redarr = -3.4
    xlen_redarr = 0.4
    ylen_redarr = 0
    
    # general
    arrow_width = 0.25
    arrow_linew = 0.1
    arrow_headwidth = 0.5
    arrow_headlength = 0.06
    
    #========== SUBPLOTS ==========#
    
    left_border = 0.15
    right_border = 0.85
    bottom_border = 0.25
    top_border = 0.925
    width_space = 0.05
    height_space = 0.05
    
    # figsize = (x,y)
    x = 25
    y = 20
    
    #==============================================================================
    #INITIALIZE PLOTTING
    #==============================================================================
    
    lon, lat = np.meshgrid(lon, lat)
    
    for rcp in rcps:
    
        data = multimodel[rcp]
        
        f, axes = plt.subplots(int(len(data)/4),4,figsize=(x,y));
        
        count = 0
        
        for ax,array in zip(axes.flatten(),data):
            
            count += 1
            
            m = Basemap(llcrnrlon=-170, llcrnrlat=-60, urcrnrlon=180, urcrnrlat=90, suppress_ticks=False);
            m.ax = ax
            m.drawcoastlines(linewidth=0.1);
            m.drawmapboundary(fill_color=col_ocean)
            m.fillcontinents(color=col_cont);
            m.pcolormesh(lon,lat,array,latlon=True,cmap=cmap,norm=norm,vmax=5,vmin=-5,zorder=3)
            ax.set_title(letters[count-1],loc='left',pad=10,fontsize=title_font,fontweight='bold')
            if count<=4:
                ax.set_title(seasons[count-1],loc='center',pad=10,fontsize=title_font)
            if count == 1:
                ax.set_ylabel(model_titles[0],fontsize=title_font,labelpad=15)
            if count == 5:
                ax.set_ylabel(model_titles[1],fontsize=title_font,labelpad=15)
            if count == 9:
                ax.set_ylabel(model_titles[2],fontsize=title_font,labelpad=15)
            if count == 13:
                ax.set_ylabel(model_titles[3],fontsize=title_font,labelpad=15)
            if count == 17:
                ax.set_ylabel(model_titles[4],fontsize=title_font,labelpad=15)
    
            ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False,
                             bottom=False, top=False, left=False, right=False, color='0.2',\
                             labelcolor='0.2',width=0.4,direction="in",length=2.5)
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
                
                
        #colorbar setup
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
        
        plt.subplots_adjust(left=left_border, right=right_border, bottom=bottom_border, top=top_border,\
                            wspace=width_space, hspace=height_space)
        plt.show()
    
        # save figure
        if flag_svplt == 0:
            None
        elif flag_svplt == 1:
            if rcp == 'rcp26':
                f.savefig(outDIR+'/si_f05.png',bbox_inches='tight',dpi=dpi)
            elif rcp == 'rcp60':
                f.savefig(outDIR+'/si_f06.png',bbox_inches='tight',dpi=dpi)
            elif rcp == 'rcp85':
                f.savefig(outDIR+'/si_f07.png',bbox_inches='tight',dpi=dpi)
            
