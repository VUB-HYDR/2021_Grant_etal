#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:18:37 2019

@author: Luke
"""


#==============================================================================
# SUMMARY
#==============================================================================


# 09 March 2020

# This script plots 100 year scaled signal maps for ice indexes in isimip lake models 

# add snippets for new models under:
#       multimodel call (line 145)
#       model titles call (line ~185)
#       model panel lettering (line ~188)
#       figsize vertical lengthening (line ~292)
#       set_title call in plotting (line ~317)


#==============================================================================
# IMPORT
#==============================================================================


import xarray as xr
import os
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


#==============================================================================
# FUNCTIONS
#==============================================================================


# reading netcdfs
def reader(filename):
    ds = xr.open_dataset(filename, decode_times=False)
    ds = ds.apply(np.fabs) # absolute vals
    if 'icestart' in filename:
        da = ds.icestart
    if 'iceend' in filename:
        da = ds.iceend
    if 'icedur' in filename:
        da = ds.icedur
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
    ds = ds.apply(np.fabs)
    ds['lon'] = clm_lon
    ds['lat'] = clm_lat
    if 'icestart' in filename:
        da = ds.icestart
    if 'iceend' in filename:
        da = ds.iceend
    if 'icedur' in filename:
        da = ds.icedur
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
    mean = aligned.mean(dim='concat_dim').values
    return mean 

# proc/plotting
def iiscalesig_proc(inDIR):

    # set directory
    os.chdir(inDIR)

    #==============================================================================
    #DEFINE LOAD SETTINGS
    #==============================================================================
    
    flag_prod=2; # 0: fldmean
                 # 1: sig
                 # 2: scaled_sig
                 # 3: eval        
    
    # list of models (filename format)
    models = ['clm45','albm','simstrat-uog','vic-lake','lakemod'] # temporarily removed "lake"
    
    # list of products
    products = ['fldmean','sig','scaled_sig','eval','timmean']
    
    # list of end variables
    endvariables = ['icestart', 'iceend', 'icedur']
    
    # list of rcps (no flagging)
    rcps = ['rcp26', 'rcp60', 'rcp85']
    
    # assert settings
    prod = products[flag_prod]
    
    #==============================================================================
    #ACCESS FILES
    #==============================================================================
    
    files = []
    for rcp in rcps:
        for var in endvariables:
            for file in [file for file in sorted(os.listdir(inDIR))\
                     if rcp in file and var in file and prod in file]:
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
    for rcp in rcps:
        arrays[rcp] = {}
        for mod in available_models:
            arrays[rcp][mod] = {}
            for var in endvariables:
                for file in files:
                    count = 0
                    if mod != 'lakemod':
                        if rcp in file and mod in file and var in file:
                            arrays[rcp][mod][var] = reader(file)
                            count += 1
                            if count == 1:
                                if mod == 'clm45': # grabbing new coord dims for LAKE data-arrays
                                    clm_lat = arrays[rcp][mod][var].lat
                                    clm_lon = arrays[rcp][mod][var].lon
                    elif mod == 'lakemod':
                        if rcp in file and mod in file and var in file:
                            arrays[rcp][mod][var] = LAKEreader(file,clm_lon,clm_lat)
                        
    
    # add key for multi-model mean
    for rcp in rcps:    
        arrays[rcp]['mmm'] = {}
        for var in endvariables:
            temp_list = []
            for mod in available_models:
                temp_list.append(arrays[rcp][mod][var])
            arrays[rcp]['mmm'][var] = ensembler(temp_list)
    
    multimodel = {}
    for rcp in rcps:    
        multimodel[rcp] = [arrays[rcp]['mmm']['icestart'],arrays[rcp]['mmm']['iceend'],arrays[rcp]['mmm']['icedur'],\
                           arrays[rcp]['clm45']['icestart'],arrays[rcp]['clm45']['iceend'],arrays[rcp]['clm45']['icedur'],\
                           arrays[rcp]['simstrat-uog']['icestart'],arrays[rcp]['simstrat-uog']['iceend'],arrays[rcp]['simstrat-uog']['icedur'],\
                           arrays[rcp]['vic-lake']['icestart'],arrays[rcp]['vic-lake']['iceend'],arrays[rcp]['vic-lake']['icedur'],\
                           arrays[rcp]['albm']['icestart'],arrays[rcp]['albm']['iceend'],arrays[rcp]['albm']['icedur'],\
                           arrays[rcp]['lakemod']['icestart'],arrays[rcp]['lakemod']['iceend'],arrays[rcp]['lakemod']['icedur']]
        
    lat = arrays['rcp26']['clm45']['icestart'].lat.values
    lon = arrays['rcp26']['clm45']['icestart'].lon.values
    
    return multimodel,lat,lon
    
def iiscalesig_plot(multimodel,lat,lon,outDIR,flag_svplt,var,prod,dpi):
    
    # list of rcps (no flagging)
    rcps = ['rcp26', 'rcp60', 'rcp85']
    
    #==============================================================================
    #DEFINE PLOT SETTINGS
    #==============================================================================
    
    # font settings
    title_font = 13
    cbtitle_font = 13
    tick_font = 11
    arrow_font = 13
    
    # continent fill
    col_cont='lightgrey'
    
    # ocean fill
    col_ocean='white'
    
    # zero change color
    col_zero='white'
    
    # list of RCPs (no flagging; for file querying and figure panel IDs)
    RCPs = ['RCP 2.6','RCP 6.0','RCP 8.5']
    
    # list of ice index titles
    ice_titles = ['Ice onset', 'Ice break-up', 'Ice duration']
    
    # list of model titles
    model_titles = ['Multi-model mean', 'CLM4.5', 'SIMSTRAT-UoG', 'VIC-LAKE', 'ALBM', 'LAKE']
    
    # list of figure panel ids
    letters = ['a', 'b', 'c',\
               'd', 'e', 'f',\
               'g', 'h', 'i',\
               'j', 'k', 'l',\
               'm', 'n', 'o',\
               'p', 'q', 'r']
    
    #========== COLORBAR ==========#
    
    # identify colors
    cmap_whole = plt.cm.get_cmap('PRGn_r')
    cmap55 = cmap_whole(0.01)
    cmap50 = cmap_whole(0.05)   #green
    cmap45 = cmap_whole(0.1)
    cmap40 = cmap_whole(0.15)
    cmap35 = cmap_whole(0.2)
    cmap30 = cmap_whole(0.25)
    cmap25 = cmap_whole(0.3)
    cmap20 = cmap_whole(0.35)
    cmap10 = cmap_whole(0.4)
    cmap5 = cmap_whole(0.45)
    cmap0 = 'white'
    cmap_5 = cmap_whole(0.55)
    cmap_10 = cmap_whole(0.6)
    cmap_20 = cmap_whole(0.65)
    cmap_25 = cmap_whole(0.7)
    cmap_30 = cmap_whole(0.75)
    cmap_35 = cmap_whole(0.8)
    cmap_40 = cmap_whole(0.85)
    cmap_45 = cmap_whole(0.9)
    cmap_50 = cmap_whole(0.95)  #purple
    cmap_55 = cmap_whole(0.99)
    
    # declare list of colors for discrete colormap of colorbar
    cmap = mpl.colors.ListedColormap([cmap40,cmap_5,cmap_10,cmap_20,cmap_30,cmap_40,cmap_50],N=7)
    
    # set color of over/under arrows in colorbar
    cmap.set_over(cmap_55)
    
    # colorbar args
    values = [0,1,5,10,15,20,25,30]
    tick_locs = [0,1,5,10,15,20,25,30]
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    # bbox (arrow plot relative to this axis)
    cb_x0 = 0.2275
    cb_y0 = 0.22
    cb_xlen = 0.55
    cb_ylen = 0.01
    
    # colorbar label
    cblabel = 'Change in ice index scaled by global mean air temperature (days/°C)'
    
    #========== ARROWS ==========#
    
    # green arrow label
    greenlabel = 'High sensitivity'
    x0_greenlab = 0.525
    y0_greenlab = -2.7
    
    # green arrow
    x0_greenarr = 0.29
    y0_greenarr = -3.3
    xlen_greenarr = 0.42
    ylen_greenarr = 0
    
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
    x = 15
    y = 30
    
    #==============================================================================
    #INITIALIZE PLOTTING
    #==============================================================================
    
    lon, lat = np.meshgrid(lon, lat)
    
    for rcp in rcps:
        
        data = multimodel[rcp]
        
        f, axes = plt.subplots(int(len(data)/3),3,figsize=(x,y));
        
        count = 0
        
        for ax,array in zip(axes.flatten(),data):
            
            count = count+1
            
            m = Basemap(projection='npaeqd',round=True,boundinglat=20,\
                        lat_0=80,lon_0=0,resolution='l');
                
            m.ax = ax
            m.drawcoastlines(linewidth=0.05);
            m.drawmapboundary(linewidth=0.15,fill_color=col_ocean);
            m.fillcontinents(color=col_cont);
            m.pcolormesh(lon,lat,array,latlon=True,cmap=cmap,norm=norm,vmax=30,vmin=0,zorder=3)
            ax.set_title(letters[count-1],loc='left',fontsize=title_font,fontweight='bold')
            if count<=3:
                ax.set_title(ice_titles[count-1],loc='center',pad=10,fontsize=title_font)
            if count == 1:
                ax.set_ylabel(model_titles[0],fontsize=title_font,labelpad=15)
            if count == 4:
                ax.set_ylabel(model_titles[1],fontsize=title_font,labelpad=15)
            if count == 7:
                ax.set_ylabel(model_titles[2],fontsize=title_font,labelpad=15)
            if count == 10:
                ax.set_ylabel(model_titles[3],fontsize=title_font,labelpad=15)
            if count == 13:
                ax.set_ylabel(model_titles[4],fontsize=title_font,labelpad=15)
            if count == 16:
                ax.set_ylabel(model_titles[5],fontsize=title_font,labelpad=15)
                
    
        values = [0,1,5,10,15,20,25,30]
        tick_locs = [0,1,5,10,15,20,25,30]
        norm = mpl.colors.BoundaryNorm(values,cmap.N)
        cbax = f.add_axes([cb_x0, cb_y0, cb_xlen, cb_ylen])
        cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,
                                       norm=norm,
                                       spacing='uniform',
                                       orientation='horizontal',
                                       extend='max',
                                       ticks=tick_locs)
        cb.set_label(cblabel,size=cbtitle_font)
        cb.ax.xaxis.set_label_position('top');
        cb.ax.tick_params(labelcolor='0.2',labelsize=tick_font,color='0.2',\
                          length=2.5,width=0.35,direction='out'); 
        cb.ax.set_xticklabels(['0','1','5','10','15','20','25',r'30$\leq$'])
        cb.outline.set_edgecolor('0.2')
        cb.outline.set_linewidth(0.4)
        
        
        # arrow setup
        plt.text(x0_greenlab, y0_greenlab, greenlabel, size=arrow_font, ha='center', va='center')
        plt.arrow(x0_greenarr, y0_greenarr, xlen_greenarr, ylen_greenarr, width=arrow_width, linewidth=arrow_linew,\
                  shape='right', head_width=arrow_headwidth, head_length=arrow_headlength,\
                  facecolor=cmap_40, edgecolor='k', clip_on=False)
        
        plt.subplots_adjust(left=left_border, right=right_border, bottom=bottom_border, top=top_border,\
                            wspace=width_space, hspace=height_space)
        plt.show()
    
        # save figure
        if flag_svplt == 0:
            None
        elif flag_svplt == 1:
            if rcp == 'rcp26':
                f.savefig(outDIR+'/si_f14.png',bbox_inches='tight',dpi=dpi)
            elif rcp == 'rcp60':
                f.savefig(outDIR+'/si_f15.png',bbox_inches='tight',dpi=dpi)
            elif rcp == 'rcp85':
                f.savefig(outDIR+'/si_f16.png',bbox_inches='tight',dpi=dpi)
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
