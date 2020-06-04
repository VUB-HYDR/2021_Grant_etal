#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:18:37 2019

@author: Luke
"""

#==============================================================================
#SUMMARY
#==============================================================================


# 17 March 2020

# This script plots sim-obs bias for icedur in isimip lake models vs era5-land 
# icedur

# produces two histograms:
#   one per panel percentile bias and one per model percentil bias


#==============================================================================
#IMPORT
#==============================================================================


import xarray as xr
import os
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt



#==============================================================================
#FUNCTIONS
#==============================================================================


# reading netcdfs
def reader(filename):
    ds = xr.open_dataset(filename, decode_times=False)
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
    mean = aligned.mean(dim='concat_dim')
    return mean

# processing function
def iieval_proc(inDIR,obsDIR):
    
    #==============================================================================
    #INITIALIZE DIRECTORIES
    #==============================================================================
    
    era5_icestart = 'era5-land_icestart_timmean_1981_2018.nc'
    era5_iceend = 'era5-land_iceend_timmean_1981_2018.nc'
    era5_icedur = 'era5-land_icedur_timmean_1981_2018.nc'
    
    # set directory
    os.chdir(inDIR)
    
    #==============================================================================
    #DEFINE LOAD SETTINGS
    #==============================================================================
    
    # settings          
    flag_prod=4; # 0: fldmean
                 # 1: sig
                 # 2: scaled_sig
                 # 3: eval
                 # 4: timmean        
    
    # list of models (filename format)
    models = ['clm45','albm','simstrat-uog','vic-lake','lakemod'] # temporarily removed "lake"
    
    # list of products
    products = ['fldmean','sig','scaled_sig','eval','timmean']
    
    # list of end variables
    endvariables = ['icestart', 'iceend', 'icedur']
    
    # list of rcps (no flagging)
    rcps = ['rcp60','rcp85']
    
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
    for mod in available_models:
        arrays[mod] = {}
        for var in endvariables:
            arrays[mod][var] = []
            for file in files:
                if mod in file and var in file:
                    arrays[mod][var].append(reader(file))
                    
    # read reference prod     
    os.chdir(obsDIR)      
    era5 = {}
    era5['icestart'] = reader(era5_icestart)
    era5['iceend'] = reader(era5_iceend)
    era5['icedur'] = reader(era5_icedur)
    
    # get dimensions
    lat = arrays['clm45']['icedur'][0].lat.values
    lon = arrays['clm45']['icedur'][0].lon.values
    
    # ensemble work on same model+var and bias   
    isimip = {}
    for mod in available_models:
        isimip[mod] = {}
        for var in endvariables:
            ensemble = ensembler(arrays[mod][var])
            array = ensemble.values - era5[var].values
            isimip[mod][var] = array
    
    # pop keyed arrays into list for multi-model eval section
    multimodel = []
    for mod in available_models:
        for var in endvariables:
            multimodel.append(isimip[mod].pop(var))
            
    return available_models,\
           era5, multimodel,lat,lon

#
def iieval_plot(outDIR,flag_svplt,\
                available_models,\
                era5,multimodel,lat,lon):
    
    #==============================================================================
    #DEFINE PLOT SETTINGS
    #==============================================================================
            
    # list of figure panel ids
    letters = ['a', 'b', 'c',\
               'd', 'e', 'f',\
               'g', 'h', 'i',\
               'j', 'k', 'l',\
               'm', 'n', 'o']
    
    # font settings
    title_font = 13
    cbtitle_font = 13
    tick_font = 11
    arrow_font = 13
    
    # continent fill
    col_cont='white'
    
    # ocean fill
    col_ocean='whitesmoke'
    
    # zero change color
    col_zero='gray'
    
    # list of ice index titles
    ice_titles = ['Ice onset', 'Ice break-up', 'Ice duration']
    
    # list of model titles
    model_titles = ['CLM4.5', 'ALBM', 'SIMSTRAT-UoG', 'VIC-LAKE', 'LAKE']
    
    
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
    cmap = mpl.colors.ListedColormap([cmap_50,cmap_40,cmap_35,cmap_30,cmap_20,cmap_5,cmap0,
                                      cmap5,cmap20,cmap30,cmap35,cmap40,cmap50],N=13)
    
    # set color of over/under arrows in colorbar
    cmap.set_over(cmap55)
    cmap.set_under(cmap_55)
    
    # colorbar args
    values = [-30,-25,-20,-15,-10,-5,5,10,15,20,25,30]
    tick_locs = [-30,-25,-20,-15,-10,0,10,15,20,25,30]
    norm = mpl.colors.BoundaryNorm(values,cmap.N)
    
    # bbox (arrow plot relative to this axis)
    cb_x0 = 0.2275
    cb_y0 = 0.22
    cb_xlen = 0.55
    cb_ylen = 0.01
    
    # colorbar label
    cblabel = 'Bias in ice onset, break-up or duration (days)'
    
    #========== ARROWS ==========#
    
    # blue arrow label
    bluelabel = 'Later date (panels a,b) or longer duration (panel c)'
    x0_bluelab = 0.825
    y0_bluelab = -2.7
    
    # blue arrow
    x0_bluearr = 0.505
    y0_bluearr = -3.3
    xlen_bluearr = 0.7
    ylen_bluearr = 0
    
    # red arrow label
    redlabel = 'Earlier date (panels a,b) or shorter duration (panel c)'
    x0_redlab = 0.175
    y0_redlab = -2.7
    
    # red arrow
    x0_redarr = 0.495
    y0_redarr = -3.3
    xlen_redarr = -0.7
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
    
    #==============================================================================
    #PLOT TEST HIST
    #==============================================================================    
        
    # second; plot hist with 5th to 95th % on xaxis per model
    sets = {}
    sets['clm45'] = multimodel[0:3]
    sets['albm'] = multimodel[3:6]
    sets['simstrat-uog'] = multimodel[6:9]
    sets['vic-lake'] =  multimodel[9:11]
    sets['lakemod'] = multimodel[12:14]
    
    mins = {}
    maxs = {}
    for mod in available_models:
        mins[mod] = []
        maxs[mod] = []
        for array in sets[mod]:
            data = array.flatten()
            data = sorted(data[np.logical_not(np.isnan(data))])
            pcts = np.percentile(data,q=[5,95])
            data_min = pcts[0]
            data_max = pcts[1]
            mins[mod].append(data_min)
            maxs[mod].append(data_max)
    
    f,axes = plt.subplots(len(available_models),3,figsize=(12,4*len(available_models)))
    count=0
    
    for mod,ax in zip(multimodel,axes.flatten()):
        
        count += 1
        
        data = mod.flatten()
        data = sorted(data[np.logical_not(np.isnan(data))])
    
        if count <= 3:
            setmin = np.min(mins['clm45'])
            setmax = np.max(maxs['clm45'])
        
        if count >=4 and count <= 6:
            setmin = np.min(mins['simstrat-uog'])
            setmax = np.max(maxs['simstrat-uog'])
            
        if count >=7 and count <= 9:
            setmin = np.min(mins['vic-lake'])
            setmax = np.max(maxs['vic-lake'])
    
        if count >=10 and count <= 12:
            setmin = np.min(mins['lakemod'])
            setmax = np.max(maxs['lakemod'])
        
        ax.hist(data,range=(setmin,setmax),rwidth=0.75)
        ax.yaxis.set_label_position('right')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(letters[count-1],loc='left',pad=10,fontsize=title_font,fontweight='bold')
        
        if count <= 3:
            ax.set_title(ice_titles[count-1],loc='center',pad=10,fontsize=title_font)
        if count == 3:
            ax.set_ylabel(model_titles[0],fontsize=title_font,rotation=-90)
        if count == 6:
            ax.set_ylabel(model_titles[1],fontsize=title_font,rotation=-90)
        if count == 9:
            ax.set_ylabel(model_titles[2],fontsize=title_font,rotation=-90)
        if count == 12:
            ax.set_ylabel(model_titles[3],fontsize=title_font,rotation=-90)
        if count == 15:
            ax.set_ylabel(model_titles[4],fontsize=title_font,rotation=-90)
    
    f.text(0.5, -0.01, '5th to 95th percentile in bias per model (days)', ha='center', fontsize=title_font)
    f.text(-0.02, 0.5, 'Grid cells [-]', va='center', rotation='vertical', fontsize=title_font)
    
    plt.tight_layout()
    plt.show()  
    
    # save figure
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        f.savefig(outDIR+'/si_f28.png',bbox_inches='tight',dpi=500)
    
    #==============================================================================
    #PLOT MAIN FIG
    #==============================================================================
    
    lon, lat = np.meshgrid(lon, lat)
        
    f, axes = plt.subplots(int(len(multimodel)/3),3,figsize=(15,25));
    
    count = 0
    
    for mod,ax in zip(multimodel,axes.flatten()):
            
        count += 1
        
        m = Basemap(projection='npaeqd',round=True,boundinglat=20,\
                    lat_0=80,lon_0=0,resolution='l');
        m.ax = ax
        m.drawcoastlines(linewidth=0.05);
        m.drawmapboundary(linewidth=0.15,fill_color=col_ocean);
        m.fillcontinents(color=col_cont);
        m.pcolormesh(lon,lat,mod,latlon=True,cmap=cmap,norm=norm,vmax=30,vmin=-30,zorder=3)
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
    cb.ax.set_xticklabels([r'$\leq$-30','-25','-20','-15','-10',\
                           '0','10','15','20','25',r'30$\leq$'])
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
    
    plt.subplots_adjust(left=left_border, right=right_border, bottom=bottom_border, top=top_border,\
                        wspace=width_space, hspace=height_space)
    plt.show()
    
    # save figure
    if flag_svplt == 0:
        None
    elif flag_svplt == 1:
        f.savefig(outDIR+'/si_f27.png',bbox_inches='tight',dpi=500)
    
