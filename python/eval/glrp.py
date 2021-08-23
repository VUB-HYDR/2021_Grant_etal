#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:03:59 2021

@author: luke
"""

# Package ID: knb-lter-ntl.10001.3 Cataloging System:https://pasta.lternet.edu.
# Data set title: Globally distributed lake surface water temperatures collected in situ and by 			satellites; 1985-2009.

# 
# This program creates numbered PANDA dataframes named dt1,dt2,dt3...,
# one for each data table in the dataset. It also provides some basic
# summaries of their contents. NumPy and Pandas modules need to be installed
# for the program to run. 



# =============================================================================
# import
# =============================================================================



import numpy as np
import pandas as pd 
import os
import xarray as xr
from scipy import stats as sts
import matplotlib.pyplot as plt
import seaborn as sb
import geopandas as gpd
from shapely.geometry import Polygon
from shapely import wkt
import os
import gdal
import copy as cp
from collections import OrderedDict
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.feature as cfeature 
cmaps = OrderedDict()
import pickle as pk



# =============================================================================
# functions
# =============================================================================



def era5l_reader(file,
                 var):
    
    da = xr.open_dataset(file,decode_times=False)
    da = da[var]
    time = pd.date_range(start='1981-01-01',
                         end='2006-01-01',
                         freq='YS')
    da['time'] = time
    da = da.rename({'longitude':'lon',
                    'latitude':'lat'})
    
    return da

def rasterize(feature_name,lon_min,lon_max,lat_min,lat_max,resolution,filename):
    """
    This function rasterizes a .shp file and saves it as a .tiff in the same directory
    Only for global extent
    input:      feature_name: Fieldname of shapefile to be burned in raster
                resolution: horizontal resolution in degrees  
                filename: input and output filename
    """
    # define command
    command = 'gdal_rasterize -a '+ feature_name\
    + ' -ot Float32 -of GTiff -te '+ str(lon_min)+' '+str(lat_min)+' '+str(lon_max)+' '+str(lat_max)+' -tr ' + str(resolution) +' '+ str(resolution)\
    + ' -co COMPRESS=DEFLATE -co PREDICTOR=1 -co ZLEVEL=6 -l '+ filename\
    + ' ' + filename+'.shp ' + filename +'.tiff'

    os.system(command)    

def read_raster(filename):
    """
    Function to read raster file
    input: file name of raster (ends in .tiff)
    output: 2D numpy array
    """
    raster = gdal.Open(filename)
    myarray = np.array(raster.GetRasterBand(1).ReadAsArray())
    myarray = np.flipud(myarray)

    return myarray

def slope_field(xarr):  
    
    # getting shapes
    m = np.prod(xarr.shape[1:]).squeeze()
    n = xarr.shape[0]
    
    # creating x and y variables for linear regression
    # x = xarr.time.to_pandas().index.to_julian_date().values[:, None]
    x = xarr.time.dt.year.values[:,None]
    y = xarr.to_masked_array().reshape(n, -1)
    
    # ############################ #
    # LINEAR REGRESSION DONE BELOW #
    xm = x.mean(0)  # mean
    ym = y.mean(0)  # mean
    ya = y - ym  # anomaly
    xa = x - xm  # anomaly
    
    # variance and covariances
    xss = (xa ** 2).sum(0) / (n - 1)  # variance of x (with df as n-1)
    yss = (ya ** 2).sum(0) / (n - 1)  # variance of y (with df as n-1)
    xys = (xa * ya).sum(0) / (n - 1)  # covariance (with df as n-1)
    # slope and intercept
    slope = xys / xss
    intercept = ym - (slope * xm)
    # statistics about fit
    df = n - 2
    r = xys / (xss * yss)**0.5
    t = r * (df / ((1 - r) * (1 + r)))**0.5
    p = sts.distributions.t.sf(abs(t), df)
    
    # preparing outputs
    out = xarr[:2].mean('time')
    # first create variable for slope and adjust meta
    xarr_slope = out.copy()
    xarr_slope.name = '_slope'
    xarr_slope.attrs['units'] = 'K / year'
    xarr_slope.values = slope.reshape(xarr.shape[1:])
    # do the same for the p value
    xarr_p = out.copy()
    xarr_p.name = '_Pvalue'
    xarr_p.attrs['info'] = "If p < 0.05 then the results from 'slope' are significant."
    xarr_p.values = p.reshape(xarr.shape[1:])
    # join these variables
    xarr_out = xarr_slope.to_dataset(name='slope')
    xarr_out['pval'] = xarr_p

    #return xarr_out
    return xarr_slope,xarr_p


def pixel(arr,
          lon,
          lat,
          out_arr = False):
    
    if out_arr == False:
        series = arr.sel(lon=lon,
                         lat=lat,
                         drop=True).squeeze().values.item()
    elif out_arr == True:
        series = arr.sel(lon=lon,
                         lat=lat,
                         drop=True).squeeze()
    
    return series


def df_indexer(slope_arr,
               series_arr,
               df,
               lon,
               lat):
    
    val = df.loc[(df['lat'] == lat) & (df['lon'] == lon),'obs'].item()
    
    latx = slope_arr.where(slope_arr == val,drop=True).squeeze().lat.values.item()
    lonx = slope_arr.where(slope_arr == val,drop=True).squeeze().lon.values.item()
    
    series = series_arr.sel(lat=latx,
                            lon=lonx,
                            drop=True).squeeze()
    
    series = series.interpolate_na(dim='time')
    
    return series
    
def arr_to_df(arr1,
              arr2):
    
    """ Take two arrays (matching ERA5L and obs). For each significant obs trend
    in arr1, take lat + lon coords, find value for this coord in ERA5L and append 
    arr1 value, arr2 value, lat and lon to dataframe.
    
    Parameters
    ----------
    arr1 : obs
    arr2 : ERA5L
    
    Returns
    ------- 
    Pandas dataframe
    """
    # fails because of d_coords for d yielding multiple locations for lat
    frame = {'obs':[],'era5l':[],'lat':[],'lon':[]}
    df = pd.DataFrame(data=frame)
    vals = arr1.values.flatten()
    data = vals[~np.isnan(vals)]
    data = np.unique(data[data != 0])
    for d in data: 
        d_coords = arr1.where(arr1==d,drop=True).squeeze()
        
        try:
            lat = np.around(d_coords.lat.values.item(),1)
        except:
            coord_len_lat = len(d_coords.lat.values)
        try:
            lon = np.around(d_coords.lon.values.item(),1)
        except:
            coord_len_lon = len(d_coords.lon.values)
            
        try:
            if coord_len_lat and coord_len_lon:
                for lo in lon:
                    for la in lat:
                        lo = np.around(lo.item(),1)
                        la = np.around(la.item(),1)
                        e = pixel(arr2,
                                  lo,
                                  la,
                                  out_arr=False)
                        df = df.append({'obs':d,'era5l':e,'lat':la,'lon':lo}, ignore_index=True)
        except:
            e = pixel(arr2,
                      lon,
                      lat,
                      out_arr=False)
            df = df.append({'obs':d,'era5l':e,'lat':lat,'lon':lon}, ignore_index=True)
        
    return df.dropna()


def ensembler(data):
    concat_dim = np.arange(len(data))
    aligned = xr.concat(data,dim=concat_dim)
    ens_mean = aligned.mean(dim='concat_dim')
    ens_std = aligned.std(dim='concat_dim')
    ens_max = aligned.max(dim='concat_dim')
    ens_min = aligned.min(dim='concat_dim')
    ens_roll = ens_mean.rolling(years=5, center=True).mean()
    dict_ens = {}
    dict_ens['mean'] = ens_mean
    dict_ens['std'] = ens_std
    dict_ens['max'] = ens_max
    dict_ens['min'] = ens_min
    dict_ens['roll'] = ens_roll
    return dict_ens


def plotter(time,
            ens_mean,
            ens_std,
            ens_max,
            ens_min,
            ens_roll,
            ax,
            lw_mean,
            lw_roll,
            col_mean,
            col_fill_a,
            col_fill_b,
            ub_alpha):
    
    ens_mean = ens_mean.values
    ens_std = ens_std.values
    ens_max = ens_max.values
    ens_min = ens_min.values
    ens_roll = ens_roll.values
    
    # plot mean line
    h = ax.plot(time, 
                ens_mean,
                lw=lw_mean, 
                color=col_mean,
                zorder=3)
    
    # plot rolling mean line
    h = ax.plot(time, 
                ens_roll,
                lw=lw_roll, 
                color=col_mean, 
                zorder=4)
            
    return h,ax


def tser_plotter(series_on,
                 series_off,
                 series_dur,
                 x,
                 y,
                 xmin,
                 xmax,
                 ymin,
                 ymax,
                 labels,
                 xticks,
                 xtick_labels,
                 tick_font,
                 title_font,
                 axis_font,
                 legend_font,
                 legend_entrylen,
                 legend_entrypad,
                 legendcols,
                 xlabel_xpos,
                 xlabel_ypos,
                 xlabel,
                 ylabel_xpos,
                 ylabel_ypos,
                 ylabel,
                 colors,
                 ub_alpha,
                 icevars,
                 x0,
                 y0,
                 xlen,
                 ylen):

    f, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(x,y))
    
    for s,c in zip(series_on,colors):
        
        time = s['mean'].years.values
        
        h,ax = plotter(time,
                       s['mean'],
                       s['std'],
                       s['max'],
                       s['min'],
                       s['roll'],
                       ax1,
                       lw_mean,
                       lw_roll,
                       c['mean'],
                       c['fill_a'],
                       c['fill_b'],
                       ub_alpha)
        
    for s,c in zip(series_off,colors):
        
        time = s['mean'].years.values
    
        h,ax = plotter(time,
                       s['mean'],
                       s['std'],
                       s['max'],
                       s['min'],
                       s['roll'],
                       ax2,
                       lw_mean,
                       lw_roll,
                       c['mean'],
                       c['fill_a'],
                       c['fill_b'],
                       ub_alpha)
        
    for s,c in zip(series_dur,colors):
        
        time = s['mean'].years.values
    
        h,ax = plotter(time,
                       s['mean'],
                       s['std'],
                       s['max'],
                       s['min'],
                       s['roll'],
                       ax3,
                       lw_mean,
                       lw_roll,
                       c['mean'],
                       c['fill_a'],
                       c['fill_b'],
                       ub_alpha)
        
    i = 0    
    for ax in (ax1,ax2,ax3):        
        ax.set_xlim(xmin,xmax)
        ax.autoscale(axis='y')
        ax.xaxis.set_ticks(xticks)
        ax.tick_params(labelsize=tick_font,axis="x",direction="in", left="off",labelleft="on")
        ax.tick_params(labelsize=tick_font,axis="y",direction="in")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_axisbelow(True)
        ax.set_ylabel(icevars[i])
        ax.set_title(letters[i],loc='left',fontsize=title_font,fontweight='bold')
        i += 1
        
    ax1.xaxis.set_ticklabels([]) 
    ax2.xaxis.set_ticklabels([]) 
    ax3.xaxis.set_ticklabels(xtick_labels)
    
    ax1.legend(handles,
              labels,
              bbox_to_anchor=(x0, y0, xlen, ylen), 
              loc=3,   #bbox: (x, y, width, height)
              ncol=3,
              fontsize=legend_font, 
              mode="expand", 
              borderaxespad=0.,\
              frameon=False, 
              columnspacing=0.05, 
              handlelength=legend_entrylen, 
              handletextpad=legend_entrypad)
    
    # labels
    f.text(xlabel_xpos, 
           xlabel_ypos, 
           xlabel, 
           ha='center', 
           fontsize=axis_font)
    
    f.savefig('si_f31_v3.png',bbox_inches='tight',dpi=200)
    
def map_plotter(x,
                y,
                proj,
                new_extent,
                lake_pts,
                col_glrp,
                x0, 
                y0, 
                xlen, 
                ylen):
    

    f, ax = plt.subplots(nrows=1,ncols=1,
                         figsize=(x,y),
                         subplot_kw=dict(projection=proj))
    ax.set_extent(new_extent,ccrs.PlateCarree())
    
    lake_pts.plot(ax=ax,
                  markersize=4,
                  color=col_glrp['mean'],
                  zorder=2,
                  transform=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.LAND, 
                   zorder=1, 
                   edgecolor='black',
                   linewidth=0.5)
    
    legend_handles = [Line2D([0], [0],
                             marker='o',
                             color='w',
                             label='GLRIP locations',
                             markerfacecolor=col_glrp['mean'])]
    
    ax.legend(handles=legend_handles,
              bbox_to_anchor=(x0, y0, xlen, ylen),
              frameon=False)
    
    f.savefig('glrp_locations_final_v2.png',bbox_inches='tight',dpi=200)


def c(x):
   col = plt.cm.Greys(x)
   fig, ax = plt.subplots(figsize=(1,1))
   fig.set_facecolor(col)
   ax.axis("off")
   plt.show()



# =============================================================================
# settings
# =============================================================================



title_font = 10
tick_font = 8
axis_font = 10
legend_font = 9

#========== LINE THICKNESS ==========#

# mean line thickness
lw_mean = 0.5
lw_roll = 2

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

col_glrp = {}
col_era = {}

col_glrp['mean'] = plt.cm.Blues(0.9)
col_glrp['fill_a'] = plt.cm.YlOrBr(0.7)
col_glrp['fill_b'] = plt.cm.YlOrBr(0.4)
col_era['mean'] = plt.cm.Greys(0.9)
col_era['fill_a'] = plt.cm.Greys(0.7)
col_era['fill_b'] = plt.cm.Greys(0.4)

# legend colors
legendcols = [col_era['mean'],
              col_glrp['mean']]

colors = [col_era,
          col_glrp]

ub_alpha = 0.5

#========== AXII ==========#

# figsize = (x,y)
x = 8
y = 8

# subplots_adjust
hspace = 0.5
top = 0.9

ymin = -15   # ymin
ymax = 15    # ymax
xmin = 1980 # xmin
xmax = 2006 # xmax

# x ticks/labels 
xticks = np.arange(1980,2010,5)
xtick_labels = [1980,None,1990,None,2000,None]

# x axis label
xlabel = 'Years'
xlabel_xpos = 0.5
xlabel_ypos = 0.05

# y axis label
ylabel = 'Ice index anomaly (days)'
ylabel_xpos = 0.05
ylabel_ypos = 0.535

# xaxis tick label sharing
axis_share = False

#========== LEGEND ==========#

# labels
lab_glrp = 'GLRIP'
lab_era = 'ERA5L'
labels = [lab_era,
          lab_glrp]

# bbox
x0 = 0.75
y0 = 1
xlen = 0.25
ylen = 0.9

# space between entries
legend_entrypad = 0.5

# length per entry
legend_entrylen = 0.75

icevars=['Ice onset anomaly',
         'Ice breakup anomaly',
         'Ice duration anomaly']

handles = [Line2D([0],[0],linestyle='-',lw=2,color=legendcols[0]),\
           Line2D([0],[0],linestyle='-',lw=2,color=legendcols[1])]
    
letters = ['a','b','c','d']

    
    
# =============================================================================
# retrieve data
# =============================================================================



os.chdir('/home/luke/documents/data/glrp/')

infile1  ="liag_ice.csv"
                 
dt1 = pd.read_csv(infile1,
                  skiprows=1,
                  sep=",",
                  names=["iceon_year",     
                         "iceon_month",     
                         "iceon_day",  
                         "iceoff_year",
                         "iceoff_month",
                         "iceoff_day",
                         "duration",
                         "froze",
                         "latitude",
                         "longitude",
                         "lakename",
                         "siteID",     
                         "value"],
                  index_col=False)

dt1 = dt1.replace(to_replace=-999,
                  value=np.nan).dropna()

# list of all lakes
lake_names = np.unique(dt1['lakename'].values)

# list of lakes with records >= 20 years
lr_lakes = []
for l in lake_names:
    record = dt1.loc[dt1['lakename']==l]
    yrs = len(record)
    if yrs >= 15:
        lr_lakes.append(record)
        
dt2 = pd.concat(lr_lakes)

lr_lakes = np.unique(dt2['lakename'].values)
nr_lr_lakes = len(lr_lakes)

day_cal_on = {}
day_cal_on[1] = 1+365
day_cal_on[2] = 32+365
day_cal_on[3] = 60+365
day_cal_on[4] = 91+365
day_cal_on[5] = 121+365
day_cal_on[6] = 152+365
day_cal_on[7] = 182+365
day_cal_on[8] = 213+365
day_cal_on[9] = 244+365
day_cal_on[10] = 274
day_cal_on[11] = 305
day_cal_on[12] = 335

day_cal_off = {}
day_cal_off[1] = 1+365
day_cal_off[2] = 32+365
day_cal_off[3] = 60+365
day_cal_off[4] = 91+365
day_cal_off[5] = 121+365
day_cal_off[6] = 152+365
day_cal_off[7] = 182+365
day_cal_off[8] = 213+365
day_cal_off[9] = 244
day_cal_off[10] = 274
day_cal_off[11] = 305
day_cal_off[12] = 335

dt2_wrk = dt2.copy(deep=True)
dt2_wrk['hyr_on_day'] = dt2_wrk.apply(lambda row: day_cal_on[int(row.iceon_month)] + row.iceon_day -1, axis=1)
dt2_wrk['hyr_off_day'] = dt2_wrk.apply(lambda row: day_cal_off[int(row.iceoff_month)] + row.iceoff_day -1, axis=1)

data = dt2_wrk.drop(columns=['iceon_month',
                             'iceon_day',
                             'iceoff_year',
                             'iceoff_month',
                             'iceoff_day',
                             'froze',
                             'lakename',
                             'value']).rename(columns={'iceon_year':'year'})

data_on = data.copy(deep=True).drop(columns=['hyr_off_day',
                                             'duration']).rename(columns={'hyr_on_day':'value'})
data_on['longitude'] = data_on.apply(lambda row: row.longitude + 180,axis=1)

data_off = data.copy(deep=True).drop(columns=['hyr_on_day',
                                              'duration']).rename(columns={'hyr_off_day':'value'})
data_off['longitude'] = data_off.apply(lambda row: row.longitude + 180,axis=1)

data_dur = data.copy(deep=True).drop(columns=['hyr_on_day',
                                              'hyr_off_day']).rename(columns={'duration':'value'})
data_dur['longitude'] = data_dur.apply(lambda row: row.longitude + 180,axis=1)



siteID_on = sorted(data_on['siteID'].unique())
siteID_off = sorted(data_off['siteID'].unique())
siteID_dur = sorted(data_dur['siteID'].unique())

# reanalysis read in
os.chdir("/home/luke/documents/data/glrp/")

era5l_on_file = "era5-land_lakes_icecover_start_1981_2006.nc"
era5l_off_file = "era5-land_lakes_icecover_end_1981_2006.nc"
era5l_dur_file = "era5-land_lakes_icecover_duration_1981_2006.nc"

var='icestart'
era5l_on = era5l_reader(era5l_on_file,
                        var)
var='iceend'
era5l_off = era5l_reader(era5l_off_file,
                         var)
var='iceduration'
era5l_dur = era5l_reader(era5l_dur_file,
                         var)

glrp_on_series = []
glrp_on_mean = []
era5l_on_series = []
era5l_on_mean = []

# =============================================================================
#     frame = {'obs':[],'era5l':[],'lat':[],'lon':[]}
#     df = pd.DataFrame(data=frame)
#     vals = arr1.values.flatten()
#     data = vals[~np.isnan(vals)]
#     data = np.unique(data[data != 0])
#     for d in data: 
#         d_coords = arr1.where(arr1==d,drop=True).squeeze()
#         
#         try:
#             lat = np.around(d_coords.lat.values.item(),1)
#         except:
#             coord_len_lat = len(d_coords.lat.values)
#         try:
#             lon = np.around(d_coords.lon.values.item(),1)
#         except:
#             coord_len_lon = len(d_coords.lon.values)
#             
#         try:
#             if coord_len_lat and coord_len_lon:
#                 for lo in lon:
#                     for la in lat:
#                         lo = np.around(lo.item(),1)
#                         la = np.around(la.item(),1)
#                         e = pixel(arr2,
#                                   lo,
#                                   la,
#                                   out_arr=False)
#                         df = df.append({'obs':d,'era5l':e,'lat':la,'lon':lo}, ignore_index=True)
#         except:
#             e = pixel(arr2,
#                       lon,
#                       lat,
#                       out_arr=False)
#             df = df.append({'obs':d,'era5l':e,'lat':lat,'lon':lon}, ignore_index=True)
#         
#     return df.dropna()
# =============================================================================

frame = {'siteID':[],'lat':[],'lon':[]}
lake_locs = pd.DataFrame(data=frame)

for s in siteID_on:
    
    # data per siteID
    data = data_on.loc[data_on['siteID']==s]
    lat = np.around(data['latitude'].unique(),1)
    lon = np.around(data['longitude'].unique(),1)
    years = sorted(data['year'].unique())
    s_d = []
    
    # handle duplicate entries (assuming these represent multiple onsets/breakups in a given year)
        # minimum for onset
        # maximum for breakup

    for y in years:
        
        d = data.loc[data['year']==y]
        
        if len(d) == 1:
            
            v = d['value'].item()
            
        elif len(d) > 1:
            
            v = d['value'].min()
        s_d.append(v)
        
    # glrp series
    years = np.asarray(years,
                       dtype=np.int64)
    s_da = xr.DataArray(s_d,
                        coords=[years],
                        dims=['years']).sel(years=slice(1981,2006))
    s_da_mean = s_da.mean(dim='years')
    s_da = s_da - s_da_mean
    
    # era5l series
    e_da_i = era5l_on.sel(lat=lat,
                          lon=lon,
                          method='nearest',
                          tolerance=0.1,
                          drop=True).squeeze()
    nan_check = e_da_i[~np.isnan(e_da_i)]
    
    if len(nan_check) > 0:
        
        e_da_v = e_da_i.values
        e_da_years = e_da_i.time.dt.year.values
        e_da_f = xr.DataArray(e_da_v,
                              coords=[e_da_years],
                              dims=['years'])
        e_da_f = e_da_f.where(e_da_f.years == s_da.years)
        e_da_f_mean = e_da_f.mean(dim='years')
        e_da_f = e_da_f - e_da_f_mean
        
        era5l_on_series.append(e_da_f)
        era5l_on_mean.append(e_da_f_mean.values.item())
        
        glrp_on_series.append(s_da)
        glrp_on_mean.append(s_da_mean.values.item())
        
        lake_locs = lake_locs.append({'siteID':s,'lat':lat.item(),'lon':lon.item()-180}, ignore_index=True)
        
    elif len(nan_check) == 0:
        
        pass

   
glrp_off_series = []
glrp_off_mean = []
era5l_off_series = []
era5l_off_mean = [] 

for s in siteID_off:
    
    # data per siteID
    data = data_off.loc[data_off['siteID']==s]
    lat = np.around(data['latitude'].unique(),1)
    lon = np.around(data['longitude'].unique(),1)
    years = sorted(data['year'].unique())
    s_d = []
    
    # handle duplicate entries (assuming these represent multiple onsets/breakups in a given year)
        # minimum for onset
        # maximum for breakup

    for y in years:
        
        d = data.loc[data['year']==y]
        
        if len(d) == 1:
            
            v = d['value'].item()
            
        elif len(d) > 1:
            
            v = d['value'].max()
        s_d.append(v)
        
    # glrp series
    years = np.asarray(years,
                       dtype=np.int64)
    s_da = xr.DataArray(s_d,
                        coords=[years],
                        dims=['years']).sel(years=slice(1981,2006))
    s_da_mean = s_da.mean(dim='years')
    s_da = s_da - s_da_mean
    
    # era5l series
    e_da_i = era5l_off.sel(lat=lat,
                          lon=lon,
                          method='nearest',
                          tolerance=0.1,
                          drop=True).squeeze()
    nan_check = e_da_i[~np.isnan(e_da_i)]
    
    if len(nan_check) > 0:
        
        e_da_v = e_da_i.values
        e_da_years = e_da_i.time.dt.year.values
        e_da_f = xr.DataArray(e_da_v,
                              coords=[e_da_years],
                              dims=['years'])
        e_da_f = e_da_f.where(e_da_f.years == s_da.years)
        e_da_f_mean = e_da_f.mean(dim='years')
        e_da_f = e_da_f - e_da_f_mean
        
        era5l_off_series.append(e_da_f)
        era5l_off_mean.append(e_da_f_mean.values.item())
        
        glrp_off_series.append(s_da)
        glrp_off_mean.append(s_da_mean.values.item())
        
        lake_locs = lake_locs.append({'siteID':s,'lat':lat.item(),'lon':lon.item()-180}, ignore_index=True)
        
    elif len(nan_check) == 0:
        
        pass
    
    
glrp_dur_series = []
glrp_dur_mean = []
era5l_dur_series = []
era5l_dur_mean = [] 

for s in siteID_dur:
    
    # data per siteID
    data = data_dur.loc[data_dur['siteID']==s]
    lat = np.around(data['latitude'].unique(),1)
    lon = np.around(data['longitude'].unique(),1)
    years = sorted(data['year'].unique())
    s_d = []
    
    # handle duplicate entries (assuming these represent multiple onsets/breakups in a given year)
        # minimum for onset
        # maximum for breakup

    for y in years:
        
        d = data.loc[data['year']==y]
        
        if len(d) == 1:
            
            v = d['value'].item()
            
        elif len(d) > 1:
            
            v = d['value'].mean()
        s_d.append(v)
        
    # glrp series
    years = np.asarray(years,
                       dtype=np.int64)
    s_da = xr.DataArray(s_d,
                        coords=[years],
                        dims=['years']).sel(years=slice(1981,2006))
    s_da_mean = s_da.mean(dim='years')
    s_da = s_da - s_da_mean
    
    # era5l series
    e_da_i = era5l_dur.sel(lat=lat,
                          lon=lon,
                          method='nearest',
                          tolerance=0.1,
                          drop=True).squeeze()
    nan_check = e_da_i[~np.isnan(e_da_i)]
    
    if len(nan_check) > 0:
        
        e_da_v = e_da_i.values
        e_da_years = e_da_i.time.dt.year.values
        e_da_f = xr.DataArray(e_da_v,
                              coords=[e_da_years],
                              dims=['years'])
        e_da_f = e_da_f.where(e_da_f.years == s_da.years)
        e_da_f_mean = e_da_f.mean(dim='years')
        e_da_f = e_da_f - e_da_f_mean
        
        era5l_dur_series.append(e_da_f)
        era5l_dur_mean.append(e_da_f_mean.values.item())
        
        glrp_dur_series.append(s_da)
        glrp_dur_mean.append(s_da_mean.values.item())
        
        lake_locs = lake_locs.append({'siteID':s,'lat':lat.item(),'lon':lon.item()-180}, ignore_index=True)
        
    elif len(nan_check) == 0:
        
        pass    
    
    

dict_era_on = ensembler(era5l_on_series)
dict_era_off = ensembler(era5l_off_series)
dict_era_dur = ensembler(era5l_dur_series)

dict_glrp_on = ensembler(glrp_on_series)
dict_glrp_off = ensembler(glrp_off_series)
dict_glrp_dur = ensembler(glrp_dur_series)

series_on = [dict_era_on,dict_glrp_on]
series_off = [dict_era_off,dict_glrp_off]
series_dur = [dict_era_dur,dict_glrp_dur]

x = 8
y = 8

tser_plotter(series_on,
             series_off,
             series_dur,
             x,
             y,
             xmin,
             xmax,
             ymin,
             ymax,
             labels,
             xticks,
             xtick_labels,
             tick_font,
             title_font,
             axis_font,
             legend_font,
             legend_entrylen,
             legend_entrypad,
             legendcols,
             xlabel_xpos,
             xlabel_ypos,
             xlabel,
             ylabel_xpos,
             ylabel_ypos,
             ylabel,
             colors,
             ub_alpha,
             icevars,
             x0,
             y0,
             xlen,
             ylen)

lake_locs = lake_locs.drop_duplicates()

# =============================================================================
# map of observations
# =============================================================================
    

# plotting data after conversions
proj = ccrs.PlateCarree()

# bbox
x0 = -0.04
y0 = 0.1
xlen = 0.25
ylen = 0.9

# figsize
x=10
y=5

# bounds
new_extent = [-180, 180, 40, 90]

# final data array with all obs-era pairs for significant obs trends
lake_pts = gpd.GeoDataFrame(lake_locs,
                            geometry=gpd.points_from_xy(lake_locs.lon, 
                                                        lake_locs.lat),
                            crs="EPSG:4326")

lake_pts = lake_pts.geometry

map_plotter(x,
            y,
            proj,
            new_extent,
            lake_pts,
            col_glrp,
            x0, 
            y0, 
            xlen, 
            ylen)

output = open('glrp_locs.pkl','wb')
pk.dump(lake_pts,output)
output.close()

