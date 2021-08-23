#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:03:59 2021

@author: luke
"""

# for ice and watertemp esa_cci netcdfs:
    # read in mask of ESA CCI for valid sample locations
    # extract series from ESA CCI and ERA5L with locations
    # use comparison processing similar to glrp_v3.py 
    
    
# saturday; add shit to interpolate nas and limit lakes to those with greater than 15 years?
# possible to do the time limit on original das before tyhe means? same with interp?

# have done timstep selection of sim via obs w/ where + interpolate
# can also do nan timsum on obs subsets to find which pixels have greater than 10-15 time steps (try on sunday)

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



def reader(file,
                 var):
    
    da = xr.open_dataset(file,decode_times=False)
    da = da[var]
    time = pd.date_range(start='1993-01-01',
                         end='2019-01-01',
                         freq='YS')
    da['time'] = time
    da = da.rename({'longitude':'lon',
                    'latitude':'lat'})
    
    return da


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
    ens_roll = ens_mean.rolling(time=5, center=True).mean()
    dict_ens = {}
    dict_ens['mean'] = ens_mean
    dict_ens['std'] = ens_std
    dict_ens['max'] = ens_max
    dict_ens['min'] = ens_min
    dict_ens['roll'] = ens_roll
    return dict_ens



def plotter(time,
            da,
            ax,
            lw_mean,
            lw_roll,
            col,
            ub_alpha):
            
        
    roll = da.rolling(time=5, 
                      center=True).mean()
    
    # plot mean line
    h = ax.plot(time, 
                da,
                lw=lw_mean, 
                color=col,
                zorder=3)
    
    # plot rolling mean line
    h = ax.plot(time, 
                roll,
                lw=lw_roll, 
                color=col, 
                zorder=4)
            
    return h,ax


def tser_plotter(subsets,
                 sim_series,
                 obs_series,
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
                 x0,
                 y0,
                 xlen,
                 ylen):

    f, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,figsize=(x,y))
    for s,ax in zip(subsets,(ax1,ax2,ax3,ax4)):
        
        # obs
        da = obs_series[s]
        time = da.time.dt.year.values
        h,ax = plotter(time,
                       da,
                       ax,
                       lw_mean,
                       lw_roll,
                       colors['obs'],
                       ub_alpha)
        
        # sim
        da = sim_series[s]
        time = da.time.dt.year.values
        h,ax = plotter(time,
                       da,
                       ax,
                       lw_mean,
                       lw_roll,
                       colors['sim'],
                       ub_alpha)
        
    i = 0    
    for ax in (ax1,ax2,ax3,ax4):        
        ax.set_xlim(xmin,xmax)
        ax.autoscale(axis='y')
        ax.xaxis.set_ticks(xticks)
        ax.tick_params(labelsize=tick_font,axis="x",direction="in", left="off",labelleft="on")
        ax.tick_params(labelsize=tick_font,axis="y",direction="in")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_axisbelow(True)
        ax.set_title(letters[i],loc='left',fontsize=title_font,fontweight='bold')
        i += 1
        
    ax1.xaxis.set_ticklabels([]) 
    ax2.xaxis.set_ticklabels([]) 
    ax3.xaxis.set_ticklabels([]) 
    ax4.xaxis.set_ticklabels(xtick_labels)
    
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
    
    f.text(ylabel_xpos, 
           ylabel_ypos, 
           ylabel, 
           va='center', 
           rotation='vertical', 
           fontsize=axis_font)
    
    f.savefig('si_f34.png',bbox_inches='tight',dpi=200)
    
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
                  marker='o',
                  markersize=1,
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
                             label='ESA CCI lakes',
                             markerfacecolor=col_glrp['mean'])]
    
    ax.legend(handles=legend_handles,
              bbox_to_anchor=(x0, y0, xlen, ylen),
              frameon=False)
    
    f.savefig('esa_cci_lswt_locations_final.png',bbox_inches='tight',dpi=200)


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
legendcols = [plt.cm.Greys(0.9),
              plt.cm.Blues(0.9)]

colors = {}
colors['obs'] = plt.cm.Blues(0.9)
colors['sim'] = plt.cm.Greys(0.9)

ub_alpha = 0.5

#========== AXII ==========#

# figsize = (x,y)
x = 8
y = 8

# subplots_adjust
hspace = 0.5
top = 0.9

ymin = -2   # ymin
ymax = 2    # ymax
xmin = 1995 # xmin
xmax = 2019 # xmax

# x ticks/labels 
xticks = np.arange(1995,2025,5)
xtick_labels = [1995,None,2005,None,2015,None]

# x axis label
xlabel = 'Years'
xlabel_xpos = 0.5
xlabel_ypos = 0.075

# y axis label
ylabel = 'Lake temperature anomaly (°C)'
ylabel_xpos = 0.05
ylabel_ypos = 0.535


# xaxis tick label sharing
axis_share = False

#========== LEGEND ==========#

# labels
lab_sim = 'ERA5L'
lab_obs = 'ESA CCI'
labels = [lab_sim,
          lab_obs]

# bbox
x0 = 0.75
y0 = 1.0
xlen = 0.25
ylen = 0.9

# space between entries
legend_entrypad = 0.5

# length per entry
legend_entrylen = 0.75

handles = [Line2D([0],[0],linestyle='-',lw=2,color=legendcols[0]),\
           Line2D([0],[0],linestyle='-',lw=2,color=legendcols[1])]
    
letters = ['a) 90°N - 23.5°N',
           'b) 23.5°N - 0°',
           'c) 0° - 23.5°S',
           'd) 23.5°S - 90°S']

    
# =============================================================================
# retrieve data
# =============================================================================


y1 = 1993
y2 = 2019

subsets = ['JAS_nh','JFM_eq','JAS_eq','JFM_sh']

os.chdir('/home/luke/documents/data/esa_cci/')

obs_jas_nh_file = "lswt_"+str(y1)+"_"+str(y2)+"_JAS_nh.nc"
obs_jas_eq_file = "lswt_"+str(y1)+"_"+str(y2)+"_JAS_eq.nc"
obs_jfm_nh_file = "lswt_"+str(y1)+"_"+str(y2)+"_JFM_sh.nc"
obs_jfm_eq_file = "lswt_"+str(y1)+"_"+str(y2)+"_JFM_eq.nc"

sim_jas_nh_file = "era5l_lmlt_"+str(y1)+"_"+str(y2)+"_JAS_nh.nc"
sim_jas_eq_file = "era5l_lmlt_"+str(y1)+"_"+str(y2)+"_JAS_eq.nc"
sim_jfm_nh_file = "era5l_lmlt_"+str(y1)+"_"+str(y2)+"_JFM_sh.nc"
sim_jfm_eq_file = "era5l_lmlt_"+str(y1)+"_"+str(y2)+"_JFM_eq.nc"

files = [obs_jas_nh_file,
         obs_jas_eq_file,
         obs_jfm_nh_file,
         obs_jfm_eq_file,
         sim_jas_nh_file,
         sim_jas_eq_file,
         sim_jfm_nh_file,
         sim_jfm_eq_file]

obs_das = {}
sim_das = {}

masks = {}
lake_locs = {}

obs_series = {}
obs_means = {}
sim_series = {}
sim_means = {}

for file in files:
    if 'lswt' in file:
        for subset in subsets:
            if subset in file:
                var='lake_surface_water_temperature'
                obs_das[subset] = reader(file,
                                         var)
    elif 'lmlt' in file:
        for subset in subsets:
            if subset in file:
                var='lmlt'
                sim_das[subset] = reader(file,
                                         var)

frame = {'lat':[],'lon':[]}
lake_locs = pd.DataFrame(data=frame)                
for subset in subsets:
    
    # generate mask + lake locations
    sim = sim_das[subset].mean(dim='time')
    sim = sim.where(sim.isnull(),1)
    obs = obs_das[subset]
    obs = obs.where(obs.isnull(),1)
    obs = obs.sum(dim='time').squeeze()
    obs = obs.where(obs >= 15)
    obs = obs.where(obs.isnull(),1)
    da = obs.where(sim == 1)
    masks[subset] = da
    msk = masks[subset]
    locs = msk.where(msk == 1,
                     drop=True).squeeze()
    locs = locs.to_dataframe().reset_index()
    locs = locs.dropna()
    locs = locs.drop(columns='lake_surface_water_temperature')
    locs = locs.drop_duplicates()
    lake_locs = lake_locs.append(locs)
    
    # observed data - interpolate missing steps
    series = obs_das[subset].where(msk == 1).interpolate_na(dim='time')
    obs_mean = series.mean(dim='time')
    tm = series - obs_mean
    obs_series[subset] = tm.mean(dim=['lat','lon'])
    
    # sim data - interpolate missing steps?
    # by masking by obs_das[subset], i think i am masking missing time steps for pixels in obs
    series = sim_das[subset].where(obs_das[subset]>0).where(msk == 1).interpolate_na(dim='time')
# =============================================================================
#     series = sim_das[subset].where(msk == 1).interpolate_na(dim='time')
# =============================================================================
    sim_mean = series.mean(dim='time')
    tm = series - sim_mean
    sim_series[subset] = tm.mean(dim=['lat','lon'])


x = 8
y = 10

tser_plotter(subsets,
             sim_series,
             obs_series,
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
             x0,
             y0,
             xlen,
             ylen)


# =============================================================================
# map of observations
# =============================================================================
    

# plotting data after conversions
proj = ccrs.PlateCarree()

# bbox
x0 = -0.04
y0 = 0.125
xlen = 0.25
ylen = 0.9

# figsize
x=10
y=5

# bounds
new_extent = [-180, 180, -60, 90]

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

# =============================================================================
# output = open('esa_cci_lswt_locs.pkl','wb')
# pk.dump(lake_pts,output)
# output.close()
# =============================================================================

