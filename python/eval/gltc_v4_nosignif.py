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




# =============================================================================
# functions
# =============================================================================



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
    
    val = df.loc[(df['lat'] == lat) & (df['lon'] == lon),'arr1'].item()
    
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
    
    frame = {'arr1':[],'arr2':[],'lat':[],'lon':[]}
    df = pd.DataFrame(data=frame)
    vals = arr1.values.flatten()
    data = vals[~np.isnan(vals)]
    for d in data: 
        d_coords = arr1.where(arr1==d,drop=True).squeeze()
        lat = round(d_coords.lat.values.item(),1)
        lon = round(d_coords.lon.values.item(),1)
        e = pixel(arr2,
                  lon,
                  lat,
                  out_arr=False)
        df = df.append({'arr1':d,'arr2':e,'lat':lat,'lon':lon}, ignore_index=True)
        
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
                zorder=4)
    
    # plot mean line
    h = ax.plot(time, 
                ens_roll,
                lw=lw_roll, 
                color=col_fill_a, 
                zorder=3)
            
    return h,ax


def tser_plotter(series_insitu,
                 series_satellite,
                 colors_insitu,
                 colors_satellite,
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
                 ub_alpha,
                 letters):

    f, (ax1,ax2) = plt.subplots(2,1,figsize=(x,y))
    
    time = np.arange(1985,2010)
    
    for s,c in zip(series_insitu,colors_insitu):
        
        h,ax1 = plotter(time,
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
        
    for s,c in zip(series_satellite,colors_satellite):
        
        h,ax2 = plotter(time,
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
    
    count = 0
    for ax in (ax1,ax2):        
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.xaxis.set_ticks(xticks)
        ax.tick_params(labelsize=tick_font,axis="x",direction="in", left="off",labelleft="on")
        ax.tick_params(labelsize=tick_font,axis="y",direction="in")
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.grid(color='0.8', linestyle='dashed', linewidth=0.5)
        ax.xaxis.grid(color='0.8', linestyle='dashed', linewidth=0.5)
        ax.set_axisbelow(True)
        ax.set_title(letters[count],loc='left',fontsize=title_font,fontweight='bold')
        count += 1
    
    
    ax1.xaxis.set_ticklabels([])
    ax2.xaxis.set_ticklabels(xtick_labels)
    
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
    f.text(xlabel_xpos, xlabel_ypos, xlabel, ha='center', fontsize=axis_font)
    f.text(ylabel_xpos, ylabel_ypos, ylabel, va='center', rotation='vertical', fontsize=axis_font)
    
# =============================================================================
#     f.savefig('gltc_tseries.png',bbox_inches='tight',dpi=200)
# =============================================================================
    
    
def map_plotter(proj,
                extent,
                insitu_pts,
                satellite_pts,
                col_insitu,
                col_satellite,
                lab_insitu,
                lab_satellite):
    
    f = plt.figure(figsize=(10,5))
    proj = ccrs.Robinson()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    insitu_pts.plot(ax=ax,
                    markersize=4,
                    color=col_insitu['mean'],
                    zorder=2,
                    transform=ccrs.PlateCarree())
    satellite_pts.plot(ax=ax,
                       markersize=4,
                       color=col_satellite['mean'],
                       zorder=2,
                       transform=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.LAND, 
                   zorder=1, 
                   edgecolor='black')
    
    legend_handles = [Line2D([0], [0],
                             marker='o',
                             color='w',
                             label=lab_insitu,
                             markerfacecolor=col_insitu['mean']),
                      Line2D([0], [0],
                             marker='o',
                             color='w',
                             label=lab_satellite,
                             markerfacecolor=col_satellite['mean'])]
    
    ax.legend(handles=legend_handles,
              frameon=False)

# =============================================================================
#     f.savefig('gltc_locations.png',bbox_inches='tight',dpi=200)
# =============================================================================

def map_plotter_test(proj,
                     extent,
                     glrp_pts,
                     col_glrp):
    
    f = plt.figure(figsize=(10,5))
    proj = ccrs.Robinson()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    glrp_pts.plot(ax=ax,
                  markersize=4,
                  color=col_glrp['mean'],
                  zorder=2,
                  transform=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.LAND, 
                   zorder=1, 
                   edgecolor='black')
    
    legend_handles = [Line2D([0], [0],
                             marker='o',
                             color='w',
                             label='GLRP locations',
                             markerfacecolor=col_glrp['mean'])]
    
    ax.legend(handles=legend_handles,
              frameon=False)

# =============================================================================
#     f.savefig('gltc_locations_og.png',bbox_inches='tight',dpi=200)
# =============================================================================

def c(x):
   col = plt.cm.Greys(x)
   fig, ax = plt.subplots(figsize=(1,1))
   fig.set_facecolor(col)
   ax.axis("off")
   plt.show()



# =============================================================================
# settings
# =============================================================================



title_font = 9
tick_font = 8
axis_font = 9
legend_font = 8

#========== LINE THICKNESS ==========#

# mean line thickness
lw_mean = 1.5
lw_roll = 0.75

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

col_insitu = {}
col_satellite = {}
col_era = {}

col_insitu['mean'] = plt.cm.YlOrBr(0.9)
col_insitu['fill_a'] = plt.cm.YlOrBr(0.7)
col_insitu['fill_b'] = plt.cm.YlOrBr(0.4)
col_satellite['mean'] = plt.cm.Greens(0.9)
col_satellite['fill_a'] = plt.cm.Greens(0.7)
col_satellite['fill_b'] = plt.cm.Greens(0.4)
col_era['mean'] = plt.cm.Greys(0.9)
col_era['fill_a'] = plt.cm.Greys(0.7)
col_era['fill_b'] = plt.cm.Greys(0.4)

colors_insitu = [col_era,
                 col_insitu]

colors_satellite = [col_era,
                    col_satellite]

ub_alpha = 0.5

#========== AXII ==========#

# figsize = (x,y)
x = 8
y = 5

# subplots_adjust
hspace = 0.5
top = 0.9

ymin = -2   # ymin
ymax = 2    # ymax
xmin = 1985 # xmin
xmax = 2010 # xmax

# x ticks/labels 
xticks = np.arange(1985,2015,5)
xtick_labels = [None,1990,None,2000,None,2010]

# x axis label
xlabel = 'Years'
xlabel_xpos = 0.5
xlabel_ypos = 0.05

# y axis label
ylabel = 'Lake temperature anomaly (°C)'
ylabel_xpos = 0.075
ylabel_ypos = 0.535

# xaxis tick label sharing
axis_share = False

#========== LEGEND ==========#

# labels
lab_insitu = 'GLTC insitu'
lab_satellite = 'GLTC satellite'
lab_era = 'ERA5L'

# bbox
x0 = 0.5
y0 = 1
xlen = 0.5
ylen = 0.9

# space between entries
legend_entrypad = 0.5

# length per entry
legend_entrylen = 0.75

# legend colors
legendcols = [col_era['mean'],
              col_insitu['mean'],
              col_satellite['mean']]

handles = [Line2D([0],[0],linestyle='-',lw=2,color=legendcols[0]),\
           Line2D([0],[0],linestyle='-',lw=2,color=legendcols[1]),\
           Line2D([0],[0],linestyle='-',lw=2,color=legendcols[2])]
# labels
lab_insitu = 'GLTC insitu'
lab_satellite = 'GLTC satellite'
lab_era = 'ERA5L'

labels= [lab_era,
         lab_insitu,
         lab_satellite]

letters = ['a','b','c','d']

# =============================================================================
# retrieve data
# =============================================================================



infile1  ="https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/10001/3/6e52deaa45c1695e7742c923ba04d16b".strip() 
infile1  = infile1.replace("https://","http://")
                 
dt1 =pd.read_csv(infile1,
                 skiprows=1,
                 sep=",",
                 names=["recordID",     
                        "variable",     
                        "year",     
                        "siteID",     
                        "value"])

# Coerce the data into the types specified in the metadata  
dt1.recordID=dt1.recordID.astype('category')  
dt1.variable=dt1.variable.astype('category') 
dt1.year=pd.to_numeric(dt1.year,errors='coerce',downcast='integer')  
dt1.siteID=dt1.siteID.astype('category') 
dt1.value=pd.to_numeric(dt1.value,errors='coerce') 
      
print("Here is a description of the data frame dt1 and number of lines\n")
print(dt1.info())
print("--------------------\n\n")                
print("Here is a summary of numerical variables in the data frame dt1\n")
print(dt1.describe())
print("--------------------\n\n")                
                         
print("The analyses below are basic descriptions of the variables. After testing, they should be replaced.\n")                 

print(dt1.recordID.describe())               
print("--------------------\n\n")
                    
print(dt1.variable.describe())               
print("--------------------\n\n")
                    
print(dt1.year.describe())               
print("--------------------\n\n")
                    
print(dt1.siteID.describe())               
print("--------------------\n\n")
                    
print(dt1.value.describe())               
print("--------------------\n\n")
                                 
                
infile2  ="https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/10001/3/6167b9938e8dc99e9ee75251c70776a9".strip() 
infile2  = infile2.replace("https://","http://")
                 
dt2 =pd.read_csv(infile2, 
                 skiprows=1,
                 sep=","  ,
                 quotechar='"' ,
                 names=["siteID",     
                        "Lake_name",     
                        "Other_names",     
                        "lake_or_reservoir",     
                        "location",     
                        "region",     
                        "latitude",     
                        "longitude",     
                        "geospatial_accuracy_km",     
                        "elevation_m",     
                        "mean_depth_m",     
                        "max_depth_m",     
                        "surface_area_km2",     
                        "volume_km3",     
                        "source",     
                        "sampling_depth",     
                        "sampling_time_of_day",     
                        "time_period",     
                        "contributor"], 
                 encoding = "unicode_escape")

# Coerce the data into the types specified in the metadata  
dt2.siteID=dt2.siteID.astype('category')  
dt2.Lake_name=dt2.Lake_name.astype('category')  
dt2.Other_names=dt2.Other_names.astype('category')  
dt2.lake_or_reservoir=dt2.lake_or_reservoir.astype('category')  
dt2.location=dt2.location.astype('category')  
dt2.region=dt2.region.astype('category') 
dt2.latitude=pd.to_numeric(dt2.latitude,errors='coerce') 
dt2.longitude=pd.to_numeric(dt2.longitude,errors='coerce') 
dt2.geospatial_accuracy_km=pd.to_numeric(dt2.geospatial_accuracy_km,errors='coerce') 
dt2.elevation_m=pd.to_numeric(dt2.elevation_m,errors='coerce') 
dt2.mean_depth_m=pd.to_numeric(dt2.mean_depth_m,errors='coerce') 
dt2.max_depth_m=pd.to_numeric(dt2.max_depth_m,errors='coerce') 
dt2.surface_area_km2=pd.to_numeric(dt2.surface_area_km2,errors='coerce') 
dt2.volume_km3=pd.to_numeric(dt2.volume_km3,errors='coerce')  
dt2.source=dt2.source.astype('category')  
dt2.sampling_depth=dt2.sampling_depth.astype('category')  
dt2.sampling_time_of_day=dt2.sampling_time_of_day.astype('category')  
dt2.time_period=dt2.time_period.astype('category')  
dt2.contributor=dt2.contributor.astype('category') 
      
print("Here is a description of the data frame dt2 and number of lines\n")
print(dt2.info())
print("--------------------\n\n")                
print("Here is a summary of numerical variables in the data frame dt2\n")
print(dt2.describe())
print("--------------------\n\n")                
                         
print("The analyses below are basic descriptions of the variables. After testing, they should be replaced.\n")                 

print(dt2.siteID.describe())               
print("--------------------\n\n")
                    
print(dt2.Lake_name.describe())               
print("--------------------\n\n")
                    
print(dt2.Other_names.describe())               
print("--------------------\n\n")
                    
print(dt2.lake_or_reservoir.describe())               
print("--------------------\n\n")
                    
print(dt2.location.describe())               
print("--------------------\n\n")
                    
print(dt2.region.describe())               
print("--------------------\n\n")
                    
print(dt2.latitude.describe())               
print("--------------------\n\n")
                    
print(dt2.longitude.describe())               
print("--------------------\n\n")
                    
print(dt2.geospatial_accuracy_km.describe())               
print("--------------------\n\n")
                    
print(dt2.elevation_m.describe())               
print("--------------------\n\n")
                    
print(dt2.mean_depth_m.describe())               
print("--------------------\n\n")
                    
print(dt2.max_depth_m.describe())               
print("--------------------\n\n")
                    
print(dt2.surface_area_km2.describe())               
print("--------------------\n\n")
                    
print(dt2.volume_km3.describe())               
print("--------------------\n\n")
                    
print(dt2.source.describe())               
print("--------------------\n\n")
                    
print(dt2.sampling_depth.describe())               
print("--------------------\n\n")
                    
print(dt2.sampling_time_of_day.describe())               
print("--------------------\n\n")
                    
print(dt2.time_period.describe())               
print("--------------------\n\n")
                    
print(dt2.contributor.describe())               
print("--------------------\n\n")



# =============================================================================
# dataframes processed by Luke
# =============================================================================
    

                
string_1 = "Lake_Temp_Summer_InSitu"
string_2 = "Lake_Temp_Summer_Satellite"

data_insitu = dt1.loc[dt1['variable'] == string_1]
data_satellite = dt1.loc[dt1['variable'] == string_2]

metadata = dt2[['siteID','latitude','longitude','surface_area_km2']]

data_insitu = pd.merge(data_insitu,metadata,on="siteID").drop(columns=['siteID','recordID','variable'])
data_satellite = pd.merge(data_satellite,metadata,on="siteID").drop(columns=['siteID','recordID','variable'])

dict_insitu = {}
dict_satellite = {}

years = sorted(data_insitu['year'].unique())

for i in years:
    dict_insitu[str(i)] = data_insitu.loc[data_insitu['year'] == i]
    dict_satellite[str(i)] = data_satellite.loc[data_satellite['year'] == i]


# reformulated target grid
res = 0.1
lons = np.arange(-180,180+res,res)
lats = np.arange(-90,90.1+res,res)



# =============================================================================
# data conversion by Inne 
# =============================================================================



# intialize empty numpy array 
values_insitu = np.empty((len(years),len(lats)-1,len(lons)-1))
values_satellite = np.empty((len(years),len(lats)-1,len(lons)-1))

# loop over years
for i,year in enumerate(years): 

    print(year)
    # select dataframe of certain year
    data_year_insitu = dict_insitu[str(year)] 
    data_year_satellite = dict_satellite[str(year)] 

    # turn pandas dataframes in geopandas with lat and lon as geometry
    gdf_insitu = gpd.GeoDataFrame(data_year_insitu, 
                                  geometry=gpd.points_from_xy(data_year_insitu.longitude, 
                                                              data_year_insitu.latitude), 
                                  crs="EPSG:4326")
    gdf_satellite = gpd.GeoDataFrame(data_year_satellite, 
                                     geometry=gpd.points_from_xy(data_year_satellite.longitude, 
                                                                 data_year_satellite.latitude), 
                                     crs="EPSG:4326")

    fn_insitu='insitu_points_'+str(year)
    fn_satellite='satellite_points_'+str(year)

    # save as shapefile to convert to raster
    gdf_insitu.to_file(fn_insitu+'.shp')
    gdf_satellite.to_file(fn_satellite+'.shp')

    # rasterize grid polygon to tiff file and read in as numpy array
    rasterize('value',lons[0],lons[-1],lats[0],lats[-1], res, fn_insitu)
    rasterize('value',lons[0],lons[-1],lats[0],lats[-1], res, fn_satellite)

    # read rasterized values into numpy array 
    year_values_insitu = read_raster(fn_insitu+'.tiff')
    year_values_satellite = read_raster(fn_satellite+'.tiff')

    # clean up
    os.remove(fn_insitu+'.shp')
    os.remove(fn_insitu+'.cpg')
    os.remove(fn_insitu+'.prj')
    os.remove(fn_insitu+'.shx')
    os.remove(fn_insitu+'.dbf')
    os.remove(fn_insitu+'.tiff')
    os.remove(fn_satellite+'.shp')
    os.remove(fn_satellite+'.cpg')
    os.remove(fn_satellite+'.prj')
    os.remove(fn_satellite+'.shx')
    os.remove(fn_satellite+'.dbf')
    os.remove(fn_satellite+'.tiff')

    # save values in  numpy array
    values_insitu[i,:,:] = year_values_insitu
    # save values in  numpy array
    values_satellite[i,:,:] = year_values_satellite

values_insitu[values_insitu == 0] = np.nan # there are no 0 temperatures in the dataset
values_satellite[values_satellite == 0] = np.nan # there are no 0 temperatures in the dataset


longitudes = np.arange(0,360,res)
latitudes = np.arange(-90,90+res,res)

time = pd.date_range(start='1985-01-01',end='2009-01-01',freq='YS')

# data arrays of insitu and satellite obs
da_insitu = xr.DataArray(values_insitu, coords=[time,latitudes,longitudes], dims=["time", "lat", "lon"])
da_satellite = xr.DataArray(values_satellite, coords=[time,latitudes,longitudes], dims=["time", "lat", "lon"])



# =============================================================================
# slope
# =============================================================================



# slope calculations; mask for significant obs trends only
slope_insitu,pval_insitu = slope_field(da_insitu)
slope_satellite,pval_satellite = slope_field(da_satellite)

# =============================================================================
# slope_insitu_signif = slope_insitu.where(pval_insitu<0.05)
# slope_satellite_signif = slope_satellite.where(pval_satellite<0.05)
# =============================================================================


# insitu obs slope
slope_insitu_jas = slope_insitu.where(slope_insitu.lat >= 23.5,drop=True).squeeze()
slope_insitu_jas_sh = slope_insitu.sel(lat=slice(-23.45,-0.05))

slope_insitu_jfm = slope_insitu.where(slope_insitu.lat <= -23.5,drop=True).squeeze()
slope_insitu_jfm_nh = slope_insitu.sel(lat=slice(0.05,23.45))

# satellite obs slope
slope_satellite_jas = slope_satellite.where(slope_satellite.lat >= 23.5,drop=True).squeeze()
slope_satellite_jas_sh = slope_satellite.sel(lat=slice(-23.45,-0.05))

slope_satellite_jfm = slope_satellite.where(slope_satellite.lat <= -23.5,drop=True).squeeze()
slope_satellite_jfm_nh = slope_satellite.sel(lat=slice(0.05,23.45))


# reanalysis read in
os.chdir("/home/luke/documents/data/gltc/knb-lter-ntl.10001.3/final/")

era5l_jas_file = "era5-land_lakes_lmlt_JAS_1985_2009.nc"
era5l_jas_sh_file = "era5-land_lakes_lmlt_JAS_sh_1985_2009.nc"

era5l_jfm_file = "era5-land_lakes_lmlt_JFM_1985_2009.nc"
era5l_jfm_nh_file = "era5-land_lakes_lmlt_JFM_nh_1985_2009.nc"

da_era5l_jas = xr.open_dataset(era5l_jas_file,decode_times=False).lmlt
da_era5l_jas['time'] = time
da_era5l_jas = da_era5l_jas.rename({'longitude':'lon',
                                    'latitude':'lat'})

da_era5l_jas_sh = xr.open_dataset(era5l_jas_sh_file,decode_times=False).lmlt
da_era5l_jas_sh['time'] = time
da_era5l_jas_sh = da_era5l_jas_sh.rename({'longitude':'lon',
                                    'latitude':'lat'})

da_era5l_jfm = xr.open_dataset(era5l_jfm_file,decode_times=False).lmlt
da_era5l_jfm['time'] = time
da_era5l_jfm = da_era5l_jfm.rename({'longitude':'lon',
                                    'latitude':'lat'})

da_era5l_jfm_nh = xr.open_dataset(era5l_jfm_nh_file,decode_times=False).lmlt
da_era5l_jfm_nh['time'] = time
da_era5l_jfm_nh = da_era5l_jfm_nh.rename({'longitude':'lon',
                                    'latitude':'lat'})


# slope calculations
slope_era5l_jas, _ = slope_field(da_era5l_jas)
slope_era5l_jas_sh, _ = slope_field(da_era5l_jas_sh)

slope_era5l_jfm, _ = slope_field(da_era5l_jfm)
slope_era5l_jfm_nh, _ = slope_field(da_era5l_jfm_nh)

# dataframes for insitu
insitu_jas = arr_to_df(slope_insitu_jas,
                       slope_era5l_jas)
insitu_jas_sh = arr_to_df(slope_insitu_jas_sh,
                          slope_era5l_jas_sh)
insitu_jfm = arr_to_df(slope_insitu_jfm,
                         slope_era5l_jfm)
insitu_jfm_nh = arr_to_df(slope_insitu_jfm_nh,
                         slope_era5l_jfm_nh)

# dataframes for satellite
satellite_jas = arr_to_df(slope_satellite_jas,
                       slope_era5l_jas)
satellite_jas_sh = arr_to_df(slope_satellite_jas_sh,
                          slope_era5l_jas_sh)
satellite_jfm = arr_to_df(slope_satellite_jfm,
                         slope_era5l_jfm)
satellite_jfm_nh = arr_to_df(slope_satellite_jfm_nh,
                         slope_era5l_jfm_nh)

frames_insitu = [insitu_jas,
                 insitu_jas_sh,
                 insitu_jfm,
                 insitu_jfm_nh]

frames_satellite = [satellite_jas,
                    satellite_jas_sh,
                    satellite_jfm,
                    satellite_jfm_nh]

# final data array with all obs-era pairs for significant obs trends
gltc_insitu = pd.concat(frames_insitu)
gltc_satellite = pd.concat(frames_satellite)



# =============================================================================
# slope scatterplot
# =============================================================================



# adding columns for source and recombining for hue
insitu_cp = cp.deepcopy(insitu_jas)
satellite_cp = cp.deepcopy(satellite_jas)
insitu_cp['source'] = 'insitu'
satellite_cp['source'] = 'satellite'
dltc_cp = pd.concat([insitu_cp,satellite_cp])

slope, intercept, r_value, p_value, std_err = sts.linregress(gltc_df['arr2'],
                                                             gltc_df['arr1'])
r_sq = r_value**2
label = r'$R^2:{0:.2f}$'.format(r_sq)

ax = sb.lmplot(data=gltc_df,
               x='arr2',
               y='arr1',
               line_kws={'label':"r'$R^2:{0:.2f}$'".format(r_sq)},
               legend=True)
ax.set(xlabel = "ERA5L slope",
       ylabel = "GLTC slope")



# =============================================================================
# time series plots
# =============================================================================



# only take series from era5_jas (no other latitude ranges contribute to significant slopes in obs)
era_insitu_series_list = []
era_satellite_series_list= []
insitu_series_list = []
satellite_series_list = []

# get series' for pixels in insitu data
for lon,lat in zip(gltc_insitu['lon'].values,gltc_insitu['lat'].values):
    
    era_series = pixel(da_era5l_jas,
                       lon,
                       lat,
                       out_arr=True)
    
    era_series = era_series - 273.15
    era_series = era_series - era_series.mean(dim='time')
    
    insitu_series =  df_indexer(slope_insitu_jas,
                                da_insitu,
                                insitu_jas,
                                lon,
                                lat)
    insitu_series = insitu_series - insitu_series.mean(dim='time')
    
    era_series = era_series.where(era_series.time == insitu_series.time)
    
    era_insitu_series_list.append(era_series)
    insitu_series_list.append(insitu_series)
 
# series' for pixels in satellite jas data
for lon,lat in zip(satellite_jas['lon'].values,satellite_jas['lat'].values):
    
    era_series = pixel(da_era5l_jas,
                       lon,
                       lat,
                       out_arr=True)
    
    era_series = era_series - 273.15
    era_series = era_series - era_series.mean(dim='time')
    
    satellite_series =  df_indexer(slope_satellite_jas,
                                   da_satellite,
                                   satellite_jas,
                                   lon,
                                   lat)
    satellite_series = satellite_series - satellite_series.mean(dim='time')
    
    era_series = era_series.where(era_series.time == satellite_series.time)
    
    era_satellite_series_list.append(era_series)
    satellite_series_list.append(satellite_series)
    
    
# dataframes for satellite
satellite_jas = arr_to_df(slope_satellite_jas,
                       slope_era5l_jas)
satellite_jas_sh = arr_to_df(slope_satellite_jas_sh,
                          slope_era5l_jas_sh)
satellite_jfm = arr_to_df(slope_satellite_jfm,
                         slope_era5l_jfm)
satellite_jfm_nh = arr_to_df(slope_satellite_jfm_nh,
                         slope_era5l_jfm_nh)

# series' for pixels in satellite jas data
for lon,lat in zip(satellite_jas['lon'].values,satellite_jas['lat'].values):
    
    era_series = pixel(da_era5l_jas,
                       lon,
                       lat,
                       out_arr=True)
    
    era_series = era_series - 273.15
    era_series = era_series - era_series.mean(dim='time')
    
    satellite_series =  df_indexer(slope_satellite_jas,
                                   da_satellite,
                                   satellite_jas,
                                   lon,
                                   lat)
    satellite_series = satellite_series - satellite_series.mean(dim='time')
    
    era_series = era_series.where(era_series.time == satellite_series.time)
    
    era_satellite_series_list.append(era_series)
    satellite_series_list.append(satellite_series)
    

dict_era_insitu = ensembler(era_insitu_series_list)
dict_era_satellite = ensembler(era_satellite_series_list)
dict_insitu = ensembler(insitu_series_list)
dict_satellite = ensembler(satellite_series_list)


series_insitu = [dict_era_insitu,
                 dict_insitu]

series_satellite = [dict_era_satellite,
                    dict_satellite]

tser_plotter(series_insitu,
             series_satellite,
             colors_insitu,
             colors_satellite,
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
             ub_alpha,
             letters)



# =============================================================================
# map of observations
# =============================================================================


os.chdir('/home/luke/documents/data/gltc/knb-lter-ntl.10001.3/final/')
proj = ccrs.Robinson()
extent = [-180, 180, -65, 90]

# dataframes for insitu
insitu_jas = arr_to_df(slope_insitu_jas,
                       slope_era5l_jas)
insitu_jas_sh = arr_to_df(slope_insitu_jas_sh,
                          slope_era5l_jas_sh)
insitu_jfm = arr_to_df(slope_insitu_jfm,
                         slope_era5l_jfm)
insitu_jfm_nh = arr_to_df(slope_insitu_jfm_nh,
                         slope_era5l_jfm_nh)

# dataframes for satellite
satellite_jas = arr_to_df(slope_satellite_jas,
                       slope_era5l_jas)
satellite_jas_sh = arr_to_df(slope_satellite_jas_sh,
                          slope_era5l_jas_sh)
satellite_jfm = arr_to_df(slope_satellite_jfm,
                         slope_era5l_jfm)
satellite_jfm_nh = arr_to_df(slope_satellite_jfm_nh,
                         slope_era5l_jfm_nh)

insitu_jas['lon'] = insitu_jas.apply(lambda row: row.lon - 180,axis=1)
insitu_jas_sh['lon'] = insitu_jas_sh.apply(lambda row: row.lon - 180,axis=1)
insitu_jfm['lon'] = insitu_jfm.apply(lambda row: row.lon - 180,axis=1)
insitu_jfm_nh['lon'] = insitu_jfm_nh.apply(lambda row: row.lon - 180,axis=1)

insitu = pd.concat([insitu_jas,
                    insitu_jas_sh,
                    insitu_jfm,
                    insitu_jfm_nh])

insitu_pts = gpd.GeoDataFrame(insitu,
                              geometry=gpd.points_from_xy(insitu.lon, 
                                                          insitu.lat),
                              crs="EPSG:4326")

insitu_pts = insitu_pts.geometry

satellite_jas['lon'] = satellite_jas.apply(lambda row: row.lon - 180,axis=1)
satellite_jas_sh['lon'] = satellite_jas_sh.apply(lambda row: row.lon - 180,axis=1)
satellite_jfm['lon'] = satellite_jfm.apply(lambda row: row.lon - 180,axis=1)
satellite_jfm_nh['lon'] = satellite_jfm_nh.apply(lambda row: row.lon - 180,axis=1)

satellite = pd.concat([satellite_jas,
                       satellite_jas_sh,
                       satellite_jfm,
                       satellite_jfm_nh])

satellite_pts = gpd.GeoDataFrame(satellite,
                                 geometry=gpd.points_from_xy(satellite.lon, 
                                                          satellite.lat),
                                 crs="EPSG:4326")

satellite_pts = satellite_pts.geometry

map_plotter(proj,
            extent,
            insitu_pts,
            satellite_pts,
            col_insitu,
            col_satellite,
            lab_insitu,
            lab_satellite)

# testing with data at original coordinates (after some subsampling)
# =============================================================================
# proj = ccrs.Robinson()
# extent = [-180, 180, 20, 90]
# data = pd.concat([data_insitu,data_satellite])
# og_pts = gpd.GeoDataFrame(data,
#                           geometry=gpd.points_from_xy(data.longitude, 
#                                                         data.latitude),
#                           crs="EPSG:4326")
# 
# og_pts = og_pts.geometry
# 
# og_pts = gdf_insitu.geometry
# map_plotter_test(proj,
#                  extent,
#                  og_pts,
#                  col_insitu)
# 
# # test locations of data once in dataframe of slopes
# insitu_jas = arr_to_df(slope_insitu_jas,
#                        slope_era5l_jas)
# satellite_jas = arr_to_df(slope_satellite_jas,
#                        slope_era5l_jas)
# shapefile_data = pd.concat([insitu_jas,satellite_jas])
# 
# shapefile_data['lon'] = shapefile_data.apply(lambda row: row.lon - 180,axis=1)
# 
# 
# shapefile_gdf = gpd.GeoDataFrame(shapefile_data,
#                                  geometry=gpd.points_from_xy(shapefile_data.lon, 
#                                                              shapefile_data.lat),
#                                  crs="EPSG:4326")
# sfile='shapefile_test_3'
# shapefile_gdf.to_file(sfile+'.shp')
# =============================================================================
