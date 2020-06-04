#!/bin/bash -l

# =======================================================================
# SUMMARY
# =======================================================================


# Processing for reconstructed mixlayertemps


# =======================================================================
# INITIALIZATION
# =======================================================================


# set output directory
outDIR=


# set starting directory
inDIR=


# set processing directory
workDIR=


# lake mask
maskfile=


# years
years=('1981_1985' '1986_1990' '1991_1995' '1996_2000' '2001_2005' '2006_2010' '2011_2015' '2016_2019')


# ==============================================================================
# PROCESSING
# ==============================================================================


# change to input directory
cd $inDIR


for year in "${years[@]}"; do

    cdo monmean era5-land_lakes_mixlayertemp_6hourly_${year}.nc $workDIR/step1_monthly_${year}.nc
    cdo yearmean era5-land_lakes_mixlayertemp_6hourly_${year}.nc $workDIR/step1_${year}.nc
    cdo ifthen $maskfile $workDIR/step1_monthly_${year}.nc $workDIR/step2_monthly_${year}.nc
    cdo ifthen $maskfile $workDIR/step1_${year}.nc $workDIR/step2_${year}.nc
    
done

cdo mergetime $workDIR/step2_monthly_*.nc $workDIR/step3_monthly.nc
cdo mergetime $workDIR/step2*.nc $workDIR/step3.nc
cdo -L setcalendar,proleptic_gregorian -setreftime,1661-01-01,00:00:00,1month -settaxis,1981-01-01,00:00:00,1month -selyear,1981/2019 $workDIR/step3_monthly.nc $outDIR/era5-land_lakes_mixlayertemp_monthly_1981_2019.nc
cdo -L setcalendar,proleptic_gregorian -setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years -selyear,1981/2019 $workDIR/step3.nc $outDIR/era5-land_lakes_mixlayertemp_1981_2019.nc
cdo timmean $outDIR/era5-land_lakes_mixlayertemp_1981_2019.nc $outDIR/era5-land_lakes_mixlayertemp_timmean_1981_2019.nc
cdo fldmean $outDIR/era5-land_lakes_mixlayertemp_1981_2019.nc $outDIR/era5-land_lakes_mixlayertemp_fldmean_1981_2019.nc

rm $workDIR/*.nc

cdo -L timmean -selyear,1981/1999 $outDIR/era5-land_lakes_mixlayertemp_fldmean_1981_2019.nc $workDIR/base.nc
cdo sub $outDIR/era5-land_lakes_mixlayertemp_fldmean_1981_2019.nc $workDIR/base.nc $outDIR/era5-land_lakes_mixlayertemp_fldmean_anomaly_1981_2019.nc


rm $workDIR/*.nc


        




