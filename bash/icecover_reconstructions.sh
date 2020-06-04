#!/bin/bash -l

# =======================================================================
# SUMMARY
# =======================================================================


# Processing for reconstructed ice cover


# =======================================================================
# INITIALIZATION
# =======================================================================


# set output directory
outDIR=


# user scratch directory
workDIR=

# isimip grid
gridfile=


# set starting directory
inDIR=


# years
years=('1981_1982' '1983_1984' '1985_1986' '1987_1988' '1989_1990' '1991_1992' '1993_1994' '1995_1996' '1997_1998' '1999_2000' '2001_2005' '2006_2007' '2008_2009' '2010_2011' '2012_2013' '2014_2015' '2016_2017' '2018_2019')


# ==============================================================================
# PROCESSING
# ==============================================================================


cd $inDIR

for year in "${years[@]}"; do

    cdo -b F64 -O -L setreftime,1981-01-01,00:00:00,1days -settaxis,$((${year:0:4}+0))-01-01,00:00:00,1days -setctomiss,0 -gtc,0.0001 -daymin era5-land_lakes_icedepth_6hourly_${year}.nc $workDIR/step1_${year}.nc

done

cdo -b F64 mergetime $workDIR/step1_*.nc $workDIR/step2.nc

rm $workDIR/step1_*.nc

# remap starting file
cdo remapnn,$gridfile $workDIR/step2.nc $workDIR/step3.nc

rm $workDIR/step2.nc


# ==============================================================================
# ICE START
# ==============================================================================


#marker for stderr ice start
echo ' '
echo 'ICE START CALC'
echo ' '


# select October to December for 1st to 2nd-last year
cdo -b F64 -O -L setctomiss,0 -muldoy -selmon,10/12 -seldate,1981-01-01T00:00:00,2019-01-01T00:00:00 $workDIR/step3.nc $workDIR/step3_a.nc
# select January to September for 2nd to last year (lag by 365)
cdo -b F64 -O -L addc,365 -setctomiss,0 -muldoy -selmon,1/9 -seldate,1982-01-01T00:00:00,2019-12-31T00:00:00 $workDIR/step3.nc $workDIR/step3_b.nc
# merge selections
cdo -b F64 mergetime $workDIR/step3_a.nc $workDIR/step3_b.nc $workDIR/step4.nc

rm $workDIR/step3_a.nc
rm $workDIR/step3_b.nc

for i in $(seq 1981 2018); do

    cdo -b F64 -O -L timmin -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/step4.nc $workDIR/step4_$i.nc

done

rm $workDIR/step4.nc

cdo -b F64 -O mergetime $workDIR/step4_*.nc $workDIR/step5.nc

cdo ifthen $gridfile $workDIR/step5.nc $workDIR/step6.nc

rm $workDIR/step4_*.nc

# set attributes
cdo -b F64 -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years -setattribute,icestart@long_name='First day of lake ice cover' -setname,'icestart' -setunit,'day of hydrological year' $workDIR/step6.nc $outDIR/era5-land_icestart_1981_2018.nc
# take fldmean
cdo -b F64 -O fldmean $outDIR/era5-land_icestart_1981_2019.nc $outDIR/era5-land_icestart_fldmean_1981_2019.nc

rm $workDIR/step5.nc
rm $workDIR/step6.nc


# ==============================================================================
# ICE END
# ==============================================================================


echo ' '
echo 'ICE END CALC'
echo ' '


# select September to December for 1st to 2nd-last year
cdo -b F64 -O -L setctomiss,0 -muldoy -selmon,9/12 -seldate,1981-01-01T00:00:00,2019-01-01T00:00:00 $workDIR/step3.nc $workDIR/step3_a.nc
# select January to August for 2nd to last year (lag by 365)
cdo -b F64 -O -L addc,365 -setctomiss,0 -muldoy -selmon,1/8 -seldate,1982-01-01T00:00:00,2019-12-31T00:00:00 $workDIR/step3.nc $workDIR/step3_b.nc
# merge selections
cdo -b F64 mergetime $workDIR/step3_a.nc $workDIR/step3_b.nc $workDIR/step4.nc

rm $workDIR/step3_a.nc
rm $workDIR/step3_b.nc

for i in $(seq 1981 2018); do

    cdo -b F64 -O -L timmax -seldate,$i-09-01T00:00:00,$(($i+1))-08-31T00:00:00 $workDIR/step4.nc $workDIR/step4_$i.nc

done

rm $workDIR/step4.nc

cdo -b F64 -O mergetime $workDIR/step4_*.nc $workDIR/step5.nc

cdo ifthen $gridfile $workDIR/step5.nc $workDIR/step6.nc

rm $workDIR/step4_*.nc

# set attributes
cdo -b F64 -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years -setattribute,iceend@long_name='Last day of lake ice cover' -setname,'iceend' -setunit,'day of hydrological year' $workDIR/step6.nc $outDIR/era5-land_iceend_1981_2019.nc
# take fldmean
cdo -b F64 -O fldmean $outDIR/era5-land_iceend_1981_2019.nc $outDIR/era5-land_iceend_fldmean_1981_2019.nc

rm $workDIR/step5.nc
rm $workDIR/step6.nc


# ==============================================================================
# DURATION
# ==============================================================================


echo ' '
echo 'ICE DURATION CALC'
echo ' '

for i in $(seq 1981 2018); do

    cdo -b F64 -O -L timsum -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/step3.nc $workDIR/step4_$i.nc

done

cdo -b F64 -O mergetime $workDIR/step4_*.nc $workDIR/step5.nc

cdo ifthen $gridfile $workDIR/step5.nc $workDIR/step6.nc

# set attributes
cdo -b F64 -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years -setattribute,icedur@long_name='Duration of lake ice cover (timsum)' -setname,'icedur' -setunit,'days' $workDIR/step6.nc $outDIR/era5-land_icedur_1981_2019.nc
# take fldmean
cdo -b F64 fldmean $outDIR/era5-land_icedur_1981_2019.nc $outDIR/era5-land_icedur_fldmean_1981_2019.nc



# ==============================================================================
# TEMPORAL MEANS
# ==============================================================================


cdo -b F64 timmean $outDIR/era5-land_icestart_1981_2019.nc $outDIR/era5-land_icestart_timmean_1981_2019.nc
cdo -b F64 timmean $outDIR/era5-land_iceend_1981_2019.nc $outDIR/era5-land_iceend_timmean_1981_2019.nc
cdo -b F64 timmean $outDIR/era5-land_icedur_1981_2019.nc $outDIR/era5-land_icedur_timmean_1981_2019.nc


# ==============================================================================
# CLEANUP
# ==============================================================================

rm $workDIR/step*.nc
