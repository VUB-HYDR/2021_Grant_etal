#!/bin/bash -l

# =======================================================================
# SUMMARY
# =======================================================================


# Daily to annual field means processing for ice thickness to ice duration


# =======================================================================
# INITIALIZATION
# =======================================================================


# set output directory
outDIR=


# set starting directory
inDIR=


# set starting directory
workDIR=


# Define settings
flag_sect=1;   	# 0: all sectors
				# 1: lakes_global


flag_models=3;  # 0: CLM45
				# 1: ALBM
				# 2: LAKE
				# 3: SIMSTRAT-UoG
				# 4: VIC-LAKE


flag_scen=2;    # 0: pre-industrial
                # 1: historical
                # 2: future


flag_tstep=0;	# 0: daily
				# 1: monthly
				# 2: annual


flag_var=2;     # 0: watertemp
                # 1: lakeicefrac
                # 2: icethick


flag_endvar=2;  # 0: icestart
                # 1: iceend
                # 2: icedur


flag_opt=0;     # 0: 2005soc_co2
                # 1: 1860soc_co2
                # 2: nosoc_co2

flag_prod=0;    # 0: fldmean
                # 1: sig
                # 2: eval


# define all possible sectors
sectors=("all" "lakes_global")


# define all possible models (top list; folder style, bottom list; file name style)
models=("CLM45" "ALBM" "LAKE" "SIMSTRAT-UoG" "VIC-LAKE")


model_fnames=("clm45" "albm" "lake" "simstrat-uog" "vic-lake")


# scenario folders
scenario_folders=("pre-industrial" "historical" "future")


# scenario filenames
scenarios=("picontrol" "rcp26" "rcp60" "rcp85")


#define all RCPs for future scenario
rcps=("rcp26" "rcp60" "rcp85")


# define lake variables
variables=("watertemp" "lakeicefrac" "icethick")


#define end product lake variables
endvariables=("icestart" "iceend" "icedur")


# define forcing
forcing=("gfdl-esm2m" "hadgem2-es" "ipsl-cm5a-lr" "miroc5")


# define timestep
timesteps=("daily" "monthly")


# define all options
options=("2005soc_co2" "1860soc_co2" "nosoc_co2")


# define end products
products=("fldmean" "sig" "eval")


# set run settings based on flags
tstep=${timesteps[$flag_tstep]}
model=${models[$flag_models]}
model_fname=${model_fnames[$flag_models]}
var=${variables[$flag_var]}
prod=${products[$flag_prod]}
scen_folder=${scenario_folders[$flag_scen]}


# ==============================================================================
# FUNCTIONS
# ==============================================================================


# first turn off history substitution (free up use of "!")
set +H

# function for setting day 1 = 1 of Cctober hydro-year
isleap_oct(){

    year=$1
            (( !(year % 4) && ( year % 100 || !(year % 400) ) )) &&
            day1=275 || day1=274
}

# function for setting day 1 = 1 of September hydro-year
isleap_sep(){

    year=$1
            (( !(year % 4) && ( year % 100 || !(year % 400) ) )) &&
            day1=245 || day1=244
}


# function for calculating ice indexes
index_calculator(){

    if [ "$endvar" == "icedur" ]
    then

            if [ "$var" == "lakeicefrac" ]
            then
                    # $1: $workDIR/step1.nc from initial file merge
                    cdo -O -L setctomiss,0 -gtc,0 -sellevidx,1 $1 $workDIR/dummy_${force}_step2.nc

            elif [ "$var" == "icethick" ]
            then

                    cdo -O -L setctomiss,0 -gtc,0 $1 $workDIR/dummy_${force}_step2.nc

            fi

            # ice dur loop (1st two elements of timeranges list)
            for i in $(seq $((${endper:0:4}+0)) $((${endper:5:4}-1))); do

                cdo -O -L fldmean -timsum -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/dummy_${force}_step2.nc $workDIR/dummy_${endvar}_${force}_step3_$i.nc

            done

            cdo mergetime $workDIR/dummy_${endvar}_${force}_step3_*.nc $2

            rm $workDIR/dummy_${force}_step2.nc
            rm $workDIR/dummy_${endvar}_${force}_step3_*.nc

    elif [ "$endvar" == "icestart" ]
    then

            if [ "$var" == "lakeicefrac" ]
            then
                    # $1: $workDIR/step1.nc from initial file merge
                    # $2: $workDIR/dummy_${endvar}_${force}_step4.nc

                    # $1: should be $workDIR/step1_${per}.nc
                    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,10/12 -seldate,$((${endper:0:4}+0))-01-01T00:00:00,$((${endper:5:4}+0))-01-01T00:00:00 -sellevidx,1 $1 $workDIR/dummy_${force}_step1_a.nc
                    # select January to September for 2nd to last year (lag by 365)
                    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/9 -seldate,$((${endper:0:4}+1))-01-01T00:00:00,$((${endper:5:4}+0))-12-31T00:00:00 -sellevidx,1 $1 $workDIR/dummy_${force}_step1_b.nc

            elif [ "$var" == "icethick" ]
            then

                    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,10/12 -seldate,$((${endper:0:4}+0))-01-01T00:00:00,$((${endper:5:4}+0))-01-01T00:00:00 $1 $workDIR/dummy_${force}_step1_a.nc
                    # select January to September for 2nd to last year (lag by 365)
                    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/9 -seldate,$((${endper:0:4}+1))-01-01T00:00:00,$((${endper:5:4}+0))-12-31T00:00:00 $1 $workDIR/dummy_${force}_step1_b.nc

            fi

            cdo -O mergetime $workDIR/dummy_${force}_step1_a.nc $workDIR/dummy_${force}_step1_b.nc $workDIR/dummy_${force}_step2.nc

            # ice dur loop (1st two elements of timeranges list)
            for i in $(seq $((${endper:0:4}+0)) $((${endper:5:4}-1))); do

                isleap_oct $i
                cdo -O -L fldmean -subc,$day1 -timmin -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/dummy_${force}_step2.nc $workDIR/dummy_${endvar}_${force}_step3_$i.nc

            done

            cdo mergetime $workDIR/dummy_${endvar}_${force}_step3_*.nc $2

            rm $workDIR/dummy_${force}_step1_a.nc
            rm $workDIR/dummy_${force}_step1_b.nc
            rm $workDIR/dummy_${force}_step2.nc
            rm $workDIR/dummy_${endvar}_${force}_step3_*.nc


    elif [ "$endvar" == "iceend" ]
    then

            if [ "$var" == "lakeicefrac" ]
            then
                    # $1: $workDIR/step1.nc from initial file merge
                    # $2: $workDIR/dummy_${endvar}_${force}_step4.nc

                    # $1: should be $workDIR/step1_${per}.nc
                    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,9/12 -seldate,$((${endper:0:4}+0))-01-01T00:00:00,$((${endper:5:4}+0))-01-01T00:00:00 -sellevidx,1 $1 $workDIR/dummy_${force}_step1_a.nc
                    # select January to September for 2nd to last year (lag by 365)
                    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/8 -seldate,$((${endper:0:4}+1))-01-01T00:00:00,$((${endper:5:4}+0))-12-31T00:00:00 -sellevidx,1 $1 $workDIR/dummy_${force}_step1_b.nc

            elif [ "$var" == "icethick" ]
            then

                    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,9/12 -seldate,$((${endper:0:4}+0))-01-01T00:00:00,$((${endper:5:4}+0))-01-01T00:00:00 $1 $workDIR/dummy_${force}_step1_a.nc
                    # select January to September for 2nd to last year (lag by 365)
                    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/8 -seldate,$((${endper:0:4}+1))-01-01T00:00:00,$((${endper:5:4}+0))-12-31T00:00:00 $1 $workDIR/dummy_${force}_step1_b.nc

            fi

            cdo -O mergetime $workDIR/dummy_${force}_step1_a.nc $workDIR/dummy_${force}_step1_b.nc $workDIR/dummy_${force}_step2.nc

            # ice dur loop (1st two elements of timeranges list)
            for i in $(seq $((${endper:0:4}+0)) $((${endper:5:4}-1))); do

                isleap_sep $i
                cdo -O -L fldmean -subc,$day1 -timmax -seldate,$i-09-01T00:00:00,$(($i+1))-08-31T00:00:00 $workDIR/dummy_${force}_step2.nc $workDIR/dummy_${endvar}_${force}_step3_$i.nc

            done

            cdo mergetime $workDIR/dummy_${endvar}_${force}_step3_*.nc $2

            rm $workDIR/dummy_${force}_step1_a.nc
            rm $workDIR/dummy_${force}_step1_b.nc
            rm $workDIR/dummy_${force}_step2.nc
            rm $workDIR/dummy_${endvar}_${force}_step3_*.nc

    fi
}

# function for setting attributes on processed files
attribute_setter(){

    if [ "$endvar" == "icedur" ]
    then

            cdo -O -L setattribute,${endvar}@standard_name="duration_of_lake_ice_cover" -setattribute,${endvar}@long_name="Duration of lake ice cover" -setname,$endvar -setunit,"days" $1 $2

    elif [ "$endvar" == "icestart" ]
    then

            cdo -O -L setattribute,${endvar}@standard_name="start_of_lake_ice_cover" -setattribute,${endvar}@long_name="Start of lake ice cover" -setname,$endvar -setunit,"days" $1 $2

    elif [ "$endvar" == "iceend" ]
    then

            cdo -O -L setattribute,${endvar}@standard_name="end_of_lake_ice_cover" -setattribute,${endvar}@long_name="End of lake ice cover" -setname,$endvar -setunit,"days" $1 $2
    fi
}

# function for setting time on annual fldmeans
date_setter(){

    cdo -O -L setcalendar,proleptic_gregorian -setreftime,1661-01-01,00:00:00,1years -settaxis,$((${endper:0:4}+0))-01-01,00:00:00,1years $1 $2
}


# ==============================================================================
# PROCESSING
# ==============================================================================


# go into climate directory
cd $inDIR
pwd


for force in "${forcing[@]}"; do
    echo $force

    if [ "$scen_folder" == "pre-industrial" ]
    then
            periods=("1661_1670" "1671_1680" "1681_1690" "1691_1700" "1701_1710" "1711_1720" "1721_1730" "1731_1740" "1741_1750" "1751_1760" "1761_1770" "1771_1780" "1781_1790" "1791_1800" "1801_1810" "1811_1820" "1821_1830" "1831_1840" "1841_1850" "1851_1860")

            endper="1661_1860"

    elif [ "$scen_folder" == "historical" ]
    then
            periods=("1861_1870" "1871_1880" "1881_1890" "1891_1900" "1901_1910" "1911_1920" "1921_1930" "1931_1940" "1941_1950" "1951_1960" "1961_1970" "1971_1980" "1981_1990" "1991_2000" "2001_2005")

            endper="1861_2005"

    elif [ "$scen_folder" == "future" ]
    then
            periods=("2006_2010" "2011_2020" "2021_2030" "2031_2040" "2041_2050" "2051_2060" "2061_2070" "2071_2080" "2081_2090" "2091_2099")

            endper="2006_2099"
    fi

    for opt in "${options[@]}"; do

        for scen in "${scenarios[@]}"; do

            if ls $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_*.nc4 1> /dev/null 2>&1;
            then

                    cdo mergetime $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_*.nc4 $workDIR/step1.nc

            fi

            for endvar in "${endvariables[@]}"; do

                index_calculator $workDIR/step1.nc $workDIR/dummy_${endvar}_${force}_step4.nc

                attribute_setter $workDIR/dummy_${endvar}_${force}_step4.nc $workDIR/dummy_${endvar}_${force}_step5.nc

                ncwa -C -v time,$endvar -a lat,lon,levlak,lev $workDIR/dummy_${endvar}_${force}_step5.nc $workDIR/dummy_${endvar}_${force}_step6.nc

                date_setter $workDIR/dummy_${endvar}_${force}_step6.nc $outDIR/${model_fname}_${force}_${scen}_${endvar}_${prod}_${endper}.nc

                rm $workDIR/dummy_${endvar}_${force}_step4.nc
                rm $workDIR/dummy_${endvar}_${force}_step5.nc
                rm $workDIR/dummy_${endvar}_${force}_step6.nc


            done

            rm $workDIR/step1.nc

        done
    done
done



