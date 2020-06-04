#!/bin/bash -l

# =======================================================================
# SUMMARY
# =======================================================================


# Daily to annual field means processing script for simstrat water temperature


# =======================================================================
# INITIALIZATION
# =======================================================================


# set output directory
outDIR=

# set starting directory
inDIR=

# set working directory
workDIR=

# set depth directory for indexing fields
deDIR=


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


flag_var=0;     # 0: watertemp
                # 1: lakeicefrac
                # 2: icethick


flag_endvar=0;  # 0: icestart
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
scenarios=("picontrol" "historical" "rcp26" "rcp60" "rcp85")


# define lake variables
variables=("watertemp" "lakeicefrac" "icethick")


# define forcing
# forcing=("gfdl-esm2m" "hadgem2-es" "ipsl-cm5a-lr" "miroc5")
forcing=("gfdl-esm2m")

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


# ==============================================================================
# FUNCTIONS
# ==============================================================================


lev_indexer(){  # $1:
                # $2:

        if [ "$model_fname" == "clm45" ]
        then

                cdo -O sellevidx,3 $1 $2

        elif [ "$model_fname" == "albm" ]
        then

                cdo intlevel3d,$deDIR/albm_og_field.nc $1 $deDIR/albm_index_field.nc $2

        elif [ "$model_fname" == "lake" ]
        then

                echo "No system for indexing yet"

        elif [ "$model_fname" == "simstrat-uog" ]
        then

                cdo -O sellevidx,3 $1 $2

        elif [ "$model_fname" == "vic-lake" ]
        then

                cdo -O sellevidx,3 $1 $2

        elif [ "$model_fname" == "gotm" ]
        then

                echo "No system for indexing yet"

        fi
}



# ==============================================================================
# PROCESSING
# ==============================================================================


# go into climate directory (indir up to CLM45)
cd $inDIR
pwd


# go into climate directory
cd $inDIR
pwd

for force in "${forcing[@]}"; do
    echo $force

    for scen_folder in "${scenario_folders[@]}"; do

        if [ "$tstep" == "daily" ]
        then
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

        elif [ "$tstep" == "monthly" ]
        then
                if [ "$scen_folder" == "pre-industrial" ]
                then
                        endper="1661_1860"

                elif [ "$scen_folder" == "historical" ]
                then
                        endper="1861_2005"

                elif [ "$scen_folder" == "future" ]
                then
                        endper="2006_2099"
                fi
        fi

        for opt in "${options[@]}"; do

            for scen in "${scenarios[@]}"; do


                if [ "$tstep" == "monthly" ]
                then
                        if test -f $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${endper}.nc4;
                        then

                                cdo yearmean $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${endper}.nc4 $workDIR/step1.nc

                                lev_indexer $workDIR/step1.nc $workDIR/step1_a.nc

                                cdo fldmean $workDIR/step1_a.nc $workDIR/step2.nc
                        fi


                elif [ "$tstep" == "daily" ]
                then

                        for per in "${periods[@]}"; do

                            if test -f $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4;
                            then

                                    cdo -L yearmean $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step1_${per}.nc

                                    lev_indexer $workDIR/step1_${per}.nc $workDIR/step1_a_${per}.nc

                                    cdo fldmean $workDIR/step1_a_${per}.nc $workDIR/step1_b_${per}.nc
                            fi
                        done


                        cdo mergetime $workDIR/step1_b_*.nc $workDIR/step2.nc

                        rm $workDIR/step1*.nc

                fi

                ncwa -C -v time,$var -a levlak,lat,lon,lev $workDIR/step2.nc $workDIR/step3.nc

                cdo -O -L setcalendar,proleptic_gregorian -setreftime,1661-01-01,00:00:00,1years -settaxis,$((${endper:0:4}+0))-01-01,00:00:00,1years $workDIR/step3.nc $outDIR/${model_fname}_${force}_${scen}_${var}_${prod}_${endper}.nc

                rm $workDIR/step2.nc
                rm $workDIR/step3.nc


            done
        done
    done
done


