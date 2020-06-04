#!/bin/bash -l

# =======================================================================
# SUMMARY
# =======================================================================


# Ice thickness files to fldmeans


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


flag_scen=2;    # 0: historical
                # 1: pre-industrial
                # 2: future


flag_tstep=0;	# 0: daily
				# 1: monthly
				# 2: annual


flag_var=2;     # 0: watertemp
                # 1: lakeicefrac
                # 2: icethick


flag_endvar=3;  # 0: icestart
                # 1: iceend
                # 2: icedur
                # 3: icethick (clm45 conversions)


flag_opt=2;     # 0: 2005soc_co2
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


#define all RCPs for future scenario
rcps=("rcp26" "rcp60" "rcp85")


# define lake variables
variables=("watertemp" "lakeicefrac" "icethick")


#define end lake variables
endvariables=("icestart" "iceend" "icedur" "icethick")


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
endvar=${endvariables[$flag_endvar]}
prod=${products[$flag_prod]}


# thickness per levidx in CLM45
clm_levthick=(0.1 1.0 2.0 3.0 4.0 5.0 7.0 7.0 10.45 10.45)


# ==============================================================================
# FUNCTIONS
# ==============================================================================


# function for setting time on annual fldmeans
date_setter(){

    cdo -O -L setcalendar,proleptic_gregorian -setreftime,1661-01-01,00:00:00,1years -settaxis,$((${endper:0:4}+0))-01-01,00:00:00,1years $1 $2
}


# ==============================================================================
# PROCESSING
# ==============================================================================


# set working directory
mkdir -p /theia/scratch/projects/climate/users/lgrant/isimip/fldmean_proc/$var/$model_fname

workDIR=/theia/scratch/projects/climate/users/lgrant/isimip/fldmean_proc/$var/$model_fname

# go into climate directory
cd $inDIR
pwd

for force in "${forcing[@]}"; do
    echo $force

        for scen_folder in "${scenario_folders[@]}"; do

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

                for per in "${periods[@]}"; do

                    if test -f $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4;
                    then
                            # time/spatial means
                            cdo -L yearmean -fldmean $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step1_${per}.nc
                    fi
                done


                cdo mergetime $workDIR/step1_*.nc $workDIR/step2.nc
                rm $workDIR/step1*.nc


                if [ "$model_fname" == "clm45" ]
                then

                        # loop through levels, multiply and ensum for icefrac -> ice thick
                        for i in $(seq 0 9); do

                            cdo sellevidx,$((i+1)) $workDIR/step2.nc $workDIR/step3_${i}.nc
                            cdo -L setattribute,${endvar}@long_name="Ice thickness" -setname,${endvar} -setunit,"m" -mulc,${clm_levthick[${i}]} $workDIR/step3_${i}.nc $workDIR/step4_${i}.nc
                            ncwa -C -v time,icethick -a levlak,lat,lon,lev $workDIR/step4_${i}.nc $workDIR/step5_${i}.nc

                            rm $workDIR/step3_${i}.nc
                            rm $workDIR/step4_${i}.nc

                        done

                        cdo enssum $workDIR/step5_*.nc $workDIR/step6.nc

                        date_setter $workDIR/step6.nc $outDIR/${model_fname}_${force}_${scen}_${endvar}_${prod}_${endper}.nc

                        rm $workDIR/step2.nc
                        rm $workDIR/step5_*.nc
                        rm $workDIR/step6.nc


                elif [ "$model_fname" != "clm45" ]
                then

                        ncwa -C -v time,icethick -a lat,lon,levlak,lev $workDIR/step2.nc $workDIR/step3.nc

                        date_setter $workDIR/step3.nc $outDIR/${model_fname}_${force}_${scen}_${endvar}_${prod}_${endper}.nc

                        rm $workDIR/step2.nc
                        rm $workDIR/step3.nc
                fi


            done
        done
    done
done


