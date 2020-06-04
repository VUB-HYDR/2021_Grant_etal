#!/bin/bash -l

# ==============================================================================
# SUMMARY
# ==============================================================================


# Daily picontrol simulated watertemp files to 1981-2019 pichunks


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

# lake mask for ERA5-Land
maskfile=


# Define settings
flag_sect=1;   	# 0: all sectors
				# 1: lakes_global


flag_models=0;  # 0: CLM45
				# 1: ALBM
				# 2: LAKE
				# 3: SIMSTRAT-UoG
				# 4: VIC-LAKE


flag_tstep=1;	# 0: daily
				# 1: monthly
				# 2: annual


flag_var=0;     # 0: watertemp
                # 1: lakeicefrac
                # 2: icethick


flag_opt=0;     # 0: 2005soc_co2
                # 1: 1860soc_co2
                # 2: nosoc_co2


flag_prod=2;    # 0: fldmean
                # 1: sig
                # 2: eval


# define all possible sectors
sectors=("all" "lakes_global")


# define all possible models (top list; folder style, bottom list; file name style)
models=("CLM45" "ALBM" "LAKE" "SIMSTRAT-UoG" "VIC-LAKE")
model_fnames=("clm45" "albm" "lake" "simstrat-uog" "vic-lake")


# scenario folders
scenario_folders=("pre-industrial" "historical" "future")


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


# define all options
products=("fldmean" "sig" "eval")


# set run settings based on flags
tstep=${timesteps[$flag_tstep]}
model=${models[$flag_models]}
model_fname=${model_fnames[$flag_models]}
var=${variables[$flag_var]}
prod=${products[$flag_prod]}


# starting years for pichunks in dictionary
chunk_ids=($(seq 0 1 10))

start=1661;

declare -A chunks

for i in "${chunk_ids[@]}"; do
    
    if [ "$i" == 0 ]
    then
        chunks[$i]=$start
    elif [ "$i" -gt 0 ]
    then
        chunks[$i]=$(( $start+39*$i ))
    fi
done
        
        
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


mean_calculator(){  # $1: $workDIR/${model_fname}_${force}_step1.nc
                    # $2: $workDIR/${model_fname}_${var}_${force}_step2_chunk$i.nc

    y1=${chunks[$i]}
    y2=$(( $y1+39 ))
    
    cdo -L yearmean -seldate,$y1-01-01T00:00:00,$y2-12-31T00:00:00 $1 $workDIR/intra_step1.nc
    cdo ifthen $maskfile $workDIR/intra_step1.nc $2
    
    rm $workDIR/intra_step1.nc
}


# function for setting attributes on processed files
attribute_setter(){         # $1: $workDIR/${model_fname}_${var}_${force}_step2_chunk$i.nc
                            # $2: $workDIR/${model_fname}_${var}_${force}_step3_chunk$i.nc

    cdo -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years $1 $2
}


# ==============================================================================
# PROCESSING
# ==============================================================================


# go into climate directory
cd $inDIR
pwd

# perform whole operation under forcing loop to limit data consumption of intermediary/copied/merged files
for force in "${forcing[@]}"; do


# ==============================================================================
# FILE ALLOCATION
# ==============================================================================


    echo "moving/levlak-selecting relevant files"

    for scen_folder in "${scenario_folders[@]}"; do
    
        for opt in "${options[@]}"; do
        
            if [ "$tstep" == "daily" ]
            then
            
                if [ "$scen_folder" == "pre-industrial" ]
                then
                
                    periods=("1661_1670" "1671_1680" "1681_1690" "1691_1700" "1701_1710" "1711_1720" "1721_1730" "1731_1740" "1741_1750" "1751_1760" "1761_1770" "1771_1780" "1781_1790" "1791_1800" "1801_1810" "1811_1820" "1821_1830" "1831_1840" "1841_1850" "1851_1860")

                elif [ "$scen_folder" == "historical" ]
                then
                
                    periods=("1861_1870" "1871_1880" "1881_1890" "1891_1900" "1901_1910" "1911_1920" "1921_1930" "1931_1940" "1941_1950" "1951_1960" "1961_1970" "1971_1980" "1981_1990" "1991_2000" "2001_2005")


                elif [ "$scen_folder" == "future" ]
                then
                        periods=("2006_2010" "2011_2020" "2021_2030" "2031_2040" "2041_2050" "2051_2060" "2061_2070" "2071_2080" "2081_2090" "2091_2099")

                        endper="2006_2099"
                fi
                
            elif [ "$tstep" == "monthly" ]
            then
            
                if [ "$scen_folder" == "pre-industrial" ]
                then

                    periods=("1661_1860")

                elif [ "$scen_folder" == "historical" ]
                then
                
                    periods=("1861_2005")

                elif [ "$scen_folder" == "future" ]
                then
                
                    periods=("2006_2099")
                
                fi 
            fi
            
            for per in "${periods[@]}"; do
            
                if test -f $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_picontrol_${opt}_${var}_global_${tstep}_${per}.nc4;
                then

                    lev_indexer $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_picontrol_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/${model_fname}_${force}_ewembi_picontrol_${opt}_${var}_global_${tstep}_${per}.nc4
                        
                fi
            done
        done
    done


# ==============================================================================
# MERGING
# ==============================================================================


    echo "merging all transferred picontrol files"

    cdo -O mergetime $workDIR/${model_fname}_${force}_ewembi_picontrol_*_${var}_global_${tstep}_*.nc4 $workDIR/${model_fname}_${force}_step1.nc
        
    echo "clean up"

    rm $workDIR/${model_fname}_${force}_ewembi_picontrol_*_${var}_global_${tstep}_*.nc4


# ==============================================================================
# ICE INDEXES
# ==============================================================================


    echo "calculating ice indexes"


    for i in "${chunk_ids[@]}"; do
    
        mean_calculator $workDIR/${model_fname}_${force}_step1.nc $workDIR/${model_fname}_${var}_${force}_step2_chunk$i.nc
        attribute_setter $workDIR/${model_fname}_${var}_${force}_step2_chunk$i.nc $workDIR/${model_fname}_${var}_${force}_step3_chunk$i.nc
        cdo fldmean $workDIR/${model_fname}_${var}_${force}_step3_chunk$i.nc $outDIR/${model_fname}_${var}_${force}_picontrol_chunk$i.nc

    done



# ==============================================================================
# CLEANUP
# ==============================================================================


    echo "clean up"

    rm $workDIR/*.nc

done
