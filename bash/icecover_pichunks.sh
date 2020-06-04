#!/bin/bash -l

# ==============================================================================
# SUMMARY
# ==============================================================================


# Daily picontrol icethickness files to 1981-2019 ice cover pichunks


# =======================================================================
# INITIALIZATION
# =======================================================================


# set output directory
outDIR=

# set starting directory
inDIR=

# set working directory
workDIR=

# lake mask for ERA5-Land
maskfile=


# Define settings
flag_sect=1;   	# 0: all sectors
				# 1: lakes_global


flag_models=3;  # 0: CLM45
				# 1: ALBM
				# 2: LAKE
				# 3: SIMSTRAT-UoG
				# 4: VIC-LAKE


flag_tstep=0;	# 0: daily
				# 1: monthly
				# 2: annual


flag_var=2;     # 0: watertemp
                # 1: lakeicefrac
                # 2: icethick


flag_opt=2;     # 0: 2005soc_co2
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
opt=${options[$flag_opt]}


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
                

# function for calculating ice start
icestart_index(){      # $1: $workDIR/${model_fname}_${force}_step1.nc 
                       # $2: $workDIR/${model_fname}_${endvar}_${force}_step2_chunk$i.nc
                       
    y1=${chunks[${int}]}
    y2=$(( $y1+39 ))
    
    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,10/12 -seldate,$y1-01-01T00:00:00,$y2-01-01T00:00:00 $1 $workDIR/intra_step1_a.nc
    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/9 -seldate,$(( $y1+1 ))-01-01T00:00:00,$y2-12-31T00:00:00 $1 $workDIR/intra_step1_b.nc
    cdo -O mergetime $workDIR/intra_step1_a.nc $workDIR/intra_step1_b.nc $workDIR/intra_step2.nc

    for i in $(seq $y1 $(( $y2-1 ))); do

        cdo -O -L timmin -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/intra_step2.nc $workDIR/intra_step3_$i.nc

    done

    cdo mergetime $workDIR/intra_step3_*.nc $workDIR/intra_step4.nc
    cdo ifthen $maskfile $workDIR/intra_step4.nc $2

    rm $workDIR/intra_step1_a.nc
    rm $workDIR/intra_step1_b.nc
    rm $workDIR/intra_step2.nc
    rm $workDIR/intra_step3_*.nc
    rm $workDIR/intra_step4.nc
}


# function for calculating ice end
iceend_index(){      # $1: $workDIR/${model_fname}_${force}_step1.nc
                     # $2: $workDIR/${model_fname}_${endvar}_${force}_step2_chunk$i.nc
    
    y1=${chunks[${int}]}
    y2=$(( $y1+39 ))
    
    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,9/12 -seldate,$y1-01-01T00:00:00,$y2-01-01T00:00:00 $1 $workDIR/intra_step1_a.nc
    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/8 -seldate,$(( $y1+1 ))-01-01T00:00:00,$y2-12-31T00:00:00 $1 $workDIR/intra_step1_b.nc
    cdo -O mergetime $workDIR/intra_step1_a.nc $workDIR/intra_step1_b.nc $workDIR/intra_step2.nc
    
    for i in $(seq $y1 $(( $y2-1 ))); do

        cdo -O -L timmax -seldate,$i-09-01T00:00:00,$(($i+1))-08-31T00:00:00 $workDIR/intra_step2.nc $workDIR/intra_step3_$i.nc

    done
    
    cdo mergetime $workDIR/intra_step3_*.nc $workDIR/intra_step4.nc
    cdo ifthen $maskfile $workDIR/intra_step4.nc $2

    rm $workDIR/intra_step1_a.nc
    rm $workDIR/intra_step1_b.nc
    rm $workDIR/intra_step2.nc
    rm $workDIR/intra_step3_*.nc
    rm $workDIR/intra_step4.nc
}


# function for calculating ice end
icedur_index(){      # $1: $workDIR/${model_fname}_${force}_step1.nc
                     # $2: $workDIR/${model_fname}_${endvar}_${force}_step2_chunk${int}.nc

    y1=${chunks[${int}]}
    y2=$(( $y1+39 ))                 
                     
    cdo -O -L setctomiss,0 -gtc,0 -seldate,$y1-01-01T00:00:00,$y2-12-31T00:00:00 $1 $workDIR/step1.nc
                     
    for i in $(seq $y1 $(( $y2-1 ))); do

        cdo -b F64 -O -L timsum -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/step1.nc $workDIR/step2_$i.nc

    done
    
    cdo -O mergetime $workDIR/step2_*.nc $workDIR/step3.nc
    cdo ifthen $maskfile $workDIR/step3.nc $2
    
    rm $workDIR/step1.nc
    rm $workDIR/step2_*.nc
    rm $workDIR/step3.nc
    
}


# function for setting attributes on processed files
attribute_setter(){         # $1: $workDIR/${model_fname}_${endvar}_${force}_step3_chunk${int}.nc
                            # $2: $workDIR/${model_fname}_${endvar}_${force}_step4_chunk${int}.nc

    if [ "$endvar" == "icedur" ]
    then
            cdo -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years -setattribute,iceend@long_name='Duration of lake ice cover (end-start)' -setname,'icedur' -setunit,'days' $1 $2
            
    elif [ "$endvar" == "icestart" ]
    then
            cdo -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years -setattribute,icestart@long_name='First day of lake ice cover' -setname,'icestart' -setunit,'day of hydrological year' $1 $2

    elif [ "$endvar" == "iceend" ]
    then
            cdo -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years -setattribute,iceend@long_name='Last day of lake ice cover' -setname,'iceend' -setunit,'day of hydrological year' $1 $2
    fi
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


    echo "moving relevant files"

    for scen_folder in "${scenario_folders[@]}"; do
    
        for opt in "${options[@]}"; do
            
            if [ "$scen_folder" == "pre-industrial" ]
            then
            
                periods=("1661_1670" "1671_1680" "1681_1690" "1691_1700" "1701_1710" "1711_1720" "1721_1730" "1731_1740" "1741_1750" "1751_1760" "1761_1770" "1771_1780" "1781_1790" "1791_1800" "1801_1810" "1811_1820" "1821_1830" "1831_1840" "1841_1850" "1851_1860")

            elif [ "$scen_folder" == "historical" ]
            then
            
                periods=("1861_1870" "1871_1880" "1881_1890" "1891_1900" "1901_1910" "1911_1920" "1921_1930" "1931_1940" "1941_1950" "1951_1960" "1961_1970" "1971_1980" "1981_1990" "1991_2000" "2001_2005")


            elif [ "$scen_folder" == "future" ]
            then
            
                periods=("2006_2010" "2011_2020" "2021_2030" "2031_2040" "2041_2050" "2051_2060" "2061_2070" "2071_2080" "2081_2090" "2091_2099")

            fi
            
            for per in "${periods[@]}"; do
            
                if test -f $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_picontrol_${opt}_${var}_global_${tstep}_${per}.nc4;
                then

                    if [ "$var" == "lakeicefrac" ]
                    then

                        cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_picontrol_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/${model_fname}_${force}_ewembi_picontrol_${opt}_${var}_global_${tstep}_${per}.nc4
                        

                    elif [ "$var" == "icethick" ]
                    then

                        cp $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_picontrol_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR
                    
                    fi
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

    for endvar in "${endvariables[@]}"; do
    
        for int in "${chunk_ids[@]}"; do
    
            if [ "$endvar" == "icestart" ]
            then
            
                icestart_index $workDIR/${model_fname}_${force}_step1.nc $workDIR/${model_fname}_${endvar}_${force}_step2_chunk${int}.nc
                attribute_setter $workDIR/${model_fname}_${endvar}_${force}_step2_chunk${int}.nc $workDIR/${model_fname}_${endvar}_${force}_step3_chunk${int}.nc
                cdo fldmean $workDIR/${model_fname}_${endvar}_${force}_step3_chunk${int}.nc $outDIR/${model_fname}_${endvar}_${force}_picontrol_chunk${int}.nc


            elif [ "$endvar" == "iceend" ]
            then
            
                iceend_index $workDIR/${model_fname}_${force}_step1.nc $workDIR/${model_fname}_${endvar}_${force}_step2_chunk${int}.nc
                attribute_setter $workDIR/${model_fname}_${endvar}_${force}_step2_chunk${int}.nc $workDIR/${model_fname}_${endvar}_${force}_step3_chunk${int}.nc
                cdo fldmean $workDIR/${model_fname}_${endvar}_${force}_step3_chunk${int}.nc $outDIR/${model_fname}_${endvar}_${force}_picontrol_chunk${int}.nc


            elif [ "$endvar" == "icedur" ]
            then

                icedur_index $workDIR/${model_fname}_${force}_step1.nc $workDIR/${model_fname}_${endvar}_${force}_step2_chunk${int}.nc
                attribute_setter $workDIR/${model_fname}_${endvar}_${force}_step2_chunk${int}.nc $workDIR/${model_fname}_${endvar}_${force}_step5_chunk${int}.nc
                cdo fldmean $workDIR/${model_fname}_${endvar}_${force}_step5_chunk${int}.nc $outDIR/${model_fname}_${endvar}_${force}_picontrol_chunk${int}.nc

            fi
        done
    done


# ==============================================================================
# CLEANUP
# ==============================================================================


    echo "clean up"

    rm $workDIR/*.nc
    rm $workDIR/*.nc4

done
