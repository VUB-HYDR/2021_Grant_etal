#!/bin/bash -l

# ==============================================================================
# SUMMARY
# ==============================================================================


# Daily ice thickness files to 1981-2019 ice cover series for evaluation and detection/attribution
# EXT response patterns


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
scenario_folders=("historical" "future")


#define all RCPs for future scenario
rcps=("rcp60" "rcp85")


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


# periods for signal
hist_periods=("1981_1990" "1991_2000" "2001_2005")
fut_periods=("2006_2010" "2011_2020")


# ==============================================================================
# FUNCTIONS
# ==============================================================================


# function for calculating ice start
icestart_index(){      # $1: $workDIR/${model_fname}_step2_rcp60_${force}.nc | $workDIR/${model_fname}_step2_rcp85_${force}.nc
                       # $2: $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc
                       
    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,10/12 -seldate,1981-01-01T00:00:00,2019-01-01T00:00:00 $1 $workDIR/step1_a.nc
    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/9 -seldate,1982-01-01T00:00:00,2019-12-31T00:00:00 $1 $workDIR/step1_b.nc
    cdo -O mergetime $workDIR/step1_a.nc $workDIR/step1_b.nc $workDIR/step2.nc

    for i in $(seq 1981 2018); do

        cdo -O -L timmin -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/step2.nc $workDIR/step3_$i.nc

    done

    cdo mergetime $workDIR/step3_*.nc $workDIR/step4.nc
    cdo ifthen $maskfile $workDIR/step4.nc $2

    rm $workDIR/step1_a.nc
    rm $workDIR/step1_b.nc
    rm $workDIR/step2.nc
    rm $workDIR/step3_*.nc
    rm $workDIR/step4.nc
}

# function for calculating ice end
iceend_index(){      # $1: $workDIR/${model_fname}_step2_rcp60_${force}.nc | $workDIR/${model_fname}_step2_rcp85_${force}.nc
                     # $2: $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc
                       
    cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,9/12 -seldate,1981-01-01T00:00:00,2019-01-01T00:00:00 $1 $workDIR/step1_a.nc
    cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/8 -seldate,1982-01-01T00:00:00,2019-12-31T00:00:00 $1 $workDIR/step1_b.nc
    cdo -O mergetime $workDIR/step1_a.nc $workDIR/step1_b.nc $workDIR/step2.nc
    
    for i in $(seq 1981 2018); do

        cdo -O -L timmax -seldate,$i-09-01T00:00:00,$(($i+1))-08-31T00:00:00 $workDIR/step2.nc $workDIR/step3_$i.nc

    done
    
    cdo mergetime $workDIR/step3_*.nc $workDIR/step4.nc
    cdo ifthen $maskfile $workDIR/step4.nc $2

    rm $workDIR/step1_a.nc
    rm $workDIR/step1_b.nc
    rm $workDIR/step2.nc
    rm $workDIR/step3_*.nc
    rm $workDIR/step4.nc
}

# function for calculating ice end
icedur_index(){      # $1: $workDIR/${model_fname}_step2_rcp60_${force}.nc | $workDIR/${model_fname}_step2_rcp85_${force}.nc
                     # $2: $workDIR/${model_fname}_icedur_${force}_${rcp}_step4.nc

    cdo -O -L setctomiss,0 -gtc,0 -seldate,1981-01-01T00:00:00,2019-12-31T00:00:00 $1 $workDIR/step1.nc
                     
    for i in $(seq 1981 2018); do

        cdo -O -L timsum -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/step1.nc $workDIR/step2_$i.nc

    done
    
    cdo -O mergetime $workDIR/step2_*.nc $workDIR/step3.nc
    cdo ifthen $maskfile $workDIR/step3.nc $2
    
    rm $workDIR/step1.nc
    rm $workDIR/step2_*.nc
    rm $workDIR/step3.nc
    
}

# function for setting attributes on processed files
attribute_setter(){         # $1: $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc
                            # $2: $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc

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


# ==============================================================================
# FILE ALLOCATION
# ==============================================================================


echo "moving relevant files"

for force in "${forcing[@]}"; do

    for scen_folder in "${scenario_folders[@]}"; do

        if [ "$scen_folder" == "historical" ]
        then scen="historical"

            for per in "${hist_periods[@]}"; do

                if [ "$var" == "lakeicefrac" ]
                then

                    cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4

                elif [ "$var" == "icethick" ]
                then

                    cp $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR

                fi
            done

        elif [ "$scen_folder" == "future" ]
        then

            for rcp in "${rcps[@]}"; do

                for per in "${fut_periods[@]}"; do

                    if [ "$var" == "lakeicefrac" ]
                    then

                        cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4

                    elif [ "$var" == "icethick" ]
                    then

                        cp $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR

                    fi
                done
            done
        fi
    done
done


# ==============================================================================
# MERGING
# ==============================================================================


echo "merging hist -> rcp60 and hist -> rcp85"

for force in "${forcing[@]}"; do

        cdo mergetime $workDIR/${model_fname}_${force}_ewembi_historical_${opt}_${var}_global_${tstep}_*.nc4 $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc
        cdo mergetime $workDIR/${model_fname}_${force}_ewembi_rcp60_${opt}_${var}_global_${tstep}_*.nc4 $workDIR/${model_fname}_step1_rcp60_${force}_2006_2020.nc
        cdo mergetime $workDIR/${model_fname}_${force}_ewembi_rcp85_${opt}_${var}_global_${tstep}_*.nc4 $workDIR/${model_fname}_step1_rcp85_${force}_2006_2020.nc
        
        cdo mergetime $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc $workDIR/${model_fname}_step1_rcp60_${force}_2006_2020.nc $workDIR/${model_fname}_step2_rcp60_${force}.nc
        cdo mergetime $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc $workDIR/${model_fname}_step1_rcp85_${force}_2006_2020.nc $workDIR/${model_fname}_step2_rcp85_${force}.nc

done


echo "clean up"

for force in "${forcing[@]}"; do

    rm $workDIR/${model_fname}_${force}_ewembi_historical_histsoc_co2_${var}_global_${tstep}_*.nc4
    rm $workDIR/${model_fname}_${force}_ewembi_rcp60_2005soc_co2_${var}_global_${tstep}_*.nc4
    rm $workDIR/${model_fname}_${force}_ewembi_rcp85_2005soc_co2_${var}_global_${tstep}_*.nc4

    rm $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc
    rm $workDIR/${model_fname}_step1_rcp60_${force}_2006_2020.nc
    rm $workDIR/${model_fname}_step1_rcp85_${force}_2006_2020.nc
    
done


# ==============================================================================
# ICE INDEXES
# ==============================================================================


echo "calculating ice indexes"

for force in "${forcing[@]}"; do

    for endvar in "${endvariables[@]}"; do
    
        for rcp in "${rcps[@]}"; do
    
            if [ "$endvar" == "icestart" ]
            then
            
                icestart_index $workDIR/${model_fname}_step2_${rcp}_${force}.nc $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc
                attribute_setter $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc
                cdo timmean $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_timmean_1981_2019.nc
                cdo fldmean $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_fldmean_1981_2019.nc


            elif [ "$endvar" == "iceend" ]
            then
            
                iceend_index $workDIR/${model_fname}_step2_${rcp}_${force}.nc $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc
                attribute_setter $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc
                cdo timmean $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_timmean_1981_2019.nc
                cdo fldmean $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_fldmean_1981_2019.nc


            elif [ "$endvar" == "icedur" ]
            then
            
                icedur_index $workDIR/${model_fname}_step2_${rcp}_${force}.nc $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc
                attribute_setter $workDIR/${model_fname}_${endvar}_${force}_${rcp}_step4.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc
                cdo timmean $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_timmean_1981_2019.nc
                cdo fldmean $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_1981_2019.nc $outDIR/${model_fname}_${force}_${rcp}_${endvar}_${prod}_fldmean_1981_2019.nc

            fi
        done
    done
done


# ==============================================================================
# CLEANUP
# ==============================================================================


echo "clean up"

rm $workDIR/*.nc
rm $workDIR/*.nc4

