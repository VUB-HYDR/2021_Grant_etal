#!/bin/bash -l

# ==============================================================================
# SUMMARY
# ==============================================================================


# Daily/monthly files to 1981-2019 watertemp for evaluation and detection/attribution
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

# set depth directory for indexing fields
deDIR=

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


flag_var=0;     # 0: watertemp
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
forcing=("ipsl-cm5a-lr" "miroc5")


# define timestep
timesteps=("daily" "monthly")


# define all options
options=("2005soc_co2" "1860soc_co2" "nosoc_co2")


# define all end products
products=("fldmean" "sig" "eval")


# set run settings based on flags
tstep=${timesteps[$flag_tstep]}
model=${models[$flag_models]}
model_fname=${model_fnames[$flag_models]}
var=${variables[$flag_var]}
opt=${options[$flag_opt]}
prod=${products[$flag_prod]}


# periods for signal
if [ "$tstep" == "daily" ]
then

    hist_periods=("1981_1990" "1991_2000" "2001_2005")
    fut_periods=("2006_2010" "2011_2020")

elif [ "$tstep" == "monthly" ]
then

    hist_periods=("1861_2005")
    fut_periods=("2006_2099")
    
fi


# ==============================================================================
# FUNCTIONS
# ==============================================================================


lev_indexer(){  # $1:
                # $2:

        if [ "$model_fname" == "clm45" ]
        then

                cdo -O sellevidx,3 $1 $workDIR/intra_step.nc

        elif [ "$model_fname" == "albm" ]
        then

                cdo intlevel3d,$deDIR/albm_og_field.nc $1 $deDIR/albm_index_field.nc $workDIR/intra_step.nc

        elif [ "$model_fname" == "lake" ]
        then

                echo "No system for indexing yet"

        elif [ "$model_fname" == "simstrat-uog" ]
        then

                cdo -O sellevidx,3 $1 $workDIR/intra_step.nc

        elif [ "$model_fname" == "vic-lake" ]
        then

                cdo -O sellevidx,3 $1 $workDIR/intra_step.nc

        elif [ "$model_fname" == "gotm" ]
        then

                echo "No system for indexing yet"

        fi

    ncks -O --mk_rec_dmn time $workDIR/intra_step.nc $2

}

# function for taking mean
mean_calculator(){  # $1: $workDIR/${model_fname}_step2_${rcp}_${force}.nc
                    # $2: $workDIR/${model_fname}_${var}_${force}_${rcp}_step4.nc

    cdo yearmean $1 $2
}

# function for setting attributes on processed files
attribute_setter(){         # $1: $workDIR/${model_fname}_${var}_${force}_${rcp}_step4.nc
                            # $2: $outDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2018.nc

    cdo -O -L setreftime,1661-01-01,00:00:00,1years -settaxis,1981-01-01,00:00:00,1years $1 $2
}

# function for removing levlak and masking
masker(){          # $1: $workDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2018.nc
                   # $2: $outDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2018.nc

    ncwa -C -v time,lat,lon,watertemp -a levlak $1 intra_step.nc
    cdo -O ifthen $maskfile intra_step.nc $2
    rm intra_step.nc
}


# ==============================================================================
# PROCESSING
# ==============================================================================


# go into climate directory (indir up to CLM45)
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

                lev_indexer $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4

            done

        elif [ "$scen_folder" == "future" ]
        then

            for rcp in "${rcps[@]}"; do

                for per in "${fut_periods[@]}"; do

                    lev_indexer $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4

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

    if [ "$tstep" == "daily" ]
    then
    
        cdo mergetime $workDIR/${model_fname}_${force}_ewembi_historical_${opt}_${var}_global_${tstep}_*.nc4 $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc

        for rcp in "${rcps[@]}"; do

            # first merge future
            cdo mergetime $workDIR/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_*.nc4 $workDIR/${model_fname}_step1_${rcp}_${force}_2006_2020.nc

            # merge historical + future
            cdo mergetime $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc $workDIR/${model_fname}_step1_${rcp}_${force}_2006_2020.nc $workDIR/${model_fname}_step2_${rcp}_${force}.nc

            cdo selyear,1981/2019 $workDIR/${model_fname}_step2_${rcp}_${force}.nc $workDIR/${model_fname}_step3_${rcp}_${force}.nc

        done

    elif [ "$tstep" == "monthly" ]
    then

        cdo selyear,1981/2005 $workDIR/${model_fname}_${force}_ewembi_historical_${opt}_${var}_global_${tstep}_1861_2005.nc4 $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc

        for rcp in "${rcps[@]}"; do

            cdo selyear,2006/2019 $workDIR/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_2006_2099.nc4 $workDIR/${model_fname}_step1_${rcp}_${force}_2006_2020.nc

            cdo mergetime $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc $workDIR/${model_fname}_step1_${rcp}_${force}_2006_2020.nc $workDIR/${model_fname}_step2_${rcp}_${force}.nc

        done

    fi
done


echo "clean-up"

for force in "${forcing[@]}"; do

    rm $workDIR/${model_fname}_step1_hist_${force}_1981_2005.nc

    for rcp in "${rcps[@]}"; do

        rm $workDIR/${model_fname}_step1_${rcp}_${force}_2006_2020.nc
        rm $workDIR/${model_fname}_step2_${rcp}_${force}.nc

    done

done


# ==============================================================================
# FINAL PRODUCTS
# ==============================================================================


echo "calculating watertemp means"

for force in "${forcing[@]}"; do
    
    for rcp in "${rcps[@]}"; do
            
        mean_calculator $workDIR/${model_fname}_step3_${rcp}_${force}.nc $workDIR/${model_fname}_${var}_${force}_${rcp}_step4.nc
        attribute_setter $workDIR/${model_fname}_${var}_${force}_${rcp}_step4.nc $workDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2019.nc
        masker $workDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2018.nc $outDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2019.nc
        cdo -O timmean $outDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2018.nc $outDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_timmean_1981_2019.nc
        cdo -O fldmean $outDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_1981_2018.nc $outDIR/${model_fname}_${force}_${rcp}_${var}_${prod}_fldmean_1981_2019.nc

    done
done


# ==============================================================================
# CLEANUP
# ==============================================================================


echo "cleanup"

rm $workDIR/*.nc
rm $workDIR/*.nc4


