#!/bin/bash -l

# ==============================================================================
# SUMMARY
# ==============================================================================


# Daily files to 1971-2000 to 2070-2099 signals and scaled signals for water temperature


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


flag_models=1;  # 0: CLM45
				# 1: ALBM
				# 2: LAKE
				# 3: SIMSTRAT-UoG
				# 4: VIC-LAKE


flag_scen=1;    # 0: pre-industrial
                # 1: historical
                # 2: future


flag_tstep=1;	# 0: daily
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
                # 2: nosoc_co2 (for simstrat hist/future)


flag_prod=1;    # 0: fldmean
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
rcps=("rcp26" "rcp60" "rcp85")


# define lake variables
variables=("watertemp" "lakeicefrac" "icethick")


#define end lake variables
endvariables=("icestart" "iceend" "icedur")


# define forcing
forcing=("gfdl-esm2m" "hadgem2-es" "ipsl-cm5a-lr" "miroc5")


# define timestep
timesteps=("daily" "monthly")


# define all options
options=("2005soc_co2" "1860soc_co2" "nosoc_co2")


# seasons
seasons=("DJF" "MAM" "JJA" "SON")


# define end products
products=("fldmean" "sig" "eval")


# scaling arrays (RCP 2.6 6.0 8.5 global mean air temperature changes per GCM)
gfdl_scalers=(0.76875 1.7713 2.76895)
hadgem2_scalers=(1.52805 3.24965 4.81365)
ipsl_scalers=(1.43065 2.7116 4.60615)
miroc5_scalers=(1.20695 2.15395 3.47205)


# set run settings based on flags
tstep=${timesteps[$flag_tstep]}
model=${models[$flag_models]}
model_fname=${model_fnames[$flag_models]}
var=${variables[$flag_var]}
endvar=${endvariables[$flag_endvar]}
opt=${options[$flag_opt]}
prod=${products[$flag_prod]}


# periods for signal
hist_periods=("1971_1980" "1981_1990" "1991_2000")
fut_periods=("2061_2070" "2071_2080" "2081_2090" "2091_2099")  


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

                cdo -O intlevel3d,$deDIR/albm_og_field.nc $1 $deDIR/albm_index_field.nc $2

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


# go into climate directory
cd $inDIR
pwd


# ==============================================================================
# HISTORICAL/FUTURE AVERAGES
# ==============================================================================


echo "calculating hist/future watertemp for seasons"

for force in "${forcing[@]}"; do
    echo $force

    for scen_folder in "${scenario_folders[@]}"; do

        if [ "$tstep" == "daily" ]
        then

                if [ "$scen_folder" == "historical" ]
                then scen="historical"

                        for per in "${hist_periods[@]}"; do

                            cdo monmean $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step_a_$per.nc

                            lev_indexer $workDIR/step_a_$per.nc $workDIR/step_b_$per.nc

                        done

                        cdo mergetime $workDIR/step_b_*.nc $workDIR/step_b.nc
                        rm $workDIR/step_a_*.nc

                        # historical seasons
                        cdo -L timmean -selmon,12,1,2 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_b.nc $workDIR/DJF_dummy_${force}_step1.nc
                        cdo -L timmean -selmon,3/5 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_b.nc $workDIR/MAM_dummy_${force}_step1.nc
                        cdo -L timmean -selmon,6/8 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_b.nc $workDIR/JJA_dummy_${force}_step1.nc
                        cdo -L timmean -selmon,9/11 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_b.nc $workDIR/SON_dummy_${force}_step1.nc

                        rm $workDIR/step_b.nc

                elif [ "$scen_folder" == "future" ]
                then

                        for rcp in "${rcps[@]}"; do

                            for per in "${fut_periods[@]}"; do

                                cdo monmean $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step_a_$per.nc

                                lev_indexer $workDIR/step_a_$per.nc $workDIR/step_b_$per.nc

                            done

                            cdo mergetime $workDIR/step_b_*.nc $workDIR/step_b.nc
                            rm $workDIR/step_a_*.nc

                            # future seasons
                            cdo -L timmean -selmon,12,1,2 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_b.nc $workDIR/DJF_dummy_${force}_${rcp}_step1.nc
                            cdo -L timmean -selmon,3/5 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_b.nc $workDIR/MAM_dummy_${force}_${rcp}_step1.nc
                            cdo -L timmean -selmon,6/8 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_b.nc $workDIR/JJA_dummy_${force}_${rcp}_step1.nc
                            cdo -L timmean -selmon,9/11 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_b.nc $workDIR/SON_dummy_${force}_${rcp}_step1.nc

                            rm $workDIR/step_b.nc

                        done

                fi


        elif [ "$tstep" == "monthly" ]
        then

                if [ "$scen_folder" == "historical" ]
                then scen="historical"
                     per="1861_2005"

                        lev_indexer $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step_a_${per}.nc

                        # historical seasons
                        cdo -L timmean -selmon,12,1,2 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_a_${per}.nc $workDIR/DJF_dummy_${force}_step1.nc
                        cdo -L timmean -selmon,3/5 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_a_${per}.nc $workDIR/MAM_dummy_${force}_step1.nc
                        cdo -L timmean -selmon,6/8 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_a_${per}.nc $workDIR/JJA_dummy_${force}_step1.nc
                        cdo -L timmean -selmon,9/11 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step_a_${per}.nc $workDIR/SON_dummy_${force}_step1.nc

                elif [ "$scen_folder" == "future" ]
                then per="2006_2099"

                    for rcp in "${rcps[@]}"; do

                            lev_indexer $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step_a_${rcp}_${per}.nc

                            # future seasons
                            cdo -L timmean -selmon,12,1,2 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_a_${rcp}_${per}.nc $workDIR/DJF_dummy_${force}_${rcp}_step1.nc
                            cdo -L timmean -selmon,3/5 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_a_${rcp}_${per}.nc $workDIR/MAM_dummy_${force}_${rcp}_step1.nc
                            cdo -L timmean -selmon,6/8 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_a_${rcp}_${per}.nc $workDIR/JJA_dummy_${force}_${rcp}_step1.nc
                            cdo -L timmean -selmon,9/11 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step_a_${rcp}_${per}.nc $workDIR/SON_dummy_${force}_${rcp}_step1.nc

                    done
                fi

        fi

    done
done


# ==============================================================================
# SIGNALS
# ==============================================================================


echo "regular and scaled signals"

for rcp in "${rcps[@]}"; do

    if [ "$rcp" == "rcp26" ]; then  flag=0
    elif [ "$rcp" == "rcp60" ]; then  flag=1
    elif [ "$rcp" == "rcp85" ]; then  flag=2
    fi
    echo $flag
    echo $rcp

    for force in "${forcing[@]}"; do

        if [ "$force" == "gfdl-esm2m" ]; then   scaler=${gfdl_scalers[$flag]}
        elif [ "$force" == "hadgem2-es" ]; then   scaler=${hadgem2_scalers[$flag]}
        elif [ "$force" == "ipsl-cm5a-lr" ]; then   scaler=${ipsl_scalers[$flag]}
        elif [ "$force" == "miroc5" ]; then   scaler=${miroc5_scalers[$flag]}
        fi
        echo $force

        for season in "${seasons[@]}"; do

            # subtract for regular watertemp signal
            cdo sub $workDIR/${season}_dummy_${force}_${rcp}_step1.nc $workDIR/${season}_dummy_${force}_step1.nc $workDIR/${season}_dummy_${force}_${rcp}_step2.nc

            # divide for scaled watertemp signal
            cdo divc,$scaler $workDIR/${season}_dummy_${force}_${rcp}_step2.nc $workDIR/${season}_dummy_${force}_${rcp}_step3.nc

            rm $workDIR/${season}_dummy_${force}_${rcp}_step1.nc

        done
    done
done


# ==============================================================================
# ENSEMBLE MEANS OF SIGNALS
# ==============================================================================


echo "ensemble means"

for rcp in "${rcps[@]}"; do
echo $rcp

    for season in "${seasons[@]}"; do

        # ensemble mean of regular watertemp signals
        cdo -O ensmean $workDIR/${season}_dummy_*_${rcp}_step2.nc $workDIR/${model_fname}_${rcp}_${var}_${season}_${prod}.nc

        ncwa -C -v lat,lon,$var -a time,levlak,lev $workDIR/${model_fname}_${rcp}_${var}_${season}_${prod}.nc $outDIR/${model_fname}_${rcp}_${var}_${season}_${prod}.nc

        # ensemble mean of scaled watertemp signals
        cdo -O ensmean $workDIR/${season}_dummy_*_${rcp}_step3.nc $workDIR/${model_fname}_${rcp}_${var}_${season}_scaled_${prod}.nc

        ncwa -C -v lat,lon,$var -a time,levlak,lev $workDIR/${model_fname}_${rcp}_${var}_${season}_scaled_${prod}.nc $outDIR/${model_fname}_${rcp}_${var}_${season}_scaled_${prod}.nc

        rm $workDIR/${season}_dummy_*_${rcp}_step2.nc
        rm $workDIR/${season}_dummy_*_${rcp}_step3.nc

    done
done


# ==============================================================================
# CLEANUP
# ==============================================================================


echo "cleanup"

rm $workDIR/*.nc


