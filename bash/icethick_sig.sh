#!/bin/bash -l

# ==============================================================================
# SUMMARY
# ==============================================================================


# Ice thickness files to 1971-2000 to 2070-2099 signals and scaled signals


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
                # 3: icethick (for clm45 conversions)


flag_opt=2;     # 0: 2005soc_co2
                # 1: 1860soc_co2
                # 2: nosoc_co2


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
endvariables=("icestart" "iceend" "icedur" "icethick")


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


# scaling arrays (RCP2.6 RCP6.0 RCP8.5 global mean air temperature changes per GCM)
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

# thickness per levidx in CLM45
clm_levthick=(0.1 1.0 2.0 3.0 4.0 5.0 7.0 7.0 10.45 10.45)


# ==============================================================================
# PROCESSING
# ==============================================================================


# go into climate directory (indir up to CLM45)
cd $inDIR
pwd


# ==============================================================================
# HISTORICAL/FUTURE AVERAGES
# ==============================================================================


echo "calculating hist/future icethick for seasons"

for force in "${forcing[@]}"; do
    echo $force

    for scen_folder in "${scenario_folders[@]}"; do

        if [ "$scen_folder" == "historical" ]
        then scen="historical"

                for per in "${hist_periods[@]}"; do

                    if test -f $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4;
                    then
                            # time/spatial means
                            cdo monmean $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step1_${per}.nc4
                    fi
                done


                cdo mergetime $workDIR/step1_*.nc4 $workDIR/step2.nc4

                rm $workDIR/step1_*.nc4


                if [ "$model_fname" == "clm45" ]
                then

                        # loop through levels, multiply and ensum for icefrac -> ice thick
                        for i in $(seq 0 9); do

                            cdo sellevidx,$((i+1)) $workDIR/step2.nc4 $workDIR/step3_${i}.nc4
                            cdo -L setattribute,${endvar}@standard_name="ice_thickness" -setattribute,${endvar}@long_name="Ice thickness" -setname,${endvar} -setunit,"m" -mulc,${clm_levthick[${i}]} $workDIR/step3_${i}.nc4 $workDIR/step4_${i}.nc4
                            ncwa -C -v time,icethick,lat,lon -a levlak $workDIR/step4_${i}.nc4 $workDIR/step5_${i}.nc4

                            rm $workDIR/step3_${i}.nc4
                            rm $workDIR/step4_${i}.nc4

                        done

                        cdo enssum $workDIR/step5_*.nc4 $workDIR/step6.nc4

                        rm $workDIR/step5_*.nc4

                elif [ "$model_fname" != "clm45" ] #other models already with icethick, but this ensures attribute unity and brings non-clm45 data up to step6
                then

                        cdo -L setattribute,${endvar}@standard_name="ice_thickness" -setattribute,${endvar}@long_name="Ice thickness" -setname,${endvar} -setunit,"m" $workDIR/step2.nc4 $workDIR/step6.nc4


                fi

                # historical seasons (identified later through no "${rcp}")
                cdo -L timmean -selmon,12,1,2 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/DJF_dummy_${force}_step7.nc4
                cdo -L timmean -selmon,3/5 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/MAM_dummy_${force}_step7.nc4
                cdo -L timmean -selmon,6/8 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/JJA_dummy_${force}_step7.nc4
                cdo -L timmean -selmon,9/11 -seldate,1971-01-01T00:00:00,2000-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/SON_dummy_${force}_step7.nc4

                rm $workDIR/step2.nc4
                rm $workDIR/step6.nc4



        elif [ "$scen_folder" == "future" ]
        then

                for rcp in "${rcps[@]}"; do

                    for per in "${fut_periods[@]}"; do

                        if test -f $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4
                        then
                                # time/spatial means
                                cdo monmean $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${per}.nc4 $workDIR/step1_${per}.nc4
                        fi
                    done

                    cdo mergetime $workDIR/step1_*.nc4 $workDIR/step2.nc4

                    rm $workDIR/step1*.nc4


                    if [ "$model_fname" == "clm45" ]
                    then

                            # loop through levels, multiply and ensum for icefrac -> ice thick
                            for i in $(seq 0 9); do

                                cdo sellevidx,$((i+1)) $workDIR/step2.nc4 $workDIR/step3_${i}.nc4
                                cdo -L setattribute,${endvar}@standard_name="ice_thickness" -setattribute,${endvar}@long_name="Ice thickness" -setname,${endvar} -setunit,"m" -mulc,${clm_levthick[${i}]} $workDIR/step3_${i}.nc4 $workDIR/step4_${i}.nc4
                                ncwa -C -v time,icethick,lat,lon -a levlak $workDIR/step4_${i}.nc4 $workDIR/step5_${i}.nc4

                                rm $workDIR/step3_${i}.nc4
                                rm $workDIR/step4_${i}.nc4

                            done

                            cdo enssum $workDIR/step5_*.nc4 $workDIR/step6.nc4

                            rm $workDIR/step5_*.nc4

                    elif [ "$model_fname" != "clm45" ]
                    then

                            cdo -L setattribute,${endvar}@standard_name="ice_thickness" -setattribute,${endvar}@long_name="Ice thickness" -setname,${endvar} -setunit,"m" $workDIR/step2.nc4 $workDIR/step6.nc4


                    fi

                    # future seasons
                    cdo -L timmean -selmon,12,1,2 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/DJF_dummy_${force}_${rcp}_step7.nc4
                    cdo -L timmean -selmon,3/5 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/MAM_dummy_${force}_${rcp}_step7.nc4
                    cdo -L timmean -selmon,6/8 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/JJA_dummy_${force}_${rcp}_step7.nc4
                    cdo -L timmean -selmon,9/11 -seldate,2070-01-01T00:00:00,2099-12-31T00:00:00 $workDIR/step6.nc4 $workDIR/SON_dummy_${force}_${rcp}_step7.nc4

                    rm $workDIR/step2.nc4
                    rm $workDIR/step6.nc4

                done

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

            # subtract for regular icethick signal
            cdo sub $workDIR/${season}_dummy_${force}_${rcp}_step7.nc4 $workDIR/${season}_dummy_${force}_step7.nc4 $workDIR/${season}_dummy_${force}_${rcp}_step8.nc4

            # divide for scaled icethick signal
            cdo divc,$scaler $workDIR/${season}_dummy_${force}_${rcp}_step8.nc4 $workDIR/${season}_dummy_${force}_${rcp}_step9.nc4

            rm $workDIR/${season}_dummy_${force}_${rcp}_step7.nc4

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

        # ensemble mean of regular icethick signals
        cdo -O ensmean $workDIR/${season}_dummy_*_${rcp}_step8.nc4 $workDIR/${model_fname}_${rcp}_${endvar}_${season}_${prod}.nc4

        ncwa -C -v icethick,lat,lon -a levlak,lev,time,time_bnds $workDIR/${model_fname}_${rcp}_${endvar}_${season}_${prod}.nc4 $outDIR/${model_fname}_${rcp}_${endvar}_${season}_${prod}.nc4

        # ensemble mean of scaled icethick signals
        cdo -O ensmean $workDIR/${season}_dummy_*_${rcp}_step9.nc4 $workDIR/${model_fname}_${rcp}_${endvar}_${season}_scaled_${prod}.nc4

        ncwa -C -v icethick,lat,lon -a levlak,lev,time,time_bnds $workDIR/${model_fname}_${rcp}_${endvar}_${season}_scaled_${prod}.nc4 $outDIR/${model_fname}_${rcp}_${endvar}_${season}_scaled_${prod}.nc4


        rm $workDIR/${season}_dummy_*_${rcp}_step8.nc4
        rm $workDIR/${season}_dummy_*_${rcp}_step9.nc4
        rm $workDIR/${model_fname}_${rcp}_${endvar}_${season}_${prod}.nc4
        rm $workDIR/${model_fname}_${rcp}_${endvar}_${season}_scaled_${prod}.nc4

    done
done


# ==============================================================================
# CLEANUP
# ==============================================================================


echo "cleanup"

for force in "${forcing[@]}"; do

    for season in "${seasons[@]}"; do

        rm $workDIR/${season}_dummy_${force}_step7.nc4

    done
done


