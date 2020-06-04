#!/bin/bash -l

# ==============================================================================
# SUMMARY
# ==============================================================================


# Ice thickness files to 1971-2000 to 2070-2099 signals and scaled signals for ice
# cover indices


# =======================================================================
# INITIALIZATION
# =======================================================================


# set output directory
outDIR=


# set starting directory
inDIR=


# set working directory
workDIR=


# Define settings
flag_sect=1;   	# 0: all sectors
				# 1: lakes_global


flag_models=2;  # 0: CLM45
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
opt=${options[$flag_opt]}
prod=${products[$flag_prod]}


# periods for signal
hist_periods=("1971_1980" "1981_1990" "1991_2000")
fut_periods=("2061_2070" "2071_2080" "2081_2090" "2091_2099")   



# ==============================================================================
# FUNCTIONS
# ==============================================================================


# function for calculating ice indexes in the historical period
histindex(){            # $1: $workDIR/dummy_${force}_step1.nc
                        # $2: $workDIR/dummy_${endvar}_${force}_step5.nc

        if [ "$endvar" == "icestart" ]
        then

                # select October to December for 1st to last year (select top level - sellevidx,1 - to reduce RAM load assuming it encompasses ice cover dur)
                cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,10/12 -seldate,$((${hist_periods[0]:0:4}+0))-01-01T00:00:00,$((${hist_periods[2]:5:4}+0))-01-01T00:00:00 $1 $workDIR/dummy_${force}_step1_a.nc
                # select January to September for 2nd to last year (lag by 365)
                cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/9 -seldate,$((${hist_periods[0]:0:4}+1))-01-01T00:00:00,$((${hist_periods[2]:5:4}+0))-12-31T00:00:00 $1 $workDIR/dummy_${force}_step1_b.nc
                cdo -O mergetime $workDIR/dummy_${force}_step1_a.nc $workDIR/dummy_${force}_step1_b.nc $workDIR/dummy_${force}_step2.nc
                # ice dur loop (1st two elements of timeranges list)
                for i in $(seq $((${hist_sig[0]:0:4}+0)) $((${hist_sig[0]:5:4}-1))); do

                    cdo -O -L timmin -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/dummy_${force}_step2.nc $workDIR/dummy_${endvar}_${force}_step3_$i.nc

                done



        elif [ "$endvar" == "iceend" ]
        then

                # select October to December for 1st to last year (select top level - sellevidx,1 - to reduce RAM load assuming it encompasses ice cover dur)
                cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,9/12 -seldate,$((${hist_periods[0]:0:4}+0))-01-01T00:00:00,$((${hist_periods[2]:5:4}+0))-01-01T00:00:00 $1 $workDIR/dummy_${force}_step1_a.nc
                # select January to September for 2nd to last year (lag by 365)
                cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/8 -seldate,$((${hist_periods[0]:0:4}+1))-01-01T00:00:00,$((${hist_periods[2]:5:4}+0))-12-31T00:00:00 $1 $workDIR/dummy_${force}_step1_b.nc
                cdo -O mergetime $workDIR/dummy_${force}_step1_a.nc $workDIR/dummy_${force}_step1_b.nc $workDIR/dummy_${force}_step2.nc
                # ice dur loop (1st two elements of timeranges list)
                for i in $(seq $((${hist_sig[0]:0:4}+0)) $((${hist_sig[0]:5:4}-1))); do

                    cdo -O -L timmax -seldate,$i-09-01T00:00:00,$(($i+1))-08-31T00:00:00 $workDIR/dummy_${force}_step2.nc $workDIR/dummy_${endvar}_${force}_step3_$i.nc

                done

        elif [ "$endvar" == "icedur" ]
        then

                # set aside ice days
                cdo -O -L setctomiss,0 -gtc,0 $1 $workDIR/dummy_${force}_step2.nc
                # ice dur loop (1st two elements of timeranges list)
                for i in $(seq $((${hist_sig[0]:0:4}+0)) $((${hist_sig[0]:5:4}-1))); do

                    cdo -O -L timsum -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/dummy_${force}_step2.nc $workDIR/dummy_${endvar}_${force}_step3_$i.nc

                done

        fi

        # remove product of this merge AFTER signal is taken
        cdo mergetime $workDIR/dummy_${endvar}_${force}_step3_*.nc $workDIR/dummy_${endvar}_${force}_step4.nc

        # average ice days per 20 year period
        cdo timmean $workDIR/dummy_${endvar}_${force}_step4.nc $2

        # remove temp files for directory cleanup
        rm $workDIR/dummy_${force}_step1_a.nc
        rm $workDIR/dummy_${force}_step1_b.nc
        rm $workDIR/dummy_${force}_step2.nc
        rm $workDIR/dummy_${endvar}_${force}_step3_*.nc
        rm $workDIR/dummy_${endvar}_${force}_step4.nc
}

# function for calculating ice indexes in the future period
futindex(){            # $1: $workDIR/dummy_${force}_${rcp}_step1.nc
                       # $2: $workDIR/dummy_${endvar}_${force}_${rcp}_step5.nc
        if [ "$endvar" == "icestart" ]
        then

                # select October to December for 1st to last year (select top level - sellevidx,1 - to reduce RAM load assuming it encompasses ice cover dur)
                cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,10/12 -seldate,$((${fut_periods[0]:0:4}+0))-01-01T00:00:00,$((${fut_periods[3]:5:4}+0))-01-01T00:00:00 $1 $workDIR/dummy_${force}_${rcp}_step1_a.nc
                # select January to September for 2nd to last year (lag by 365)
                cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/9 -seldate,$((${fut_periods[0]:0:4}+1))-01-01T00:00:00,$((${fut_periods[3]:5:4}+0))-12-31T00:00:00 $1 $workDIR/dummy_${force}_${rcp}_step1_b.nc
                cdo -O mergetime $workDIR/dummy_${force}_${rcp}_step1_a.nc $workDIR/dummy_${force}_${rcp}_step1_b.nc $workDIR/dummy_${force}_${rcp}_step2.nc
                # ice dur loop (1st two elements of timeranges list)
                for i in $(seq $((${fut_sig[0]:0:4}+0)) $((${fut_sig[0]:5:4}-1))); do

                    cdo -O -L timmin -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/dummy_${force}_${rcp}_step2.nc $workDIR/dummy_${endvar}_${force}_${rcp}_step3_$i.nc

                done



        elif [ "$endvar" == "iceend" ]
        then

                # select October to December for 1st to last year (select top level - sellevidx,1 - to reduce RAM load assuming it encompasses ice cover dur)
                cdo -O -L setctomiss,0 -muldoy -gtc,0 -selmon,9/12 -seldate,$((${fut_periods[0]:0:4}+0))-01-01T00:00:00,$((${fut_periods[3]:5:4}+0))-01-01T00:00:00 $1 $workDIR/dummy_${force}_${rcp}_step1_a.nc
                # select January to September for 2nd to last year (lag by 365)
                cdo -O -L addc,365 -setctomiss,0 -muldoy -gtc,0 -selmon,1/8 -seldate,$((${fut_periods[0]:0:4}+1))-01-01T00:00:00,$((${fut_periods[3]:5:4}+0))-12-31T00:00:00 $1 $workDIR/dummy_${force}_${rcp}_step1_b.nc
                cdo -O mergetime $workDIR/dummy_${force}_${rcp}_step1_a.nc $workDIR/dummy_${force}_${rcp}_step1_b.nc $workDIR/dummy_${force}_${rcp}_step2.nc
                # ice dur loop (1st two elements of timeranges list)
                for i in $(seq $((${fut_sig[0]:0:4}+0)) $((${fut_sig[0]:5:4}-1))); do

                    cdo -O -L timmax -seldate,$i-09-01T00:00:00,$(($i+1))-08-31T00:00:00 $workDIR/dummy_${force}_${rcp}_step2.nc $workDIR/dummy_${endvar}_${force}_${rcp}_step3_$i.nc

                done


        elif [ "$endvar" == "icedur" ]
        then

                # set aside ice days
                cdo -O -L setctomiss,0 -gtc,0 $1 $workDIR/dummy_${force}_${rcp}_step2.nc
                # ice dur loop (1st two elements of timeranges list)
                for i in $(seq $((${fut_sig[0]:0:4}+0)) $((${fut_sig[0]:5:4}-1))); do

                    cdo -O -L timsum -seldate,$i-10-01T00:00:00,$(($i+1))-09-31T00:00:00 $workDIR/dummy_${force}_${rcp}_step2.nc $workDIR/dummy_${endvar}_${force}_${rcp}_step3_$i.nc

                done

        fi

        # merge individual years
        cdo mergetime $workDIR/dummy_${endvar}_${force}_${rcp}_step3_*.nc $workDIR/dummy_${endvar}_${force}_${rcp}_step4.nc

        # average ice days per 30 year period
        cdo timmean $workDIR/dummy_${endvar}_${force}_${rcp}_step4.nc $2

        # remove temp files for directory cleanup
        rm $workDIR/dummy_${force}_${rcp}_step1_a.nc
        rm $workDIR/dummy_${force}_${rcp}_step1_b.nc
        rm $workDIR/dummy_${force}_${rcp}_step2.nc
        rm $workDIR/dummy_${endvar}_${force}_${rcp}_step3_*.nc
        rm $workDIR/dummy_${endvar}_${force}_${rcp}_step4.nc
}

# function for setting attributes on processed files
attribute_setter(){         # $1: $workDIR/${endvar}_${rcp}_step8.nc
                            # $2: $outDIR/${model_fname}_${rcp}_${endvar}_${prod}.nc

    if [ "$endvar" == "icedur" ]
    then

            cdo -O -L setattribute,${endvar}@standard_name="duration_of_lake_ice_cover" -setattribute,${endvar}@long_name="Duration of lake ice cover" -setname,$endvar -setunit,"days" $1 $workDIR/intra_step.nc

    elif [ "$endvar" == "icestart" ]
    then

            cdo -O -L setattribute,${endvar}@standard_name="start_of_lake_ice_cover" -setattribute,${endvar}@long_name="Start of lake ice cover" -setname,$endvar -setunit,"days" $1 $workDIR/intra_step.nc

    elif [ "$endvar" == "iceend" ]
    then

            cdo -O -L setattribute,${endvar}@standard_name="end_of_lake_ice_cover" -setattribute,${endvar}@long_name="End of lake ice cover" -setname,$endvar -setunit,"days" $1 $workDIR/intra_step.nc
    fi

    ncwa -C -v $endvar,lat,lon -a levlak,lev,time,time_bnds $workDIR/intra_step.nc $2
}



# ==============================================================================
# PROCESSING
# ==============================================================================


# go into climate directory (indir up to CLM45)
cd $inDIR
pwd


# ==============================================================================
# HISTORICAL/FUTURE AVERAGES
# ==============================================================================


echo "calculating hist/future ice duration"

for force in "${forcing[@]}"; do
    echo $force

    for scen_folder in "${scenario_folders[@]}"; do

            if [ "$scen_folder" == "historical" ]
            then scen="historical"
                 hist_sig="1971_2000"

                    # merge relevant years (if var is lakeicefrac, first select first vertical level to reduce data load)
                    if [ "$var" == "lakeicefrac" ]
                    then

                            cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${hist_periods[0]}.nc4 $workDIR/dummy_${force}_${hist_periods[0]}_step0.nc
                            cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${hist_periods[1]}.nc4 $workDIR/dummy_${force}_${hist_periods[1]}_step0.nc
                            cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${hist_periods[2]}.nc4 $workDIR/dummy_${force}_${hist_periods[2]}_step0.nc
                            cdo mergetime $workDIR/dummy_${force}_*_step0.nc $workDIR/dummy_${force}_step1.nc

                            rm $workDIR/dummy_${force}_*_step0.nc

                    # merge relevant years
                    elif [ "$var" == "icethick" ]
                    then

                            cdo mergetime $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${hist_periods[0]}.nc4 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${hist_periods[1]}.nc4 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${scen}_${opt}_${var}_global_${tstep}_${hist_periods[2]}.nc4 $workDIR/dummy_${force}_step1.nc

                    fi

                    for endvar in "${endvariables[@]}"; do

                        histindex $workDIR/dummy_${force}_step1.nc $workDIR/dummy_${endvar}_${force}_step5.nc

                    done

                    # remove temp files for directory cleanup
                    rm $workDIR/dummy_${force}_step1.nc


            elif [ "$scen_folder" == "future" ]
            then fut_sig="2070_2099"

                    for rcp in "${rcps[@]}"; do

                        # merge relevant years (if var is lakeicefrac, first select first vertical level to reduce data load)
                        if [ "$var" == "lakeicefrac" ]
                        then

                                cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[0]}.nc4 $workDIR/dummy_${force}_${fut_periods[0]}_step0.nc
                                cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[1]}.nc4 $workDIR/dummy_${force}_${fut_periods[1]}_step0.nc
                                cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[2]}.nc4 $workDIR/dummy_${force}_${fut_periods[2]}_step0.nc
                                cdo sellevidx,1 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[3]}.nc4 $workDIR/dummy_${force}_${fut_periods[2]}_step0.nc
                                cdo mergetime $workDIR/dummy_${force}_*_step0.nc $workDIR/dummy_${force}_${rcp}_step1.nc

                                rm $workDIR/dummy_${force}_*_step0.nc

                        # merge relevant years
                        elif [ "$var" == "icethick" ]
                        then

                                cdo mergetime $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[0]}.nc4 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[1]}.nc4 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[2]}.nc4 $inDIR/$model/$force/$scen_folder/${model_fname}_${force}_ewembi_${rcp}_${opt}_${var}_global_${tstep}_${fut_periods[3]}.nc4 $workDIR/dummy_${force}_${rcp}_step1.nc

                        fi

                        for endvar in "${endvariables[@]}"; do

                            futindex $workDIR/dummy_${force}_${rcp}_step1.nc $workDIR/dummy_${endvar}_${force}_${rcp}_step5.nc

                        done

                        # remove temp files for directory cleanup
                        rm $workDIR/dummy_${force}_${rcp}_step1.nc



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

        for endvar in "${endvariables[@]}"; do

            # subtract for regular ice days signal
            cdo sub $workDIR/dummy_${endvar}_${force}_${rcp}_step5.nc $workDIR/dummy_${endvar}_${force}_step5.nc $workDIR/dummy_${endvar}_${force}_${rcp}_step6.nc

            # divide for scaled icedays signal
            cdo divc,$scaler $workDIR/dummy_${endvar}_${force}_${rcp}_step6.nc $workDIR/dummy_${endvar}_${force}_${rcp}_step7.nc

            rm $workDIR/dummy_${endvar}_${force}_${rcp}_step5.nc

        done
    done
done


# ==============================================================================
# ENSEMBLE MEANS OF SIGNALS
# ==============================================================================


echo "ensemble means"

for rcp in "${rcps[@]}"; do
echo $rcp

    for endvar in "${endvariables[@]}"; do

        # ensemble mean of regular icedays signals (wildcard on force)
        cdo ensmean $workDIR/dummy_${endvar}_*_${rcp}_step6.nc $workDIR/${endvar}_${rcp}_step8.nc
        attribute_setter $workDIR/${endvar}_${rcp}_step8.nc $outDIR/${model_fname}mod_${rcp}_${endvar}_${prod}.nc

        # ensemble mean of scaled icedays signals (wildcard on force)
        cdo ensmean $workDIR/dummy_${endvar}_*_${rcp}_step7.nc $workDIR/${endvar}_${rcp}_step9.nc
        attribute_setter $workDIR/${endvar}_${rcp}_step9.nc $outDIR/${model_fname}mod_${rcp}_${endvar}_scaled_${prod}.nc

        rm $workDIR/dummy_${endvar}_*_${rcp}_step6.nc
        rm $workDIR/${endvar}_${rcp}_step8.nc

        rm $workDIR/dummy_${endvar}_*_${rcp}_step7.nc
        rm $workDIR/${endvar}_${rcp}_step9.nc

    done
done


# ==============================================================================
# CLEANUP
# ==============================================================================


echo "cleanup"

for force in "${forcing[@]}"; do

    for endvar in "${endvariables[@]}"; do

        rm $workDIR/dummy_${endvar}_${force}_step5.nc

    done
done
