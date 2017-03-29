#!/bin/bash
#SBATCH -J enkf_3km_NODA_GFS
#SBATCH -N 16 -n 256
#SBATCH -p normal
#SBATCH -t 06:00:00

#load configuration files, functions, parameters

### CM ### 
cd $HOME/WRF_DA/scripts
export CONFIG_FILE=$HOME/WRF_DA/scripts/enkf_NODA_GFS_RI
##########

. $CONFIG_FILE
. util.sh

if [[ ! -d $WORK_DIR ]]; then mkdir -p $WORK_DIR; fi
cd $WORK_DIR

date
#start cycling
export CYCLENUM=0
export DATE=$DATE_START
export PREVDATE=$DATE_START
export NEXTDATE=$DATE_START
while [ $NEXTDATE -le $DATE_END ]; do

#calculate time
#cycle period (CP) for the current cycle
  if [[ $DATE == $DATE_START ]]; then
    export FIRST_CYCLE=true
    if $TWO_WAY_NESTING; then
      export CP=`min ${FIRST_CYCLE_PERIOD[@]}`
    else
      export CP=${FIRST_CYCLE_PERIOD[$RUN_DOMAIN-1]}
    fi
  else
    export FIRST_CYCLE=false
    if $TWO_WAY_NESTING; then
      export CP=`min ${CYCLE_PERIOD[@]}`
    else
      export CP=${CYCLE_PERIOD[$RUN_DOMAIN-1]}
    fi

    # Final WRF/WRF_ENS forecasts
    if [ $CYCLENUM -eq $NUMBER_OF_CYCLES ]; then
      export CP=$((${FORECAST_MINUTES[$RUN_DOMAIN-1]}-${FIRST_CYCLE_PERIOD[$RUN_DOMAIN-1]}-($CYCLENUM-1)*${CYCLE_PERIOD[$RUN_DOMAIN-1]}))
    fi

    # Quit after final WRF/WRF_ENS forecast
    if [ $CYCLENUM -gt $NUMBER_OF_CYCLES ]; then
      break
    fi

  fi

#some other domain-specific parameters
  if $TWO_WAY_NESTING; then
    export OWMAX=`min ${OBS_WIN_MAX[@]}`
    export OWMIN=`min ${OBS_WIN_MIN[@]}`
    export DT=${TIME_STEP[0]}
    #export DT=`min ${TIME_STEP[@]}`
    export MPS=`min ${MINUTE_PER_SLOT[@]}`
    export FCSTM=`min ${FORECAST_MINUTES[@]}`
  else
    export OWMAX=${OBS_WIN_MAX[$RUN_DOMAIN-1]}
    export OWMIN=${OBS_WIN_MIN[$RUN_DOMAIN-1]}
    export DT=${TIME_STEP[0]}
    #export DT=${TIME_STEP[$RUN_DOMAIN-1]}
    export MPS=${MINUTE_PER_SLOT[$RUN_DOMAIN-1]}
    export FCSTM=${FORECAST_MINUTES[$RUN_DOMAIN-1]}
  fi

#calculate start_date and run_minutes, used by namelist_wrf.sh to generate correct
#time in namelist.input
  if $RUN_4DVAR; then
    export start_date_cycle=$DATE
    export run_minutes_cycle=`echo $CP+$OWMAX |bc`
    if ! $FIRST_CYCLE; then
      export start_date_cycle=`advance_time $start_date_cycle $OWMIN`
      export run_minutes_cycle=`echo $run_minutes+$OWMIN |bc`
    fi
  else
    export start_date_cycle=$DATE
    export run_minutes_cycle=$CP
  fi

  if $RUN_DETERMINISTIC; then
    export run_minutes_forecast=`diff_time $DATE $DATE_END`
  else
    export run_minutes_forecast=$FCSTM
  fi

#calculate time when LBC is available, for nested domain runs, the LBC interval is 
#the wrfout_interval from parent domain run.
#if CP < LBINT, the small-period cycle within LBINT will share the same wrfbdy file
#and use wrfinput files from previous wrf run instead if not available directly from real.
  if $TWO_WAY_NESTING || [ $RUN_DOMAIN == 1 ]; then
    export LBINT=$LBC_INTERVAL
  else
    export LBINT=${WRFOUT_INTERVAL[${PARENT_ID[$RUN_DOMAIN-1]}-1]}
  fi
  export minute_off=`echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$LBINT" |bc`
  export LBDATE=`advance_time $start_date_cycle -$minute_off`

#calculate time when parent domain (and outermost domain) output is available
  export PARENTDATE=$DATE_START
  PARENTCP=${CYCLE_PERIOD[${PARENT_ID[$RUN_DOMAIN-1]}-1]}
  while [ `advance_time $PARENTDATE $PARENTCP` -le $DATE ]; do
    export PARENTDATE=`advance_time $PARENTDATE $PARENTCP`
  done
  export OUTERMOSTDATE=$DATE_START
  while [ `advance_time $OUTERMOSTDATE ${CYCLE_PERIOD[0]}` -le $DATE ]; do
    export OUTERMOSTDATE=`advance_time $OUTERMOSTDATE ${CYCLE_PERIOD[0]}`
  done

  export NEXTDATE=`advance_time $DATE $CP`
  echo "----------------------------------------------------------------------"
  echo "CURRENT CYCLE: `wrf_time_string $DATE` => `wrf_time_string $NEXTDATE`"
  echo "CYCLE NUMBER: $CYCLENUM"

  mkdir -p {run,rc,fc,output,obs}/$DATE

#clear error tags
  for d in `ls run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` != "complete" ]]; then
      echo > run/$DATE/$d/stat
    fi
  done

#check if ICBC can be skipped and spin-up file used
  if $SKIP_SPINUP && $FIRST_CYCLE; then
    for d in `ls $WORK_DIR/run/$DATE/`; do
        mkdir -p $WORK_DIR/run/$DATE/$d
        echo "complete" > $WORK_DIR/run/$DATE/$d/stat
    done
    $SCRIPT_DIR/module_skipSpin.sh
    #check dependency
    wait_for_module $WORK_DIR/run/$DATE/skipSpin
  fi

#run components
  if ! $INIT_INNER_DOM; then
    $SCRIPT_DIR/module_wps.sh &
    $SCRIPT_DIR/module_real.sh &
  else
    mkdir -p $WORK_DIR/run/$DATE/real
    echo "complete" > $WORK_DIR/run/$DATE/real/stat
    mkdir -p $WORK_DIR/run/$DATE/wps
    echo "complete" > $WORK_DIR/run/$DATE/wps/stat
  fi

  if $RUN_ENKF && $FIRST_CYCLE; then
    ### CM ### 
    if $RUN_NDOWN_ONLY; then
      $SCRIPT_DIR/module_ndown.sh &
      #check dependency
      wait_for_module $WORK_DIR/run/$DATE/ndown
      mkdir -p $WORK_DIR/run/$DATE/perturb_ic
      echo "complete" > $WORK_DIR/run/$DATE/perturb_ic/stat
    elif  $INIT_INNER_DOM; then
      echo 'Initializing inner domain(s) in WRF'
      mkdir -p $WORK_DIR/run/$DATE/perturb_ic
      echo "complete" > $WORK_DIR/run/$DATE/perturb_ic/stat
    else
      $SCRIPT_DIR/module_perturb_ic.sh &
    fi
    ##########
  fi

  if $RUN_ENKF || $RUN_4DVAR && ! $FIRST_CYCLE; then
    $SCRIPT_DIR/module_obsproc.sh &
  fi

  if $RUN_ENKF && ! $FIRST_CYCLE; then
    $SCRIPT_DIR/module_enkf.sh &
  fi

  if $RUN_4DVAR && ! $FIRST_CYCLE; then
    $SCRIPT_DIR/module_4dvar.sh &
  fi

  if $RUN_ENKF || $RUN_4DVAR && ! $INIT_INNER_DOM; then
    $SCRIPT_DIR/module_update_bc.sh &
  else
    mkdir -p $WORK_DIR/run/$DATE/update_bc
    echo "complete" > $WORK_DIR/run/$DATE/update_bc/stat
  fi

  if $RUN_ENKF && [ $CYCLENUM -lt $NUMBER_OF_CYCLES ]; then
    $SCRIPT_DIR/module_wrf_ens.sh &
  fi

  if $RUN_DETERMINISTIC; then
    $SCRIPT_DIR/module_wrf.sh &
  fi

  wait
  date

#check errors  
  for d in `ls -t run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` == "error" ]]; then
      echo CYCLING STOP DUE TO FAILED COMPONENT: $d
      exit 1
    fi
  done

#advance to next cycle
  export PREVDATE=$DATE
  export DATE=$NEXTDATE

#add count to cycle number
  export CYCLENUM=$(($CYCLENUM+1))

done
echo CYCLING COMPLETE

