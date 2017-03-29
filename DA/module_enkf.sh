#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/enkf
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../../$PREVDATE/wrf_ens ../obsproc 


#Run EnKF
echo running > stat
echo "  Running EnKF..."

domlist=`seq 1 $MAX_DOM`
#link files
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  lfs setstripe -c 1 $rundir/$dm

  echo "    Linking files for domain $dm"
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id
    ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE`
    cp -L fort.`expr 80010 + $NE` fort.`expr 90010 + $NE` >> link.log 2>&1 &
  done
  wait

# linking mean by Minamide 2015.2.25
#  cp -L $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean
  cp -L $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.`expr 80011 + $NUM_ENS`
  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 90011 + $NUM_ENS` 
  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 60011 + $NUM_ENS`

#  cp -L fort.80011 fort.`expr 80011 + $NUM_ENS`
#  cp -L fort.80011 fort.`expr 90011 + $NUM_ENS`
  ln -fs $WRF_DIR/run/LANDUSE.TBL .
  ln -fs $CRTM_DIR/crtm_wrf/coefficients .
  ln -sf $CRTM_DIR/crtm_wrf/gt_lut .
  #Observations
  #LITTLE_R format from obsproc
  ln -fs $DATA_DIR/obs/${DATE:0:4}/obs_gts_`wrf_time_string $DATE`.3DVAR obs_3dvar_${DATE}00
  #airborne radar superobs
  #ln -fs $DATA_DIR/so/2010/${DATE}_all.so_ass airborne_${DATE}_so
  ln -fs $DATA_DIR/so/${DATE:0:4}/${DATE}_all.so_ass airborne_${DATE}_so
  ln -fs $DATA_DIR/GOES/radiance_${dm}_${DATE}_so radiance_${DATE}_so



  #different strategy with minSLP
  if [[ ${DATE:10:2} == '00'  ]]; then
    ln -fs $ENKF_DIR/enkf.mpi .
  else
    #ln -fs $ENKF_DIR/enkf_minslp.mpi enkf.mpi
    ln -fs $ENKF_DIR/enkf_hydro.mpi enkf.mpi
  fi

  $SCRIPT_DIR/namelist_enkf.sh $n > namelist.enkf
  cd ..
done

# replacing prior mean with determinstic forecast
if $REPLACE_MEAN; then
 if [[ $REPLACE_MEAN_WITH == "prior_forecast" ]]; then
 #if [[ ${DATE:10:2} == '00'  ]]; then
  tid=0
  nn=$((($enkf_ntasks+$HOSTPPN-$enkf_ntasks%$HOSTPPN)/$HOSTPPN))
  nt=$(($total_ntasks/$HOSTPPN/$nn))
#  for n in $domlist; do
   n=3
   dm=d`expr $n + 100 |cut -c2-`
   cd $dm
   if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
    cd replace_mean
    echo "  Replacing ens mean with $REPLACE_MEAN_WITH for domain $dm"
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv ../fort.`expr 80010 + $NE` fort.`expr 80010 + $NE`
      cp -L fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
    done
    if [[ $REPLACE_MEAN_WITH == "prior_forecast" ]]; then
      #ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.70010
      ln -fs $WORK_DIR/rc/$DATE/wrfinput_$dm fort.70010
    fi
    ln -fs $ENKF_DIR/replace_mean.exe .
    ln -sf $CRTM_DIR/crtm_wrf/gt_lut .
    export SLURM_TASKS_PER_NODE="$HOSTPPN(x$SLURM_NNODES)"
    ibrun -n $enkf_ntasks -o $((tid*$enkf_ntasks)) ./replace_mean.exe $NUM_ENS >& replace_mean.log &
    export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
    tid=$((tid+1))
    if [[ $tid == $nt ]]; then
      tid=0
      wait
    fi
    #./replace_mean.exe $NUM_ENS >& replace_mean.log
    #watch_log replace_mean.log Successful 1 $rundir
    #for NE in `seq 1 $((NUM_ENS+1))`; do
    #  mv fort.`expr 90010 + $NE` ../fort.`expr 80010 + $NE`
    #  cp ../fort.`expr 80010 + $NE` ../fort.`expr 90010 + $NE` 
    #done
    cd ../..
#  done

#  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm/replace_mean/
    watch_log replace_mean.log Successful 1 $rundir
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv fort.`expr 90010 + $NE` ../fort.`expr 80010 + $NE`
      cp ../fort.`expr 80010 + $NE` ../fort.`expr 90010 + $NE`
    done
    cd ../..
#  done
 #fi
 fi
fi

#run enkf.mpi
tid=0
nn=$((($enkf_ntasks+$enkf_ppn-$enkf_ntasks%$enkf_ppn)/$enkf_ppn))
nt=$(($total_ntasks/$HOSTPPN/$nn))
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  echo "    Running enkf.mpi for domain $dm"
  $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./enkf.mpi >& enkf.log &
  tid=$((tid+1))
  if [[ $tid == $nt ]]; then
    tid=0
    wait
  fi
  cd ..
done
wait

#Check output
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
#  watch_log $dm/enkf.log Successful 5 $rundir
  watch_log $dm/${DATE}.finish_flag _ 5 $rundir
done

# replacing mean with first guess (GFS/FNL) reanalysis
if $REPLACE_MEAN; then
 if [[ $REPLACE_MEAN_WITH != "prior_forecast" ]]; then
  tid=0
  nn=$((($enkf_ntasks+$HOSTPPN-$enkf_ntasks%$HOSTPPN)/$HOSTPPN))
  nt=$(($total_ntasks/$HOSTPPN/$nn))
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    if [[ ! -d replace_mean ]]; then mkdir -p replace_mean; fi
    cd replace_mean
    echo "  Replacing ens mean with $REPLACE_MEAN_WITH for domain $dm"
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv ../fort.`expr 90010 + $NE` fort.`expr 80010 + $NE`
      cp fort.`expr 80010 + $NE` fort.`expr 90010 + $NE`
    done
    if [[ $REPLACE_MEAN_WITH == "forecast" ]]; then
      ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.70010
    elif [[ $REPLACE_MEAN_WITH == "gfs" ]]; then
      ln -fs $WORK_DIR/rc/$DATE/wrfinput_$dm fort.70010
    fi
    ln -fs $ENKF_DIR/replace_mean.exe .
    export SLURM_TASKS_PER_NODE="$HOSTPPN(x$SLURM_NNODES)"
    ibrun -n $enkf_ntasks -o $((tid*$enkf_ntasks)) ./replace_mean.exe $NUM_ENS >& replace_mean.log &
    export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS/$SLURM_NNODES))(x$SLURM_NNODES)"
    tid=$((tid+1))
    if [[ $tid == $nt ]]; then
      tid=0
      wait
    fi
    cd ../..
  done
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm/replace_mean/
    #./replace_mean.exe $NUM_ENS >& replace_mean.log
    watch_log replace_mean.log Successful 1 $rundir
    for NE in `seq 1 $((NUM_ENS+1))`; do
      mv fort.`expr 90010 + $NE` ../
    done
    cd ..
    cd ..
  done
 fi
fi

for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -sf $rundir/$dm/fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id
#    cp $dm/fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id
###    ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id $dm/fort.`expr 90010 + $NE`
    ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
    if [ ! -f $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id ]; then
      ln -fs $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id
    fi
  done
  cp $dm/fort.`expr 80011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean
  cp $dm/fort.`expr 90011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean
  ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean $WORK_DIR/fc/$DATE/wrfinput_${dm}
done


echo complete > stat

