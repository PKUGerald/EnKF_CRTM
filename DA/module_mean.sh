#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/enkf
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi
if [[ ${DATE:10:2} == '00'  ]]; then exit; fi

#Check dependency
wait_for_module ../../$PREVDATE/enkf 

#Run EnKF
echo running > stat
echo "  Calculating Ensemble Mean..."

domlist=`seq 1 $MAX_DOM`
#link files
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
  if [ -f $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean ]; then continue; fi
  cd $dm
  lfs setstripe -c 1 $rundir/$dm

  echo "    Linking files for domain $dm"
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    ln -fs $WORK_DIR/run/${PREVDATE:0:10}00/wrf_ens/$id/wrfout_${dm}_`wrf_time_string $DATE` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id
    ln -fs $WORK_DIR/run/${PREVDATE:0:10}00/wrf_ens/$id/wrfout_${dm}_`wrf_time_string $DATE` fort.`expr 80010 + $NE`
  done
  wait

  cp -L $WORK_DIR/run/${PREVDATE:0:10}00/wrf_ens/$id/wrfout_${dm}_`wrf_time_string $DATE` fort.`expr 80011 + $NUM_ENS`
  ln -fs $ENKF_DIR/cal_mean.mpi .
  $SCRIPT_DIR/namelist_enkf.sh $n > namelist.enkf
  cd ..
done

#run enkf.mpi
tid=0
nn=$((($enkf_ntasks+$enkf_ppn-$enkf_ntasks%$enkf_ppn)/$enkf_ppn))
nt=$(($total_ntasks/$HOSTPPN/$nn))
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  echo "    Running cal_mean.mpi for domain $dm"
  $SCRIPT_DIR/job_submit.sh $enkf_ntasks $((tid*$enkf_ntasks)) $enkf_ppn ./cal_mean.mpi >& cal_mean.log &
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

for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  cp $dm/fort.`expr 80011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean
  ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean
  ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean $WORK_DIR/fc/$DATE/wrfinput_${dm}
done


echo complete > stat

