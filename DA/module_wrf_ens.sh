#!/bin/bash
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/wrf_ens
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../update_bc 

#Setup for wrf run
echo running > stat
echo "  Running WRF ensemble..."

tid=1  #does not start from 0, because the wrf forecast runs with ens at the same time.
nt=$((total_ntasks/$wrf_ntasks))
for re_run in `seq 1 3`; do
  tid=1
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    if [[ ! -d $id ]]; then mkdir $id; fi
    touch $id/rsl.error.0000
    if [[ `tail -n1 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

    cd $id
    lfs setstripe -c 1 $rundir/$id
    ln -fs $WRF_DIR/run/* .
    rm -f namelist.*

    for n in `seq 1 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id wrfinput_$dm
    done
    ln -fs $WORK_DIR/fc/$DATE/wrfbdy_d01_$id wrfbdy_d01

    if [[ $SST_UPDATE == 1 ]]; then
      ln -fs $WORK_DIR/rc/$LBDATE/wrflowinp_d?? .
    fi

    export start_date=$start_date_cycle
    export run_minutes=$run_minutes_cycle 
    export inputout_interval=$run_minutes
    export inputout_begin=0
    export inputout_end=$run_minutes
    export GET_PHYS_FROM_FILE=false
    $SCRIPT_DIR/namelist_wrf.sh wrfw $RUN_DOMAIN > namelist.input

#    if [[ $re_run == 1 ]]; then
    $SCRIPT_DIR/job_submit.sh $wrf_ntasks $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
    tid=$((tid+1))
    if [[ $tid == $nt ]]; then
      tid=1
      wait
    fi
#    else
#
### for time-saving
#      export t_rerun=`echo "2^$((re_run-1))"|bc`
#      export tid_rerun=`min $t_rerun $nt`
#      $SCRIPT_DIR/job_submit.sh $((tid_rerun*$wrf_ntasks)) $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
#      tid=$((tid+$tid_rerun))
#      if [[ $tid >= $nt ]]; then
#        tid=$tid_rerun
#        wait
#      fi
### end time-saveing
#    fi
    cd ..
  done
wait
## next line is for time-saving
#export t_rerun=`echo "2^$re_run"|bc`
#export tid=`min $t_rerun $nt`
done

for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  watch_log $id/rsl.error.0000 SUCCESS 1 $rundir
  outfile=$id/wrfinput_d01_`wrf_time_string $NEXTDATE`
#  watch_file $outfile 1 $rundir
  mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $NEXTDATE`_$id
  if [ $MAX_DOM -gt 1 ]; then
    for n in `seq 2 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      outfile=$id/wrfout_${dm}_`wrf_time_string $NEXTDATE`
      watch_file $outfile 1 $rundir
###      mv $outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $NEXTDATE`_$id
      ln -fs $rundir/$outfile $WORK_DIR/fc/$DATE/wrfinput_${dm}_`wrf_time_string $NEXTDATE`_$id
    done
  fi
done

if $CLEAN; then rm $rundir/$id/wrfout* ; fi
echo complete > stat
