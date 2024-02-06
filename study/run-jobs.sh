
# check and wait so as not to exceed slurm job limit
slurm_job_wait() {
    # njobs setter
    if [ -z "$NJOBS_LIMIT" ]
    then
        export NJOBS_LIMIT=1000
    fi
    # run the job with set options
    if [ -z "$SKIP_SBATCH" ]
    then
        n_jobs=$(squeue -u walkanth | wc -l)
        stime=1
        while [ $n_jobs -ge $NJOBS_LIMIT ]
        do
            echo "Number of jobs over $NJOBS_LIMIT, waiting $stime seconds..."
            sleep $stime
            n_jobs=$(squeue -u walkanth | wc -l)
            # adjust sleep time
            if [ $stime -lt 10 ]
            then
                ((stime=stime+1))
            fi
        done
        echo "Currently queued jobs: $n_jobs"
    else
        sleep $SLEEP_TIMER
    fi
}

count_jobs() {
    squeue -h -o "%A" -u $(whoami) | wc -l
}

run_all_cases() {
    for ((i=70; i<=250; i+=10)); do
        export EPZ=$(echo "scale=2; $i/100" | bc)
        export EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            export FARNE=$(echo "scale=2; $j/100" | bc)
            export FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                # echo "$JPATH/thermo-states-$EPZ-$FARNE.yaml"
                sbatch -J "combustor-$EPZ-$FARNE" batch.sh
            fi
        done
    done
}

run_restarts() {
    for ((i=70; i<=250; i+=10)); do
        export EPZ=$(echo "scale=2; $i/100" | bc)
        export EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            export FARNE=$(echo "scale=2; $j/100" | bc)
            export FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                # echo "$JPATH/thermo-states-$EPZ-$FARNE.yaml"
                sbatch -J "combustor-$EPZ-$FARNE" restarted.sh
            fi
        done
    done
}
