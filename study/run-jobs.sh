
# check and wait so as not to exceed slurm job limit
slurm_job_wait() {
    # njobs setter
    if [ -z "$NJOBS_LIMIT" ]
    then
        export NJOBS_LIMIT=100
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

run_takeoff_cases() {
    # All of the standard cases
    for ((i=70; i<=250; i+=10)); do
        export EPZ=$(echo "scale=2; $i/100" | bc)
        export EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            export FARNE=$(echo "scale=2; $j/100" | bc)
            export FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                sbatch -J "combustor-$EPZ-$FARNE" batch.sh
            fi
        done
    done
}


run_lowtol() {
    for ((i=70; i<=250; i+=10)); do
        export EPZ=$(echo "scale=2; $i/100" | bc)
        export EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            export FARNE=$(echo "scale=2; $j/100" | bc)
            export FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir lowtol --npz 1 --precon_off"
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                # echo "$JPATH/thermo-states-$EPZ-$FARNE.yaml"
                sbatch -J "combustor-$EPZ-$FARNE" batch.sh
            fi
        done
    done
}

run_more_cases() {
    # cruise
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir cruise --thrust 0.554 --emodel 5B"
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                sbatch -J "combustor-$EPZ-$FARNE-cruise" batch.sh
            fi
        done
    done

    # Adjust ambient water
    epz=("0.7" "1.5")
    export FARNE="0.10"
    for cepz in "${epz[@]}"; do
        export EPZ=$(printf "%.1f\n" $cepz | sed 's/^0+//')
        for ((j=0; j<=4; j+=1)); do
            export WATER=$(echo "scale=2; $j/100" | bc)
            export WATER=$(printf "%.2f\n" $WATER | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir water --thrust 0.554 --emodel 5B --h2o $WATER --name $WATER"
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                # echo "$JPATH/thermo-states-$EPZ-$FARNE.yaml"
                sbatch -J "combustor-$EPZ-$FARNE-$WATER-cruise" batch.sh
            fi
        done
    done

    # idle
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir idle --thrust 0.06 --emodel 7B --complex_ambient"
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                sbatch -J "combustor-$EPZ-$FARNE-idle" batch.sh
            fi
        done
    done

    #  restarted for start time
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir restarted --restdir cruise --thrust 0.554 --emodel 5B --stime=12.0"
            if ! [ -f "$JPATH/thermo-states-$EPZ-$FARNE.yaml" ]; then
                sbatch -J "combustor-$EPZ-$FARNE-res" batch.sh
            fi
        done
    done

    # run perturbed entrainment cases for idle
    export COMBUSTOR_OPTIONS="--equiv_ratio 1.5 --farnesane 0.1 --outdir entrain --restdir idle --thrust 0.06 --emodel 7B --complex_ambient --entf 0.1"
    sbatch -J "combustor-$EPZ-$FARNE-idle-pet-0.1" batch.sh

    export COMBUSTOR_OPTIONS="--equiv_ratio 1.5 --farnesane 0.1 --outdir entrain --restdir idle --thrust 0.06 --emodel 7B --complex_ambient --entf 10.0"
    sbatch -J "combustor-$EPZ-$FARNE-idle-pet-10" batch.sh
}
