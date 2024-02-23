
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

run_more_cases() {

    run_cruise

    run_perturbed_ambient_water

    run_idle

    run_time_restarted

    run_perturbed_entrainment

}

run_cruise() {
    # cruise
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir cruise --thrust 0.554 --emodel 5B"
            if ! [ -f "cruise/thermo-states-$EPZ-$FARNE.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    echo "Executing on slurm"
                    sbatch -J "combustor-$EPZ-$FARNE-cruise" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    echo "Executing on local"
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-cruise-$EPZ-$FARNE.log" 2>&1 &
                fi
            else
                echo "SKIPPING $COMBUSTOR_OPTIONS"
            fi
        done
    done
}

run_perturbed_ambient_water() {
    # Adjust ambient water
    epz=("0.7" "1.5")
    export FARNE="0.10"
    for cepz in "${epz[@]}"; do
        export EPZ=$(printf "%.1f\n" $cepz | sed 's/^0+//')
        for ((j=0; j<=4; j+=1)); do
            export WATER=$(echo "scale=2; $j/100" | bc)
            export WATER=$(printf "%.2f\n" $WATER | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir water --thrust 0.554 --emodel 5B --h2o $WATER --name $WATER"
            if ! [ -f "water/thermo-states-$EPZ-$FARNE-$WATER.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    echo "Executing on slurm"
                    sbatch -J "combustor-$EPZ-$FARNE-water-$WATER" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    echo "Executing on local"
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-water-$WATER.log" 2>&1 &
                fi
            fi
        done
    done
}

run_idle() {
    # idle
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir idle --thrust 0.06 --emodel 7B --complex_ambient"
            if ! [ -f "idle/thermo-states-$EPZ-$FARNE.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    echo "Executing on slurm"
                    sbatch -J "combustor-$EPZ-$FARNE-idle" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    echo "Executing on local"
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-idle-$EPZ-$FARNE.log" 2>&1 &
                fi
            fi
        done
    done
}

run_time_restarted() {
    #  restarted for start time
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir restarted --restdir cruise --thrust 0.554 --emodel 5B --stime=12.0"
            if ! [ -f "restarted/thermo-states-$EPZ-$FARNE.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    sbatch -J "combustor-$EPZ-$FARNE-res" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-trs-$EPZ-$FARNE.log" 2>&1 &
                fi
            fi
        done
    done
}

run_perturbed_entrainment() {
    # run perturbed entrainment cases for idle
    for j in "0.1" "1.0" "5.0" "10.0"; do
        export COMBUSTOR_OPTIONS="--equiv_ratio 1.5 --farnesane 0.1 --outdir entrain --restdir minimal --thrust 1.0 --emodel 7B --complex_ambient --entf $j --name $j"
        if ! [ -f "entrain/thermo-states-$EPZ-$FARNE.yaml" ]; then
            if [ -f "CLUSTER.txt" ]; then
                echo "Executing on slurm"
                sbatch -J "combustor-$EPZ-$FARNE-idle-pet-0.1" batch.sh
            elif [ -f "TEST.txt" ]; then
                echo "Executing test"
                echo "Running $COMBUSTOR_OPTIONS"
            else
                echo "Executing on local"
                nohup combustor $COMBUSTOR_OPTIONS > "ptb-ent-$j.log" 2>&1 &
            fi
        fi
    done
}

run_takeoff_cases() {
    # All of the standard cases
    for ((i=70; i<=250; i+=10)); do
        export EPZ=$(echo "scale=2; $i/100" | bc)
        export EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            export FARNE=$(echo "scale=2; $j/100" | bc)
            export FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            if ! [ -f "takeoff/thermo-states-$EPZ-$FARNE.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    echo "Executing on slurm"
                    sbatch -J "combustor-$EPZ-$FARNE-takeoff" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    echo "Executing on local"
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-takeoff-$EPZ-$FARNE.log" 2>&1 &
                fi
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
            if ! [ -f "lowtol/thermo-states-$EPZ-$FARNE.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    echo "Executing on slurm"
                    sbatch -J "combustor-$EPZ-$FARNE-lowtol" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    echo "Executing on local"
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-lowtol-$EPZ-$FARNE.log" 2>&1 &
                fi
            fi
        done
    done
}

my_slurm_jobs() {
    squeue --format="%.12i %.35j" -u walkanth # time  %.10M
}

list_slurm_jobs_by_name() {
    if [ -z "$1" ]; then
        echo "No name given, listing all"
    fi
    squeue --format="%.12i %.35j %.10M" -u walkanth | grep -- "$1"
}

cancel_slurm_jobs_by_name() {
    if [ -z "$1" ]; then
        echo "No name given"
        return 1
    fi

    scancel $(my_slurm_jobs | grep -- "$1" | grep -o '\b[0-9]\{7\}\b')
}

run_cruise_nonox() {
    # cruise
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir nonox --thrust 0.554 --emodel 5B --fmodel ../cac/data/combustor.yaml"
            if ! [ -f "nonox/thermo-states-$EPZ-$FARNE.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    echo "Executing on slurm"
                    sbatch -J "combustor-$EPZ-$FARNE-nonox" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    echo "Executing on local"
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-cruise-nonox-$EPZ-$FARNE.log" 2>&1 &
                fi
            else
                echo "SKIPPING $COMBUSTOR_OPTIONS"
            fi
        done
    done
}

run_nox_equilibrium() {
    # cruise
    for ((i=70; i<=250; i+=10)); do
        EPZ=$(echo "scale=2; $i/100" | bc)
        EPZ=$(printf "%.1f\n" $EPZ | sed 's/^0+//')
        for ((j=0; j<=20; j+=1)); do
            FARNE=$(echo "scale=2; $j/100" | bc)
            FARNE=$(printf "%.2f\n" $FARNE | sed 's/^0+//')
            export COMBUSTOR_OPTIONS="--equiv_ratio $EPZ --farnesane $FARNE --outdir noxeq --restdir cruise --thrust 0.554 --emodel 5B --nox -1.0"
            if ! [ -f "noxeq/thermo-states-$EPZ-$FARNE.yaml" ]; then
                if [ -f "CLUSTER.txt" ]; then
                    echo "Executing on slurm"
                    sbatch -J "combustor-$EPZ-$FARNE-noxeq" batch.sh
                elif [ -f "TEST.txt" ]; then
                    echo "Executing test"
                    echo "Running $COMBUSTOR_OPTIONS"
                else
                    echo "Executing on local"
                    nohup combustor $COMBUSTOR_OPTIONS > "ptb-cruise-noxeq-$EPZ-$FARNE.log" 2>&1 &
                fi
            else
                echo "SKIPPING $COMBUSTOR_OPTIONS"
            fi
        done
    done
}
