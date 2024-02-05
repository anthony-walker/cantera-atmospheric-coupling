
export EPZ=2.5
export FARNE=0.15
sbatch  -J "combustor-$EPZ-$FARNE" batch.sh

export EPZ=2.5
export FARNE=0.12
sbatch  -J "combustor-$EPZ-$FARNE" batch.sh

export EPZ=2.4
export FARNE=0.13
sbatch  -J "combustor-$EPZ-$FARNE" batch.sh

export EPZ=2.4
export FARNE=0.04
sbatch  -J "combustor-$EPZ-$FARNE" batch.sh

export EPZ=0.7
export FARNE=0.19
sbatch  -J "combustor-$EPZ-$FARNE" batch.sh
