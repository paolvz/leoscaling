#!/bin/bash

# Array of node counts
export NODE_COUNTS=(1 2 4 8 16)

export METHOD=$1
export MATRIX_SIZE=$2

# Set the number of iterations for the average 
if [ $METHOD == "naive" ]; then
    export NUM_AVG=1
else 
    export NUM_AVG=5
fi

#OUTPUT DIRECTORY
# Check if directory exists
if [  -d "OUT_$1_$2" ]; then
    rm -r OUT_$1_$2
fi

mkdir OUT_$1_$2

if [ $1 != "cublas" ]; then
    export NTASKS=1
    export CPUS_PER_TASK=112
    export PARTITION=dcgp_usr_prod
    export ACCOUNT=ict24_dssc_cpu
    if [ $3 == "normal" ]; then
    export QOS=normal
    else
    export QOS=dcgp_qos_$3
    fi
   
else
    export NTASKS=4
    export CPUS_PER_TASK=8
    export PARTITION=boost_usr_prod
    export ACCOUNT=ict24_dssc_gpu
    export GPUS_PER_NODE=4
    if [ $3 == "normal" ]; then
    export QOS=normal
    else
    export QOS=boost_qos_$3
    fi

    
fi

echo -e "\nMATRIX MULTIPLICATION:\n"
echo PARTITION: $PARTITION
echo MATRIX_SIZE: $MATRIX_SIZE
echo METHOD: $METHOD
echo -e "CPUS_PER_TASK: $CPUS_PER_TASK"
echo -e "TASKS PER NODE: $NTASKS\n"





# Loop through the node counts
for NODES in "${NODE_COUNTS[@]}"
do  
    echo Number of nodes=$NODES
    if [ $1 != "cublas" ]; then
        sbatch -A $ACCOUNT -p $PARTITION --qos=$QOS --output=OUT_$1_$2/$NODES.out --job-name=rest_$1_$2_$NODES --nodes=$NODES --ntasks-per-node=$NTASKS  -c $CPUS_PER_TASK --gres=tmpfs:10g --exclusive --mem=481G ./run_$1.sh 
    else
        sbatch -A $ACCOUNT -p $PARTITION --qos=$QOS --output=OUT_$1_$2/$NODES.out --job-name=rest_$1_$2_$NODES --nodes=$NODES --ntasks-per-node=$NTASKS  -c $CPUS_PER_TASK --gres=gpu:$GPUS_PER_NODE --mem=481G --exclusive ./run_$1.sh 
    fi
    
done


