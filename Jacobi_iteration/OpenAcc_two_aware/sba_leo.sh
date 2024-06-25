#!/bin/bash

# Array of node counts
# Array of node counts
export NODE_COUNTS=(1 2 4 8 16)


export MATRIX_SIZE=$1
export ITERATIONS=$2

#OUTPUT DIRECTORY
# Check if directory exists
if [  -d "OUT_$1_$2" ]; then
    rm -r OUT_$1_$2
fi

mkdir OUT_$1_$2


export NTASKS=4
export CPUS_PER_TASK=8
export PARTITION=boost_usr_prod

if [ $3 == "normal" ]; then
    export QOS=normal
else
    export QOS=boost_qos_$3
fi

export ACCOUNT=ict24_dssc_gpu
export GPUS_PER_NODE=4

    


echo -e "\nGPU JACOBI ITERATION:\n"
echo PARTITION: $PARTITION
echo MATRIX_SIZE: $MATRIX_SIZE
echo ITERATIONS: $ITERATIONS
echo -e "CPUS_PER_TASK: $CPUS_PER_TASK"
echo -e "TASKS PER NODE: $NTASKS"
echo -e "GPUS PER NODE: $GPUS_PER_NODE\n"





# Loop through the node counts
for NODES in "${NODE_COUNTS[@]}"
do  
    echo Number of nodes=$NODES
   
    sbatch -A $ACCOUNT -p $PARTITION --qos=$QOS --output=OUT_$1_$2/$NODES.out --job-name=jacobi_acc_$1_$2_$NODES --nodes=$NODES --ntasks-per-node=$NTASKS  -c $CPUS_PER_TASK --gres=gpu:$GPUS_PER_NODE --mem=40G --exclusive ./run_acc.sh 
 
 
    
done