#!/bin/bash

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


export NTASKS=1
export CPUS_PER_TASK=112
export PARTITION=dcgp_usr_prod
if [ $3 == "normal" ]; then
    export QOS=normal
else
    export QOS=dcgp_qos_$3
fi

export ACCOUNT=ict24_dssc_cpu


    


echo -e "\nCPU JACOBI ITERATION:\n"
echo PARTITION: $PARTITION
echo QOS: $QOS
echo MATRIX_SIZE: $MATRIX_SIZE
echo ITERATIONS: $ITERATIONS    
echo -e "CPUS_PER_TASK: $CPUS_PER_TASK"
echo -e "TASKS PER NODE: $NTASKS\n"






# Loop through the node counts
for NODES in "${NODE_COUNTS[@]}"
do  
    echo Number of nodes=$NODES
   
    sbatch -A $ACCOUNT -p $PARTITION --qos=$QOS --output=OUT_$1_$2/$NODES.out --job-name=jacobi_cont_$1_$2_$NODES --nodes=$NODES --ntasks-per-node=$NTASKS  --gres=tmpfs:10g -c $CPUS_PER_TASK --mem=160G --exclusive ./run_cont.sh 
 
 
    
done



