#!/bin/bash
echo -e "\n###############################################"
echo job_name: $SLURM_JOB_NAME
echo nodes: $SLURM_NNODES
echo ntasks-per-node: $SLURM_NTASKS_PER_NODE
echo cpus-per-task: $SLURM_CPUS_PER_TASK
echo gpu-per-task: $SLURM_GPUS_PER_TASK

#Check on which machine the job is running
export CLUSTER_NAME=$(hostname)




#Check if "leonardo" is in the cluster name
if [[ $CLUSTER_NAME == *"leonardo"* ]]; then
    echo "Running on Leonardo"
    # Load any necessary modules
    module load openmpi/4.1.6--nvhpc--23.11    
fi

echo -e "###############################################\n"


# Set the number of threads to 64
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


#COMPILE

mpirun -np 1 mpicc -O3 -D MATRIX_SIZE=$MATRIX_SIZE -D ITERATIONS=$ITERATIONS -o jacobi_acc_io_$SLURM_NNODES.x  -Minfo=all -acc  jacobi_acc_io.c


# RUN
for i in 1 2 3 4 5
do
mpirun  --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$OMP_NUM_THREADS --bind-to core jacobi_acc_io_$SLURM_NNODES.x >>  OUT_${MATRIX_SIZE}_${ITERATIONS}/avg_${SLURM_NNODES}_${i}.out
done 

rm -f jacobi_acc_io_$SLURM_NNODES.x


#GENERATE DATA AND PLOT

if [ $SLURM_NNODES -eq 1 ]; then

sh run_data.sh

fi
