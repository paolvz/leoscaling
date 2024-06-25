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
    # Set the CUDA+MPI path
    export CUDAMPI="-I/leonardo/prod/opt/libraries/openmpi/4.1.6/nvhpc--23.11/include -L/leonardo/prod/opt/libraries/openmpi/4.1.6/nvhpc--23.11/lib -lmpi -lcublas"
    export ARCH=sm_80
fi

echo -e "###############################################\n"


# Set the number of threads to 64
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


#COMPILE


mpirun -np 1 nvcc -O3 -D MATRIX_SIZE=$MATRIX_SIZE -arch $ARCH $CUDAMPI rest_new_cublas.cu -o rest_new_cublas_$SLURM_NNODES.x -Xcompiler -fopenmp


# RUN
for i in $(seq 1 $NUM_AVG)
do
mpirun  --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$OMP_NUM_THREADS --bind-to core rest_new_cublas_$SLURM_NNODES.x >> OUT_${METHOD}_${MATRIX_SIZE}/avg_${SLURM_NNODES}_${i}.out
done


rm -f rest_new_cublas_$SLURM_NNODES.x


#GENERATE DATA AND PLOT

if [ $SLURM_NNODES -eq 1 ]; then

sh run_data.sh

fi
