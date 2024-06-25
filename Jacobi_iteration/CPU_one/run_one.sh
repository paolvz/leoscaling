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
else
    echo "Running on Orfeo"
    module load cuda/12.1
    module load openMPI/4.1.5/gnu/12.2.1
    
fi

echo -e "###############################################\n"



export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


#COMPILE

mpirun -np 1 mpicc -O3 -D MATRIX_SIZE=$MATRIX_SIZE -D ITERATIONS=$ITERATIONS jacobi_one_$METHOD.c -o jacobi_one_${METHOD}_${SLURM_NNODES}.x -fopenmp


# RUN
for i in 1 2 3 4 5
do
mpirun -np $SLURM_NNODES --map-by ppr:1:node:pe=$OMP_NUM_THREADS --bind-to core jacobi_one_${METHOD}_${SLURM_NNODES}.x  >> OUT_${MATRIX_SIZE}_${ITERATIONS}_${METHOD}/avg_${SLURM_NNODES}_${i}.out
done

rm -f jacobi_one_${METHOD}_${SLURM_NNODES}.x


#GENERATE DATA AND PLOT




if [ $SLURM_NNODES -eq 1 ]; then

sh run_data.sh

fi


   









