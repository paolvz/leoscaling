#!/bin/bash
echo -e "\n###############################################"
echo job_name: $SLURM_JOB_NAME
echo nodes: $SLURM_NNODES
echo ntasks-per-node: $SLURM_NTASKS_PER_NODE
echo cpus-per-task: $SLURM_CPUS_PER_TASK

export CLUSTER_NAME=$(hostname)

#Check if "leonardo" is in the cluster name
if [[ $CLUSTER_NAME == *"leonardo"* ]]; then
    echo "Running on Leonardo"

    # Load any necessary modules
    module load openmpi/4.1.6--nvhpc--23.11  
    #module load openblas/0.3.24--nvhpc--23.11 
    module load intel-oneapi-mkl/2023.2.0 
fi
echo -e "###############################################\n"


# Set the number of threads to 112
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#COMPILE
export WITHMKL="-L${MKLROOT}/lib/intel64  -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"

#mpirun -np 1 mpicc -O3 -march=native -D MATRIX_SIZE=$MATRIX_SIZE -I${OPENBLASROOT}/include -L/${OPENBLASROOT}/lib -lopenblas -o rest_new_blas_$SLURM_NNODES.x rest_new_blas.c -fopenmp
mpirun -np 1 mpicc -O3 -march=native -D MATRIX_SIZE=$MATRIX_SIZE $WITHMKL -o rest_new_blas_$SLURM_NNODES.x rest_new_blas.c -fopenmp



# RUN
for i in $(seq 1 $NUM_AVG)
do
mpirun -np $SLURM_NNODES --map-by ppr:1:node:pe=$OMP_NUM_THREADS  --bind-to core rest_new_blas_$SLURM_NNODES.x >> OUT_${METHOD}_${MATRIX_SIZE}/avg_${SLURM_NNODES}_${i}.out
done

rm -f rest_new_blas_$SLURM_NNODES.x


if [ $SLURM_NNODES -eq 1 ]; then

sh run_data.sh

fi

