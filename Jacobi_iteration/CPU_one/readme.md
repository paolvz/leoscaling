# Implementation of Jacobi Iteration using One side MPI comunication on CPU

There are two version of this implementation `jacobi_one_put.c` using the MPI_Put routine and the `jacobi_one_get.c` using the 
MPI_Get routine.

## Launch a single instance of the algorithm in local 

Requirements: OpenMPI and OpenMP installation (latest)

1) Define enviroment variables

- To use MPI_Get
```
export METHOD=put
```
or 
- To use MPI_Put

```
export METHOD=put
```
- Select the MATRIX SIZE (62 suggested to print, matrix with size > 100 will not saved)
```
export MATRIX_SIZE=62
```

- Select the number of Jacobi Iterations (2000 suggested)
```
export ITERATIONS=2000
```

2) Compile using mpicc command
```
mpicc -O3 -D MATRIX_SIZE=$MATRIX_SIZE -D ITERATIONS=$ITERATIONS jacobi_one_$METHOD.c -o jacobi_one_${METHOD}.x -fopenmp

```
3) Run the algorithm in parallel using mpirun with num_proc processes
```
mpirun -np num_proc jacobi_one_${METHOD}.x
```
Results are saved in matrix.bin file.

4) Plot the matrix with gnuplot
```
gnuplot plot_io.plt --persist
```

## Launch scaling on CPU nodes of Leonardo Cluster (DCGP PARTITION)

1) Set your cpu account

Modify the `sba_leo.sh` file by setting 
```
export $ACCOUNT=my_cpu_account
```
with the name of you leonardo CPU account

2) Run the scaling

Running with a MATRIX_SIZE=30000, ITERATIONS=10 and METHOD=put:
```
sh sba_leo.sh 30000 10 put normal
```

























