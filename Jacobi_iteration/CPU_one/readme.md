# Implementation of Jacobi Iteration using One side MPI comunication on CPU

There are two version of this implementation `jacobi_one_put.c` using the MPI_Put routine and the `jacobi_one_get.c` using the 
MPI_Get routine.

The default scaling settings are:
- 1 2 4 8 16 scaling nodes
- 1 MPI task mapped per node
- 112 OpenMP thread per node (1 per core)

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
NOTE: To plot matrix with size != 62, change the following line with your size in the `plot_io.plt` file before running gnuplot:

```
array=62x62
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

The scaling will be performed on 1 2 4 8 and 16 nodes, you can modify this setting by changing the following line in `sba_leo.sh`:
```
export NODE_COUNTS=(1 2 4 8 16)
```

After the scaling the results will be stored in a folder of the form `OUT_size_iterations_method` (i.e. `OUT_30000_10_put`).

3) Plot the results
In the output folder run:
  ```
  sh gen_data.sh
  gnuplot plot_data.plt
  ```
Will be generated 5 plots:
- plot_comm.png (Communication scaling)
- plot_init.png (Initialization scaling)
- plot_comp.plt (Computation scaling)
- plot.png (Total stacked scaling)
- speedup.png (Speedup plot)

The raw data results are in the `OUT_size_iterations_method.dat` file (i.e. `OUT_30000_10_put.dat`)

























