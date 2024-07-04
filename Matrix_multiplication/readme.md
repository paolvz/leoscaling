# Implementation of Parallel Matrix Multiplication (CPU and GPU)

There are three version of this implementation:
- `rest_new_naive.c` naive parallel matrix multiplication on CPU
- `rest_new_naive.c` blas parallel matrix multiplication on CPU (using MKL or OpenBLAS)
- `rest_new_naive.c` CUDA cublas parallel matrix multiplication on GPU
- 
The default scaling settings are:
- CPU  
  - 1 2 4 8 16 scaling nodes
  - 1 MPI task mapped per node
  - 112 OpenMP thread per node (1 per core)

- GPU 
  - 1 2 4 8 16 scaling nodes
  - 4 MPI task mapped per node
  - 8 OpenMP threads per task (32 per node, 1 per core)
 
    
## Launch a single instance of the algorithm in local on CPU

Requirements: OpenMPI and OpenMP installation (latest)
NOTE: *These are instructions only for the naive version. 
The steps for the blas version are the same but you have to compile the `rest_new_blas.c` file by linking the Intel MKL or OpenBLAS library.*

1) Define enviroment variables

Seleect the MATRIX SIZE 
```
export MATRIX_SIZE=100
```



2) Compile using mpicc command (use -DSAVE if you want to save the result)
```
mpicc rest_new_naive.c -o rest_new_naive.x -fopenmp -DMATRIX_SIZE=$MATRIX_SIZE -DSAVE

```
3) Run the algorithm in parallel using mpirun with num_proc processes
```
mpirun -np num_proc rest_new_naive.c
```
Results are saved in `result.bin` file.

4) Read the result (choose your rows and columns)
```
gcc read_matrix.c -o read_matrix.x -DROWS=100 -DCOLS=100
./read_matrix.x
```


## Launch scaling on nodes of Leonardo Cluster (DCGP AND BOOSTER PARTITION)

1) Set your cpu account

Modify the `sba_leo.sh` file by setting 
```
export $ACCOUNT=my_cpu_account
```
with the name of you leonardo CPU account

2) Run the scaling
NOTE: *By default the blas version uses MKL, you can change to OpenBLAS by de-comment the rows related in `sbaleo.sh` and comment the MKL related rows*

Running with METHOD=blas and MATRIX_SIZE=30000:
```
sh sba_leo.sh 30000 blas normal
```

Other options for method: `naive` and `cublas`

The scaling will be performed on 1 2 4 8 and 16 nodes, you can modify this setting by changing the following line in `sba_leo.sh`:
```
export NODE_COUNTS=(1 2 4 8 16)
```

After the scaling the results will be stored in a folder of the form `OUT_method_size` (i.e. `OUT_blas_30000`).

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
- plot_cpu_gpu.png (Communication From/To GPU, only for cublas)

The raw data results are in the `OUT_size_method.dat` file (i.e. `OUT_blas_30000.dat`)

NOTE: *For blas and cublas implementation for each node we run 5 times and take the average times.*



