# Implementation of Jacobi Iteration using Two side MPI (cuda-aware) communication and OpenACC on GPU

Classical Parallel Jacobi Iteration using (cuda-aware) MPI communications routine and OpenACC.

The default scaling settings are:
- 1 2 4 8 16 scaling nodes
- 4 MPI task mapped per node (1 per GPU)
- 8 OpenMP thread per task (32 per node, 1 per core)


## Launch scaling on GPU nodes of Leonardo Cluster (BOOSTER PARTITION)

1) Set your cpu account

Modify the `sba_leo.sh` file by setting 
```
export $ACCOUNT=my_cpu_account
```
with the name of you leonardo CPU account

2) Run the scaling

Running with a MATRIX_SIZE=30000, ITERATIONS=10:
```
sh sba_leo.sh 30000 10 normal
```

The scaling will be performed on 1 2 4 8 and 16 nodes, you can modify this setting by changing the following line in `sba_leo.sh`:
```
export NODE_COUNTS=(1 2 4 8 16)
```

After the scaling the results will be stored in a folder of the form `OUT_size_iterations` (i.e. `OUT_30000_10`).

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

The raw data results are in the `OUT_size_iterations.dat` file (i.e. `OUT_30000_10.dat`)

NOTE: *For each node we run 5 times and take the average performances*
