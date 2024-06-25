# Scaling Experiments on Cineca Leonardo Cluster

Code and scripts to run and perform scaling on some classical parallel algorithms (Jacobi Iteration on Heat Equation and Matrix Multiplication) on the Leonardo Supercomputer of Cineca (both on CPU and GPU nodes).

Programming Language: C
Library used: MPI (Cuda - Aware), OpenMP, Cuda and OpenACC

# Tested Algorithms and Folders

- `Jacobi_iteration/`
  - `CPU_one/` - Jacobi Iteration One Side Communication on CPU
  - `CPU_two/` - Jacobi Iteration Two Side Communication on CPU
  - `OpenAcc_two_aware/` - Jacobi Iteration Two Side Communication on GPU (Cuda Aware MPI)
 
- `Matrix_multiplication/` 
  - `*naive*` files - Parallel Mat mul using naive implementation on CPU
  - `*blas*` files - Parallel Mat mul using BLAS routines on CPU
  - `*cublas*` files - Parallel Mat mul using CUBLAS on GPU (non aware MPI)


In each folder (are present the source .c codes, the runners (`run*` files) for loading the correct modules and compiling and the sbatch script (`sba_leo.sh`) to run the scaling.

More detailed info provided in the readme of each implementation.

# Setup 

1) Clone this repo in your leonardo $HOME or $SCRATCH area (assuming you have acces to the cluster) using
```
git clone https://github.com/paolvz/leoscaling.git
```
2) Go in the folder of the algorithm that you prefer and follow the local readme instructions.




