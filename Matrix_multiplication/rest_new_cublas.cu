#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <omp.h>


void check_result(double*C, double*B, int n_row, int n_col, int rank)
{
       int sum = 0;
            for (int i = 0; i < n_row; i++)
            {
                for (int j = 0; j < n_col; j++)
                {
                    if (C[i * n_col + j] == B[i * n_col + j])
                    {   
                        //printf("C[%d][%d] = %g, B[%d][%d] = %g\n", i, j, C[i * n_col + j], i, j, B[i * n_col + j]);
                        sum++;
                    }
                }
            }
            if (sum == n_row * n_col)
        {
            printf("Matrix multiplication is correct, rank: %d\n", rank);
        }
        else
        {
            printf("Matrix multiplication is incorrect, rank: %d\n", rank);
        }
        
}

void print_matrix(double*A, int n_row, int n_col)
{
    for (int i = 0; i < n_row; i++)
    {
        for (int j = 0; j < n_col; j++)
        {
            printf("%g ", A[i * n_col + j]);
        }
        printf("\n");
    }
    printf("\n");
}


    

int main(int argc, char **argv)
{     
    
    int N = MATRIX_SIZE;

    int rank, size;

     
    int thread_level;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

     if (rank == 0){
    printf("GPU CUBLAS Matrix Multiplication\n");
    printf("Matrix Size: %d\n", N);
    }

    double max_comp_time;
    double max_comm_time;
    double max_init_time;
    double max_gpu_cpu_time;
    double final_comp = 0.0;
    double final_comm = 0.0;
    double final_init = 0.0;
    double final_gpu_cpu = 0.0;
    
    char name[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name( name, &len );

    
    

    /////// GPU /////////////
    int gpu_count;
    cudaGetDeviceCount(&gpu_count);

    int gpunumber = rank % gpu_count;
    int device;
    cudaSetDevice(gpunumber);
    cudaGetDevice(&device);
    printf("Rank: %d, GPU: %d, Node: %s\n", rank, device, name);
    ///////////////////////////

    int N_LOC = N / size;
    int rest = N % size;
    int offset = 0;


    if (rank < rest) {N_LOC++;}else{offset = rest;}
 
    // allocate array with malloc
    double*A_LOC = (double*)calloc(N_LOC * N, sizeof(double));
    double*B_LOC = (double*)malloc(N * N_LOC * sizeof(double));
    
    

    // initialize A_LOC and B_LOC
    MPI_Barrier(MPI_COMM_WORLD);
    double start_init = MPI_Wtime();
    #pragma omp parallel for
    for (int i = 0; i < N_LOC; i++)
    {
        for (int j = 0; j < N; j++)
       
    {   if (j == i + N_LOC * rank +offset)
        {

        A_LOC[i*N+j] = 1;

        }
    }
    }

    #pragma omp parallel for
    for (int i = 0; i < N * N_LOC; i++)
    {
        B_LOC[i] = 10 + rank + 2;
    }
    
    double end_init = MPI_Wtime();
    final_init = end_init - start_init;
    

    

    // All Gather

    double* B_TEMP_N = (double*)malloc((long long int)(N / size + 1) * N * sizeof(double));
    double* B_TEMP = (double*)malloc((long long int)(N / size + 1) * N_LOC * sizeof(double));
    double*C_TEMP_N = (double*)malloc((long long int)N_LOC * N * sizeof(double));
    
    
    
    
    ////// GPU ///////
    double *cu_A_LOC;
    double *cu_B_TEMP_N;
    double *cu_C_TEMP_N;
    
    cudaMalloc((void**)&cu_A_LOC, (long long int)N_LOC * N * sizeof(double));
    cudaMalloc((void**)&cu_C_TEMP_N, (long long int)N_LOC * N * sizeof(double));
    cudaMalloc((void**)&cu_B_TEMP_N, (long long int)(N / size + 1) * N * sizeof(double));



    double start_cpu_gpu;
    double end_cpu_gpu;

    start_cpu_gpu = MPI_Wtime();
    
    cudaMemcpy(cu_A_LOC, A_LOC, N_LOC * N * sizeof(double), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    
    end_cpu_gpu = MPI_Wtime();
    final_gpu_cpu += end_cpu_gpu - start_cpu_gpu;
    
    double *cu_C_LOC = cu_C_TEMP_N;
    ////// GPU ///////

    int N_COL = N_LOC;
   
   

    


    for (int k = 0; k < size; k++)
    {
        
        MPI_Barrier(MPI_COMM_WORLD);
        double start_comm = MPI_Wtime();
        
        
        if (k < rest)
        {
            N_COL = N / size + 1;

            cu_C_LOC = cu_C_TEMP_N + (k * N_COL);
            
            

          
        }
        else
        {
            N_COL = N / size;
            
            cu_C_LOC = cu_C_TEMP_N + (k * N_COL + rest);
          
        }
        


        
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N_LOC; i++)
        {
            for (int j = 0; j < N_COL; j++)
            {

                if (k >= rest)
                {

                    B_TEMP[i * N_COL + j] = B_LOC[(k * N_COL + rest) + i * N + j];
                }
                else
                {

                    B_TEMP[i * N_COL + j] = B_LOC[(k * N_COL) + i * N + j];
                }
            }
        }

        int *receivecounts = (int *)malloc(size * sizeof(int));
        int *displs = (int *)malloc(size * sizeof(int));

        int sub_receivecounts = N_COL * N_LOC;
        MPI_Allgather(&sub_receivecounts, 1, MPI_INT, receivecounts, 1, MPI_INT, MPI_COMM_WORLD);

        displs[0] = 0;

        for (int i = 1; i < size; i++)
        {
            displs[i] = displs[i - 1] + receivecounts[i - 1];
        }
        

      
        MPI_Allgatherv(B_TEMP, N_COL * N_LOC, MPI_DOUBLE, B_TEMP_N, receivecounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        
        double end_comm = MPI_Wtime();
        final_comm += end_comm - start_comm;

       
        
        
        /////////// GPU /////////////


        
        
        
        
        
        MPI_Barrier(MPI_COMM_WORLD);

        start_cpu_gpu = MPI_Wtime();

       
        cudaMemcpy(cu_B_TEMP_N, B_TEMP_N, N_COL * N * sizeof(double), cudaMemcpyHostToDevice);
        
        cudaDeviceSynchronize();

        end_cpu_gpu = MPI_Wtime();
        final_gpu_cpu += end_cpu_gpu - start_cpu_gpu;


        
        const double alpha = 1.0;
        const double beta = 0.0;

        cublasHandle_t handle;
        cublasCreate(&handle);
        
    

        double start_comp = MPI_Wtime();

        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N_COL, N_LOC, N, &alpha, cu_B_TEMP_N, N_COL, cu_A_LOC, N, &beta, cu_C_LOC, N);
        cudaDeviceSynchronize();
        
        
        double end_comp = MPI_Wtime();
        final_comp += end_comp - start_comp;


        
        cublasDestroy(handle);
        

  
    }

    /////// GPU /////////////
    MPI_Barrier(MPI_COMM_WORLD);
    start_cpu_gpu = MPI_Wtime();

    cudaMemcpy(C_TEMP_N, cu_C_TEMP_N, N_LOC * N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    end_cpu_gpu = MPI_Wtime();
    final_gpu_cpu += end_cpu_gpu - start_cpu_gpu;
    /////// GPU /////////////
    
    
    //check_result(C_TEMP_N, B_LOC, N, N_LOC, rank);

    #ifdef SAVE
    
        MPI_File file;
        MPI_Offset displacement;
        MPI_Status status;
        
        

        
        displacement = (rank * N_LOC  + offset )* N * sizeof(double);
        
        MPI_File_open(MPI_COMM_WORLD, "result.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
        MPI_File_set_view(file, displacement, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
        MPI_File_write(file, C_TEMP_N, N_LOC * N * sizeof(double), MPI_CHAR, &status);
        MPI_File_close(&file);
    
    #endif
    

    MPI_Reduce(&final_comp, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&final_comm, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&final_init, &max_init_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&final_gpu_cpu, &max_gpu_cpu_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

   
      if (rank == 0)
        {
            printf("\nComputation Time (Max): %f seconds\n", max_comp_time);
            printf("Communication Time (Max): %f seconds\n", max_comm_time);
            printf("Initialization Time (Max): %f seconds\n", max_init_time);
            printf("GPU-CPU Time (Max): %f seconds\n", max_gpu_cpu_time);

            printf("\n%f %f %f %f %s\n", max_init_time, max_comm_time, max_comp_time, max_gpu_cpu_time, getenv("SLURM_NNODES"));
        }
    
    
    free(A_LOC);
    free(B_LOC);
    free(B_TEMP);
    free(B_TEMP_N);
    free(C_TEMP_N);

    cudaFree(cu_A_LOC);
    cudaFree(cu_B_TEMP_N);
    cudaFree(cu_C_TEMP_N);
   
    




 

    MPI_Finalize();





    return 0;
}
