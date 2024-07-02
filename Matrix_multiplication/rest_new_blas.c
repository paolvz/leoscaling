#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <mkl.h>





void check_result(double*C, double*B, int n_row, int n_col, int rank)
{
       int sum = 0;
            for (long long int i = 0; i < n_row; i++)
            {
                for (long long int j = 0; j < n_col; j++)
                {
                    if (C[i * n_col + j] == B[i * n_col + j])
                    {
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
    for (long long int i = 0; i < n_row; i++)
    {
        for (long long int j = 0; j < n_col; j++)
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
    printf("Naive Matrix Multiplication\n");
    printf("Matrix Size: %d\n", N);
    }

    double max_comp_time;
    double max_comm_time;
    double max_init_time;
    double final_comp = 0.0;
    double final_comm = 0.0;
    double final_init = 0.0;

    if (rank == 0){
    #pragma omp parallel
    {
    #pragma omp single
    {
    printf("OMP Threads: %d\n", omp_get_num_threads());
    }}}

    int N_LOC = N / size;
    int rest = N % size;
    int offset = 0;


    if (rank < rest) {N_LOC++;}else{offset = rest;}
 
    // allocate array with malloc
    double*A_LOC = (double*)calloc((long long int)N_LOC * N, sizeof(double));
    double*B_LOC = (double*)malloc((long long int)N * N_LOC * sizeof(double));
    
    

    // initialize A_LOC and B_LOC
    MPI_Barrier(MPI_COMM_WORLD);
    double start_init = MPI_Wtime();

    #pragma omp parallel for collapse(2)
    for (long long int i = 0; i < N_LOC; i++)
    {
        for (long long int j = 0; j < N; j++)
       
    {   if (j == i + N_LOC * rank +offset)
        {

        A_LOC[i*N+j] = 1;

        }
    }
    }
    
    #pragma omp parallel for
    for (long long int i = 0; i < N * N_LOC; i++)
    {
        B_LOC[i] = 10 + rank + 2;
    }

    double end_init = MPI_Wtime();
    final_init = end_init - start_init;
    

    

    // All Gather

    double* B_TEMP_N = (double*)malloc((long long int)(N / size + 1) * N * sizeof(double));
    double* B_TEMP = (double*)malloc((long long int)(N / size + 1) * N_LOC * sizeof(double));
    double*C_TEMP_N = (double*)malloc((long long int)N_LOC * N * sizeof(double));

    
    double* C_LOC;
    int N_COL = N_LOC;
    
    if (size > 1){
    for (int k = 0; k < size; k++)
    {
        
        MPI_Barrier(MPI_COMM_WORLD);    
        double start_comm = MPI_Wtime();
  
        if (k < rest)
        {
            N_COL = N / size + 1;
            C_LOC = C_TEMP_N + (N_COL * k);
        }
        else
        {
            N_COL = N / size;
            C_LOC = C_TEMP_N + (N_COL * k + rest);
        }
        


        
        #pragma omp parallel for collapse(2)
        for (long long int i = 0; i < N_LOC; i++)
        {
            for (long long int j = 0; j < N_COL; j++)
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

  

        MPI_Barrier(MPI_COMM_WORLD);
        double start_comp = MPI_Wtime();
        
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_LOC, N_COL, N, 1, A_LOC, N, B_TEMP_N, N_COL, 0, C_LOC, N);
        double end_comp = MPI_Wtime();
        final_comp += end_comp - start_comp; 
        
        
    }

    free(B_TEMP);
    free(B_TEMP_N);

    } else {
        double start_comp = MPI_Wtime();
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_LOC, N, N, 1, A_LOC, N, B_LOC, N, 0, C_TEMP_N, N);
        double end_comp = MPI_Wtime();
        final_comp += end_comp - start_comp; 
    }

    //check_result(C_TEMP_N, B_LOC, N, N_LOC, rank);




    MPI_Reduce(&final_comp, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&final_comm, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&final_init, &max_init_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    

    if (rank == 0)
        {
            printf("\nComputation Time (Max): %f seconds\n", max_comp_time);
            printf("Communication Time (Max): %f seconds\n", max_comm_time);
            printf("Initialization Time (Max): %f seconds\n", max_init_time);

            printf("\n%f %f %f %s\n", max_init_time, max_comm_time, max_comp_time, getenv("SLURM_NNODES"));
        }
    

    free(A_LOC);
    free(B_LOC);
    free(C_TEMP_N);
    

    MPI_Finalize();





    return 0;
}
