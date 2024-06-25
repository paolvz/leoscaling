#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>




int main(int argc, char **argv)
{
    
       
    int thread_level;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);
    
    int rank, size;

    int N = MATRIX_SIZE;
    int NUM_ITER = ITERATIONS;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N_LOC = N / size ;
    int rest = N % size;
    int offset = 0;
    
    double max_init_time;
    double max_comm_time;
    double max_comp_time;

    double final_comm = 0;
    double final_comp = 0;
    
    if (rank == 0)
    {
        printf("\nMatrix size: %d\n", N);
        printf("Iterations: %d\n", NUM_ITER);
    }
   


   if (rank == 0){
    #pragma omp parallel
    {
    #pragma omp single
    {
    printf("OMP Threads: %d\n", omp_get_num_threads());
    }}}

   

    if (rank != 0 && rank != size - 1) {N_LOC += 2;}
    else {N_LOC += 1;}


    if (rank < rest) {N_LOC++;} else {offset = rest;}
    
    double *MAT_LOC = (double *)malloc((long long int)N_LOC * N * sizeof(double));
    double *MAT_LOC_NEW = (double *)malloc((long long int) N_LOC * N * sizeof(double));
    
    MPI_Barrier(MPI_COMM_WORLD);
    double start_init = MPI_Wtime();

    //// INNER INITIALIZATION ////
    #pragma omp parallel for
    for (long long int i = 0; i < N_LOC; i++)
    {
        for (long long int j = 0; j < N; j++)
        {
            MAT_LOC[i * N + j] = 0.5;
            MAT_LOC_NEW[i * N + j] = 0.5;
        }
    }

    // TOP ROW = 0
    if (rank == 0){
    
    #pragma omp parallel for
    for (long long int i = 0; i < N; i++)
        {
            MAT_LOC[i] = 0;
            MAT_LOC_NEW[i] = 0;
        }
    
    }

    
    

  

    //LINEARY SPACED BOTTOM COLUMN
    if (rank == size-1){
    #pragma omp parallel for
    for (long long int i = 0; i < N; i++)
    {   double step = (100.0 - 0) / (N - 1);
        MAT_LOC[(long long int)(N_LOC - 1) * N + i] = 100 - i * step;
        MAT_LOC_NEW[(long long int)(N_LOC - 1) * N + i] = 100 - i * step;
    }
    }

    
    //LINEARY SPACED LEFT COLUMN (BY RANKS)
    #pragma omp parallel for
    for (long long int i = 1; i < N_LOC-1; i++)
    {   
        long long int global_idx = (N_LOC-2)*rank+(i-1)+offset;
        
        if (rank == 0){
            global_idx = (N_LOC)*rank+(i)+offset;
        } 
        else if (rank == size-1){
            global_idx = (N_LOC-2)*rank+(i+1)+offset;
        }

        double step = (100.0 - 0) / (N  - 1);
        MAT_LOC[(long long int)i*N] = global_idx  * step ;
        MAT_LOC_NEW[(long long int)i*N] = global_idx  * step ;
    }


    //RIGHT COLUMN = 0 (BY RANKS)
    #pragma omp parallel for
    for (long long int i = 0; i < N_LOC; i++)
    {
        MAT_LOC[(long long int)i*N + N - 1] = 0;
        MAT_LOC_NEW[(long long int)i*N + N - 1] = 0;
    }

    double end_init = MPI_Wtime();
    double final_init = end_init - start_init;

    /////////////////////////
   

 



for (int iter= 0; iter <= NUM_ITER; iter++)

{   
     //// COMMUNICATION ////
    
    

    
    
    double * ghost_up = MAT_LOC;
    double * ghost_down = MAT_LOC + (N_LOC-1) * N;

    
    double* first_row_point = MAT_LOC + N;
    double* last_row_point = MAT_LOC + (N_LOC-2) * N;
   
    
  

   
    

    
   
    MPI_Barrier(MPI_COMM_WORLD);
    double start_comm = MPI_Wtime();

    MPI_Request request_up, request_down;
    if (rank != 0) {
        MPI_Isend(first_row_point, N, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &request_up);
        MPI_Recv(ghost_up, N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if (rank != size - 1) {
        MPI_Isend(last_row_point, N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request_down);
        MPI_Recv(ghost_down, N, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if (rank != 0) {
        MPI_Wait(&request_up, MPI_STATUS_IGNORE);
    }
    
    if (rank != size - 1) {
        MPI_Wait(&request_down, MPI_STATUS_IGNORE);
    }

    double end_comm = MPI_Wtime();
    final_comm += end_comm - start_comm;
    
    ///////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);
    double start_comp = MPI_Wtime();

    ///// COMPUTATION /////
    #pragma omp parallel for collapse(2)
    for (long long int i = 1; i < N_LOC-1; i++)
    {
        for (long long int j = 1; j < N-1; j++)
        {   
            MAT_LOC_NEW[(long long int)i * N + j] = 0.25 * (MAT_LOC[(long long int)(i - 1) * N + j] + MAT_LOC[(long long int)(i + 1) * N + j] + MAT_LOC[(long long int)i * N + j - 1] + MAT_LOC[(long long int)i * N + j + 1]);
        }
    }
    
    double end_comp = MPI_Wtime();
    final_comp += end_comp - start_comp;

    ///////////////////////////
    
    
    // SWAP ARRAYS USING POINTERS
    double *temp = MAT_LOC;
    MAT_LOC = MAT_LOC_NEW;
    MAT_LOC_NEW = temp;



}
    
  

   

   
   
    
   if (N < 100){

 
   MPI_File file;
   MPI_Status status;
   
   int FINAL_N_LOC = N_LOC;
    if (rank != 0 & rank != size - 1) {FINAL_N_LOC -=2;}
    else {FINAL_N_LOC-=1;}

   int displacement = (rank * FINAL_N_LOC + offset) * N * sizeof(double);

    double* final_pointer;

    if (rank == 0){
        final_pointer = MAT_LOC;
    }
    else{
        final_pointer = MAT_LOC + N ;
    }
  
    
    
    
    MPI_File_open(MPI_COMM_WORLD, "matrix.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_File_set_view(file, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_write(file, final_pointer, FINAL_N_LOC * N, MPI_DOUBLE, &status);
    MPI_File_close(&file);
    
   }

    MPI_Reduce(&final_comp, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&final_comm, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&final_init, &max_init_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {

        printf("\nInitialization Time (Max): %f seconds\n", max_init_time);
        printf("Communication Time (Max): %f seconds\n", max_comm_time);
        printf("Computation Time (Max): %f seconds\n", max_comp_time);

        printf("\n%f %f %f %s\n", max_init_time, max_comm_time, max_comp_time, getenv("SLURM_NNODES"));
    }
    
    
    
    free(MAT_LOC);
    free(MAT_LOC_NEW);




 
    
   

    MPI_Finalize();

    return 0;
}
