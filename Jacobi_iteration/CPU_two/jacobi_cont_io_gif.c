#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 62
#define NUM_ITER 2000


int main(int argc, char **argv)
{
    
    MPI_Init(&argc, &argv);
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N_LOC = N / size ;
    int rest = N % size;
    int offset = 0;
    
    if (rank != 0 && rank != size - 1) {N_LOC += 2;}
    else {N_LOC += 1;}


    if (rank < rest) {N_LOC++;} else {offset = rest;}
    
    double *MAT_LOC = (double *)malloc(N_LOC * N * sizeof(double));
    double *MAT_LOC_NEW = (double *)malloc(N_LOC * N * sizeof(double));
    
    //// INNER INITIALIZATION ////
    for (int i = 0; i < N_LOC; i++)
    {
        for (int j = 0; j < N; j++)
        {
            MAT_LOC[i * N + j] = 0.5;
            MAT_LOC_NEW[i * N + j] = 0.5;
        }
    }

    // TOP ROW = 0
    if (rank == 0){
    
   
    for (int i = 0; i < N; i++)
        {
            MAT_LOC[i] = 0;
            MAT_LOC_NEW[i] = 0;
        }
    
    }

    
    

  

    //LINEARY SPACED BOTTOM COLUMN
    if (rank == size-1){
    
    for (int i = 0; i < N; i++)
    {   double step = (100.0 - 0) / (N - 1);
        MAT_LOC[(N_LOC - 1) * N + i] = 100 - i * step;
        MAT_LOC_NEW[(N_LOC - 1) * N + i] = 100 - i * step;
    }
    }

    
    //LINEARY SPACED LEFT COLUMN (BY RANKS)
    
    for (int i = 1; i < N_LOC-1; i++)
    {   
        int global_idx = (N_LOC-2)*rank+(i-1)+offset;
        
        if (rank == 0){
            global_idx = (N_LOC)*rank+(i)+offset;
        } 
        else if (rank == size-1){
            global_idx = (N_LOC-2)*rank+(i+1)+offset;
        }

        double step = (100.0 - 0) / (N  - 1);
        MAT_LOC[i*N] = global_idx  * step ;
        MAT_LOC_NEW[i*N] = global_idx  * step ;
    }


    //RIGHT COLUMN = 0 (BY RANKS)
    for (int i = 0; i < N_LOC; i++)
    {
        MAT_LOC[i*N + N - 1] = 0;
        MAT_LOC_NEW[i*N + N - 1] = 0;
    }

    ///////////////////////////
   

 
    
  




    





for (int iter= 0; iter <= NUM_ITER; iter++)

{   
     //// COMMUNICATION ////
    
    

    
    
    double * ghost_up = MAT_LOC;
    double * ghost_down = MAT_LOC + (N_LOC-1) * N;

    
    double* first_row_point = MAT_LOC + N;
    double* last_row_point = MAT_LOC + (N_LOC-2) * N;
   
    
  

   
    

    
   
    
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


    
    ///////////////////////////


    ///// COMPUTATION /////
    
    for (int i = 1; i < N_LOC-1; i++)
    {
        for (int j = 1; j < N-1; j++)
        {   
            MAT_LOC_NEW[i * N + j] = 0.25 * (MAT_LOC[(i - 1) * N + j] + MAT_LOC[(i + 1) * N + j] + MAT_LOC[i * N + j - 1] + MAT_LOC[i * N + j + 1]);
        }
    }

    ///////////////////////////
    
    
    // SWAP ARRAYS USING POINTERS
    double *temp = MAT_LOC;
    MAT_LOC = MAT_LOC_NEW;
    MAT_LOC_NEW = temp;




    

   
    


 
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
  
    if (iter % 40 == 0)
    {
    // Open the file
    char filename[100];
    sprintf(filename, "./mat_files/matrix_%d.bin", iter);
    
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);


    //file seek
    MPI_File_set_view(file, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

    // Write the local portion of the matrix to the file
    MPI_File_write(file, final_pointer, FINAL_N_LOC * N, MPI_DOUBLE, &status);

    // Close the file
    MPI_File_close(&file);
    }
    
}

    
    





 
    
   

    MPI_Finalize();

    return 0;
}
