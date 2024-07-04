#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    FILE *file;
    int rows, cols;
    double **matrix;
    
    rows = ROWS;
    cols = COLS;
    // Open the binary file for reading
    file = fopen("result.bin", "rb");
    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }


    // Allocate memory for the matrix
    matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));
    }

    // Read the matrix data
    for (int i = 0; i < rows; i++) {
        fread(matrix[i], sizeof(double), cols, file);
    }

    // Close the file
    fclose(file);

    // Print the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }

    // Free the memory
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}