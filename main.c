#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "timer.h"
#include "dgesv.h"

double *generate_matrix(unsigned int size, unsigned int seed) {
    double *matrix = (double *) malloc(sizeof(double) * size * size);
    srand(seed);
    for (unsigned int i = 0; i < size * size; i++)
        matrix[i] = (double)(rand() % 1000) / 10.0 + 0.1;

    for (unsigned int i = 0; i < size; i++) {
        double sum = 0.0;
        for (unsigned int j = 0; j < size; j++)
            if (i != j) sum += fabs(matrix[i * size + j]);
        matrix[i * size + i] = sum + 100.0;
    }
    return matrix;
}

double *generate_vector(unsigned int size, unsigned int seed) {
    double *vector = (double *) malloc(sizeof(double) * size);
    srand(seed);
    for (unsigned int i = 0; i < size; i++)
        vector[i] = (double)(rand() % 1000) / 10.0 + 0.1;
    return vector;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <matrix_size>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int nrhs = 1;

    double *A = generate_matrix(N, 1);
    double *B = generate_vector(N, 2);

 
    
    timeinfo start, end;
    timestamp(&start);

   
    int info = my_dgesv(N, nrhs, A, B);
    
    timestamp(&end);

    if (info != 0) {
        printf("Error: matrix is singular\n");
    } else {
        printf("Solved linear system of size %d\n", N);
        printf("Time taken by my_dgesv: %ld ms\n", diff_milli(&start, &end));
    }

    free(A);
    free(B);
    return info;
}