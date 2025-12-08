#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    int n = (argc > 1) ? atoi(argv[1]) : 2048;
    int nrhs = 1, info;
    int *ipiv = (int *)malloc(n * sizeof(int));
    double *A = (double *)malloc(n * n * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));

   
    for (int i = 0; i < n * n; i++) A[i] = (double)rand() / RAND_MAX;
    for (int i = 0; i < n; i++) b[i] = (double)rand() / RAND_MAX;

    double start = omp_get_wtime();
    
    
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, 1);
    
    double end = omp_get_wtime();
    printf("MKL_Ref Size %d: %.4f seconds (%d ms)\n", n, end - start, (int)((end-start)*1000));

    free(A); free(b); free(ipiv);
    return 0;
}
