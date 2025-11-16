#include "dgesv.h"
#include <stdlib.h>
#include <math.h>

int my_dgesv(int n, int nrhs, double *a, double *b)
{
    int i, j, k, max_row;
    double tmp;

    // Gaussian elimination with partial pivoting
    for (k = 0; k < n; k++) {
        // Find pivot
        max_row = k;
        for (i = k+1; i < n; i++)
            if (fabs(a[i*n + k]) > fabs(a[max_row*n + k]))
                max_row = i;

        // Swap rows in A
        if (max_row != k) {
            for (j = 0; j < n; j++) {
                tmp = a[k*n + j];
                a[k*n + j] = a[max_row*n + j];
                a[max_row*n + j] = tmp;
            }
            // Swap rows in B
            for (j = 0; j < nrhs; j++) {
                tmp = b[k*nrhs + j];
                b[k*nrhs + j] = b[max_row*nrhs + j];
                b[max_row*nrhs + j] = tmp;
            }
        }

        // Eliminate below pivot
        for (i = k+1; i < n; i++) {
            double factor = a[i*n + k] / a[k*n + k];
            a[i*n + k] = 0.0;
            for (j = k+1; j < n; j++)
                a[i*n + j] -= factor * a[k*n + j];
            for (j = 0; j < nrhs; j++)
                b[i*nrhs + j] -= factor * b[k*nrhs + j];
        }
    }

    // Back substitution
    for (j = 0; j < nrhs; j++) {
        for (i = n-1; i >= 0; i--) {
            tmp = b[i*nrhs + j];
            for (k = i+1; k < n; k++)
                tmp -= a[i*n + k] * b[k*nrhs + j];
            b[i*nrhs + j] = tmp / a[i*n + i];
        }
    }

    return 0;
}
