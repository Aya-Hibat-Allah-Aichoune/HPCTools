#include "dgesv.h"
#include <stdlib.h>
#include <math.h>

int my_dgesv(int n, int nrhs, double * restrict a, double * restrict b)
{
    int i, j, k, max_row;
    double tmp;

    for (k = 0; k < n; k++) {
        
        max_row = k;
        for (i = k+1; i < n; i++)
            if (fabs(a[i*n + k]) > fabs(a[max_row*n + k]))
                max_row = i;

        if (max_row != k) {
            for (j = 0; j < n; j++) {
                tmp = a[k*n + j];
                a[k*n + j] = a[max_row*n + j];
                a[max_row*n + j] = tmp;
            }
            for (j = 0; j < nrhs; j++) {
                tmp = b[k*nrhs + j];
                b[k*nrhs + j] = b[max_row*nrhs + j];
                b[max_row*nrhs + j] = tmp;
            }
        }

        double *row_k = &a[k*n];

        for (i = k+1; i < n; i++) {
            double factor = a[i*n + k] / row_k[k];
            double *row_i = &a[i*n];

            row_i[k] = 0.0;

            #pragma GCC ivdep 
            for (j = k+1; j < n; j++) {
                row_i[j] -= factor * row_k[j];
            }
        }
        
        for (i = k+1; i < n; i++) {
            double factor = a[i*n + k] / a[k*n + k];
            #pragma GCC ivdep
            for (j = 0; j < nrhs; j++)
                b[i*nrhs + j] -= factor * b[k*nrhs + j];
        }
    }

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
