#include "dgesv.h"
#include <stdlib.h>
#include <math.h>

// Task 2: Optimized with restrict and Loop Fusion
int my_dgesv(int n, int nrhs, double * restrict a, double * restrict b)
{
    int i, j, k, max_row;
    double tmp;

    for (k = 0; k < n; k++) {

        // 1. Pivot search
        max_row = k;
        for (i = k+1; i < n; i++)
            if (fabs(a[i*n + k]) > fabs(a[max_row*n + k]))
                max_row = i;

        // 2. Swap rows if needed
        if (max_row != k) {
            // Swap row in A
            for (j = k; j < n; j++) { // Start from k (optimization)
                tmp = a[k*n + j];
                a[k*n + j] = a[max_row*n + j];
                a[max_row*n + j] = tmp;
            }
            // Swap row in B
            for (j = 0; j < nrhs; j++) {
                tmp = b[k*nrhs + j];
                b[k*nrhs + j] = b[max_row*nrhs + j];
                b[max_row*nrhs + j] = tmp;
            }
        }

        // 3. Gaussian Elimination (Fused Loop)
        // Hna drna optimization: hsebna factor merra whda
        // w khdmna 3la A w B f nafs l wa9t bach nstaghlo Cache
        
        double pivot = a[k*n + k];
        
        for (i = k+1; i < n; i++) {
            double factor = a[i*n + k] / pivot;
            
            // Update Row A (Vectorized)
            #pragma GCC ivdep
            for (j = k+1; j < n; j++) {
                a[i*n + j] -= factor * a[k*n + j];
            }
            
            // Update Row B (Vectorized)
            #pragma GCC ivdep
            for (j = 0; j < nrhs; j++) {
                b[i*nrhs + j] -= factor * b[k*nrhs + j];
            }
        }
    }

    // 4. Back Substitution
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
