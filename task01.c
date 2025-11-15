#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 3
double augmented_matrix[N][N + 1] = {
    {2.0, 1.0, -1.0, 8.0},
    {-3.0, -1.0, 2.0, -11.0},
    {-2.0, 1.0, 2.0, -3.0}
};

double solution[N];

void gaussian_elimination() {
    printf("--- 1. Forward Elimination ---\n");
    int i, j, k;
    double factor;
    for (k = 0; k < N - 1; k++) {
        int max_row = k;
        for (i = k + 1; i < N; i++) {
            if (fabs(augmented_matrix[i][k]) > fabs(augmented_matrix[max_row][k])) {
                max_row = i;
            }
        }
        if (max_row != k) {
            for (j = k; j < N + 1; j++) {
                double temp = augmented_matrix[k][j];
                augmented_matrix[k][j] = augmented_matrix[max_row][j];
                augmented_matrix[max_row][j] = temp;
            }
            printf(" Row %d swapped with Row %d\n", k, max_row);
        }
        if (fabs(augmented_matrix[k][k]) < 1e-9) {
            printf("\n!!! Singular Matrix! Cannot solve! !!!\n");
            return;
        }
        for (i = k + 1; i < N; i++) {
            factor = augmented_matrix[i][k] / augmented_matrix[k][k];

            for (j = k; j < N + 1; j++) {
                augmented_matrix[i][j] = augmented_matrix[i][j] - factor * augmented_matrix[k][j];
            }
        }
    }
}

void back_substitution() {
    printf("\n--- 2. Back Substitution ---\n");
    double sum;

    for (int i = N - 1; i >= 0; i--) {
        sum = 0.0;
        
        for (int j = i + 1; j < N; j++) {
            sum = sum + augmented_matrix[i][j] * solution[j];
        }

        solution[i] = (augmented_matrix[i][N] - sum) / augmented_matrix[i][i];
    }
}

void print_results() {
    printf("\n--- 3. Solution Vector (X) ---\n");
    for (int i = 0; i < N; i++) {
        char var = 'X'; 
        if (i == 0) var = 'x';
        if (i == 1) var = 'y';
        if (i == 2) var = 'z';

        printf("Variable %c: %.4f\n", var, solution[i]);
    }
}

int main() {
    
    printf("Gaussian Elimination Solver for A * X = B\n");
    printf("Expected solution for 3x3 example: (x=2, y=3, z=-1)\n\n");
    printf("Initial Augmented Matrix [A | B]:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++) {
            printf("%.2f\t", augmented_matrix[i][j]);
        }
        printf("\n");
    }
    gaussian_elimination();
    printf("\nMatrix after Elimination [U | C] :\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N + 1; j++) {
            printf("%.4f\t", augmented_matrix[i][j]);
        }
        printf("\n");
    }

    back_substitution();

    print_results();

    return 0;
}