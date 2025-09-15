#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void multi(double **A, double **B, double **C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}

int main() {
    int tam[] = {100, 200, 500, 1000};
    int ntma = sizeof(tam) / sizeof(tam[0]);

    for (int s = 0; s < ntma; s++) {
        int N = tam[s];
        printf("\nMultiplicación de matrices %dx%d\n", N, N);

        double **A = (double **)malloc(N * sizeof(double *));
        double **B = (double **)malloc(N * sizeof(double *));
        double **C = (double **)malloc(N * sizeof(double *));
        for (int i = 0; i < N; i++) {
            A[i] = (double *)malloc(N * sizeof(double));
            B[i] = (double *)malloc(N * sizeof(double));
            C[i] = (double *)malloc(N * sizeof(double));
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = (i + j) % 10;
                B[i][j] = (i - j) % 10;
                C[i][j] = 0.0;
            }
        }

        double ini = get_time();
        multi(A, B, C, N);
        double tot = get_time() - ini;

        printf("Tiempo: %f segundos\n", tot);

        for (int i = 0; i < N; i++) {
            free(A[i]);
            free(B[i]);
            free(C[i]);
        }
        free(A);
        free(B);
        free(C);
    }

    return 0;
}
                                                                                                                                                                                                                                                               