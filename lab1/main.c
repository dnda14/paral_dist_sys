#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX 1000   

double A[MAX][MAX];
double x[MAX], y[MAX];

void init() {
    for (int i = 0; i < MAX; i++) {
        x[i] = 1.0;  
        y[i] = 0.0;
        for (int j = 0; j < MAX; j++) 
            A[i][j] = (double)(i + j); 
        
    }
}

double gtiempo() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

int main() {
    init();

    double ini = gtiempo();
    for (int i = 0; i < MAX; i++) {
        for (int j = 0; j < MAX; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
    double t1 = gtiempo() - ini;
    printf("Tiempo primer bucle, fila: %f segundos\n", t1);

    for (int i = 0; i < MAX; i++) y[i] = 0.0;

    ini = gtiempo();
    for (int j = 0; j < MAX; j++) {
        for (int i = 0; i < MAX; i++) {
            y[i] += A[i][j] * x[j];
        }
    }
    double t2 = gtiempo() - ini;
    printf("Tiempo segundo bucle ,col: %f segundos\n", t2);

    return 0;
}
