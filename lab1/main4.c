#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

double gtiempo() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

double *gmem(int N) {
    double *m = (double*)aligned_alloc(64, sizeof(double) * N * N);
    if (!m) { perror("alloc"); exit(EXIT_FAILURE); }
    return m;
}

void mul_cla(const double *A, const double *B, double *C, int N) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k)
                sum += A[i*N + k] * B[k*N + j];
            C[i*N + j] = sum;
        }
}

void mul_blo(const double *A, const double *B, double *C, int N, int bs) {
    for (int i = 0; i < N*N; ++i) C[i] = 0.0;

    for (int ii = 0; ii < N; ii += bs)
        for (int jj = 0; jj < N; jj += bs)
            for (int kk = 0; kk < N; kk += bs) {
                int i_max = (ii + bs < N) ? ii + bs : N;
                int j_max = (jj + bs < N) ? jj + bs : N;
                int k_max = (kk + bs < N) ? kk + bs : N;
                for (int i = ii; i < i_max; ++i)
                    for (int k = kk; k < k_max; ++k) {
                        double a_ik = A[i*N + k];
                        for (int j = jj; j < j_max; ++j)
                            C[i*N + j] += a_ik * B[k*N + j];
                    }
            }
}

void gen_rd(double *A, int N, unsigned seed) {
    srand(seed);
    for (int i = 0; i < N*N; ++i) A[i] = (double)(rand() % 10);
}

double prom(const double *A, const double *B, double *C, int N, int reps) {
    double total = 0.0;
    for (int r=0; r<reps; r++) {
        double t0 = gtiempo();
        mul_cla(A,B,C,N);
        total += gtiempo() - t0;
    }
    return total / reps;
}

double promblo(const double *A, const double *B, double *C, int N, int bs, int reps) {
    double total = 0.0;
    for (int r=0; r<reps; r++) {
        double t0 = gtiempo();
        mul_blo(A,B,C,N,bs);
        total += gtiempo() - t0;
    }
    return total / reps;
}

int main() {
    int sizes[] = {128, 256, 512, 1024}; 
    int blocktam[] = {8, 16, 32, 64};  
    int reps = 3; 

    printf("N,bs,tcla,t_block\n");

    for (int si = 0; si < sizeof(sizes)/sizeof(sizes[0]); ++si) {
        int N = sizes[si];

        double *A = gmem(N);
        double *B = gmem(N);
        double *C = gmem(N);

        gen_rd(A, N, 1234);
        gen_rd(B, N, 4321);

        double tcla = prom(A,B,C,N,reps);

        for (int bi = 0; bi < sizeof(blocktam)/sizeof(blocktam[0]); ++bi) {
            int bs = blocktam[bi];
            double t_block = promblo(A,B,C,N,bs,reps);
            printf("%d,%d,%.6f,%.6f\n", N, bs, tcla, t_block);
        }

        free(A); free(B); free(C);
    }
    return 0;
}
