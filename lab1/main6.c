#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

long long mem_reads = 0;
long long mem_writes = 0;

double ftiempo() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

double *gmem(int N) {
    double *m = (double*)aligned_alloc(64, sizeof(double) * N * N);
    if (!m) { perror("alloc"); exit(EXIT_FAILURE); }
    return m;
}

void mult_cla(const double *A, const double *B, double *C, int N) {
    printf("=== ALGORITMO CLASICO - ANALISIS DE ACCESOS ===\n");
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k) {
                // Contar accesos
                mem_reads += 2;  // A[i*N + k] + B[k*N + j]
                sum += A[i*N + k] * B[k*N + j];
                
                // Mostrar patrón para matrices pequeñas
                if (N <= 8) {
                    printf("i=%d, j=%d, k=%d: A[%d], B[%d]\n", 
                           i, j, k, i*N + k, k*N + j);
                }
            }
            mem_writes += 1;  // C[i*N + j] = sum
            C[i*N + j] = sum;
        }
    }
}

void mul_blo(const double *A, const double *B, double *C, int N, int bs) {
    printf("=== ALGORITMO POR BLOQUES - ANALISIS DE ACCESOS ===\n");
    printf("Tamaño de bloque: %d x %d\n", bs, bs);
    
    for (int i = 0; i < N*N; ++i) {
        C[i] = 0.0;
        mem_writes++;
    }
    
    int block_cont = 0;
    
    for (int ii = 0; ii < N; ii += bs) {
        for (int jj = 0; jj < N; jj += bs) {
            for (int kk = 0; kk < N; kk += bs) {
                block_cont++;
                
                int i_max = (ii + bs < N) ? ii + bs : N;
                int j_max = (jj + bs < N) ? jj + bs : N;
                int k_max = (kk + bs < N) ? kk + bs : N;
                
                printf("Bloque %d: ii=%d-%d, jj=%d-%d, kk=%d-%d\n", 
                       block_cont, ii, i_max-1, jj, j_max-1, kk, k_max-1);
                
                for (int i = ii; i < i_max; ++i) {
                    for (int k = kk; k < k_max; ++k) {
                        double a_ik = A[i*N + k];
                        mem_reads++;  // A[i*N + k]
                        
                        for (int j = jj; j < j_max; ++j) {
                            mem_reads += 2;   // B[k*N + j] + C[i*N + j]
                            mem_writes += 1;  // C[i*N + j] +=
                            C[i*N + j] += a_ik * B[k*N + j];
                        }
                    }
                }
            }
        }
    }
    
    printf("Total de bloques procesados: %d\n", block_cont);
}

void gen_rd(double *A, int N, unsigned seed) {
    srand(seed);
    for (int i = 0; i < N*N; ++i) {
        A[i] = (double)(rand() % 10);
        mem_writes++;
    }
}

void revisar_cache(int N, int bs) {
    printf("\n=== ANALISIS TEORICO DE CACHE ===\n");
    printf("Tamaño de matriz: %d x %d\n", N, N);
    printf("Elementos totales: %d\n", N*N);
    printf("Memoria por matriz: %.2f KB\n", (N*N*sizeof(double))/1024.0);
    printf("Memoria total (3 matrices): %.2f KB\n", (3*N*N*sizeof(double))/1024.0);
    
    // Asumiendo cache L1 típica de 32KB
    int cache_size = 32 * 1024;  // 32KB
    int elements_in_cache = cache_size / sizeof(double);  // ~4096 elementos
    
    printf("\nCache L1 típica: 32KB (~%d elementos double)\n", elements_in_cache);
    
    if (N*N > elements_in_cache) {
        printf("⚠️  Las matrices NO caben completamente en L1\n");
    } else {
        printf("✅ Las matrices caben en L1\n");
    }
    
    printf("\nAlgoritmo por bloques - Tamaño de bloque: %d x %d\n", bs, bs);
    int elements_per_block = bs * bs;
    printf("Elementos por bloque: %d\n", elements_per_block);
    
    if (3 * elements_per_block < elements_in_cache) {
        printf("✅ Los 3 bloques (A, B, C) caben juntos en L1\n");
    } else {
        printf("⚠️  Los 3 bloques NO caben juntos en L1\n");
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <classic|blocked> <N> [bs]\n", argv[0]);
        return 1;
    }
    
    const char *mode = argv[1];
    int N = atoi(argv[2]);
    int bs = (argc > 3) ? atoi(argv[3]) : 32;
    
    if (N > 16) {
        printf("Para análisis detallado, usa N <= 16\n");
        printf("Ejecutando análisis teórico solamente...\n");
        revisar_cache(N, bs);
        return 0;
    }
    
    double *A = gmem(N);
    double *B = gmem(N);
    double *C = gmem(N);
    
    gen_rd(A, N, 1234);
    gen_rd(B, N, 4321);
    
    mem_reads = 0;
    mem_writes = 0;
    
    printf("Matrices %d x %d inicializadas\n", N, N);
    revisar_cache(N, bs);
    
    double t0 = ftiempo();
    
    if (strcmp(mode, "classic") == 0) {
        mult_cla(A, B, C, N);
    } else {
        mul_blo(A, B, C, N, bs);
    }
    
    double tiempo = ftiempo() - t0;
    
    printf("\n=== ESTADISTICAS ===\n");
    printf("Tiempo de ejecución: %.6f segundos\n", tiempo);
    printf("Lecturas de memoria: %lld\n", mem_reads);
    printf("Escrituras de memoria: %lld\n", mem_writes);
    printf("Total accesos: %lld\n", mem_reads + mem_writes);
    
    long long operaciones_teoricas = 2LL * N * N * N; // multiply + add por cada elemento
    printf("Operaciones aritméticas teóricas: %lld\n", operaciones_teoricas);
    printf("Ratio accesos/operaciones: %.2f\n", 
           (double)(mem_reads + mem_writes) / operaciones_teoricas);
    
    free(A); 
    free(B); 
    free(C);
    return 0;
}