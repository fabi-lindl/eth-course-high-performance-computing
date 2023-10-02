#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

const int MIN_N = 1000 / sizeof(int);      // From 1 KB
const int MAX_N = 20000000 / sizeof(int);  // to 20 MB.
const int NUM_SAMPLES = 100;
const int M = 100000000;    // Operations per sample.
int a[MAX_N];               // Permutation array.


void sattolo(int *p, int N) {
    /*
     * Generate a random single-cycle permutation using Satollo's algorithm.
     *
     * https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#Sattolo's_algorithm
     * https://danluu.com/sattolo/
     */
    for (int i = 0; i < N; ++i)
        p[i] = i;
    for (int i = 0; i < N - 1; ++i)
        std::swap(p[i], p[i + 1 + rand() % (N - i - 1)]);
}

double measure(int N, int mode) {

    // Stop execution time of this function. 
    double times = 0;
    auto tstart = std::chrono::system_clock::now();

    // Initialize pointer of the array. 
    int arr[N];
    int *ptr = &arr[0];

    if (mode == 0) {
        // TODO: Question 1b: Use the sattolo function to generate a random one-cycle permutation.
        
        // Pass the provided function the pointer and N as arguments. 
        sattolo(ptr, N);

    } else if (mode == 1) {
        // TODO: Question 1c: Initialize the permutation such that k jumps by 1 item every step (cyclically).
        
        // The array needs to contain elements in ascending order, where elmt_k + 1 == elmt_k+1. 
        // The first element of the array needs to be int 1. Zero would not allow the first jump. 
        for (int i = 0; i < N; i++) {
            ptr[i] = i + 1;
        }
        // The last element of the array needs to be zero so that the traversing starts again. 
        ptr[N-1] = 0;

    } else if (mode == 2) {
        // TODO: Question 1d: Initialize the permutation such that k jumps by 64 bytes (cyclically).

        // The first array element needs to be 16 so that a jump by 16 elements (64 bytes / 4 bytes (=int))
        // can be guaranteed. All other elements are then incremented by 1. 
        for (int i = 0; i < N; i++) {
            ptr[i] = i + 16;
        }        
    }

    // Print the array. 
    // printf("Array setup\n");
    // for (int i = 0; i < N; i++) {
    //     printf("%d  : %d\n", i, arr[i]);
    // }

    // TODO: Question 1b: Traverse the list (make M jumps, starting from k = 0) and measure the execution time.
    
    // Number of steps. 
    int step = *ptr;
    // printf("\nJump values\n");
    if (mode == 0 || mode == 1) {
        // printf("---------------------------\n");    
        for (int i = 0; i < M; i++) {
            step = *(ptr+step);

            // if (i < 10) {
            //     printf("step: %d\n", step);
            // }

            // if (i > 240 && i < 260) {
            //     printf("step: %d\n", step);
            // }
            
        }
        // printf("---------------------------\n");
    }
    else {
        // Mode 3.
        // Max value is stored at array position N-1. 
        
        int maxVal = ptr[N-1]; 
        // printf("N = %d    maxVal = %d\n---------------------\n", N, maxVal);

        for (int i = 0; i < M; i++) {
            step = *(ptr+step);

            // if (i < 30) {
            //     printf("step: %d\n", step);
            // }

            // Check boundary. 
            if (step+16 > maxVal) {
                // printf("boundary check: step = %d\n", step);
                step = (step+16) % N;
                // printf("step = %d\n", step);
            }

        }
    }
    
    // TODO: Question 1b: Return execution time in seconds.
    
    auto tend = std::chrono::system_clock::now();
    times += std::chrono::duration<double>(tend-tstart).count();
    
    return times;
}

void run_mode(int mode) {
    /*
     * Run the measurement for many different values of N and output in a
     * format compatible with the plotting script.
     */
    // printf("%9s  %9s  %7s  %7s\n", "N", "size[kB]", "t[s]", "op_per_sec[10^9]");
    // printf("%9s  %9s  %7s\n", "N", "size[kB]", "op_per_sec[10^9]");
    printf("%9s    %7s\n", "N", "op_per_sec[10^9]");
    for (int i = 0; i < NUM_SAMPLES; ++i) {

    // for (int i = 0; i < 1; ++i) {

        // Generate N in a logarithmic scale.
        int N = (int)(MIN_N * std::pow((double)MAX_N / MIN_N,
                                       (double)i / (NUM_SAMPLES - 1)));
        double t = measure(N, mode);
        // printf("%9d  %9.1lf  %7.5lf  %7.6lf\n",
        //        N, N * sizeof(int) / 1024., t, M / t * 1e-9);

        // Print with size. 
        // printf("%9d  %9.1lf  %7.6lf\n",
        //        N, N * sizeof(int) / 1024., M / t * 1e-9);

        // Print only performance. 
        printf("%9d    %7.6lf\n",
               N, M / t * 1e-9);

        // // Print only the size. 
        // printf("%9.1lf\n", N * sizeof(int) / 1024.);

        fflush(stdout);
    }
    printf("Mode done!\n");
    printf("\n\n");
    printf("Mode done!\n");
}

int main() {
    // Question 1b:
    // printf("MODE 0\n");
    // run_mode(0);   // Random.

    // TODO: Enable for Question 1c:
    // printf("MODE 1\n");
    // run_mode(1);   // Sequential (jump by sizeof(int) bytes).

    // TODO: Enable for Question 1d:
    printf("MODE 2\n");
    run_mode(2);   // Sequential (jump by cache line size, i.e. 64 bytes).

    return 0;
}

