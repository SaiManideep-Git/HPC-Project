#include <stdio.h>
#include <omp.h>

int main() {
    long long num_steps = 10000000000;   // number of intervals (higher = more accurate)
    double step = 1.0 / (double) num_steps;
    double pi = 0.0;

    int num_threads = 12;   // use 12 threads
    omp_set_num_threads(num_threads);

    double start = omp_get_wtime();

    #pragma omp parallel
    {
        double x, sum = 0.0;
        int tid = omp_get_thread_num();

        // Parallel for loop with reduction to avoid race conditions
        #pragma omp for
        for (long long i = 0; i < num_steps; i++) {
            x = (i + 0.5) * step;   // midpoint rule
            sum += 4.0 / (1.0 + x * x);
        }

        // Accumulate results safely
        #pragma omp atomic
        pi += sum;
    }

    pi *= step;

    double end = omp_get_wtime();

    printf("Approximate value of PI = %.15f\n", pi);
    printf("Time taken with %d threads = %f seconds\n", num_threads, end - start);

    return 0;
}
