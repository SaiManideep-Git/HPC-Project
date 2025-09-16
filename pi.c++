#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main() {
    long long num_points = 10000000;   // total random points
    long long inside_circle = 0;
    int num_threads = 12;

    omp_set_num_threads(num_threads);

    double start = omp_get_wtime();

    #pragma omp parallel
    {
        unsigned int seed = 1234 + omp_get_thread_num(); // thread-specific seed
        long long local_count = 0;

        #pragma omp for
        for (long long i = 0; i < num_points; i++) {
            double x = (double)rand_r(&seed) / RAND_MAX;
            double y = (double)rand_r(&seed) / RAND_MAX;

            if (x * x + y * y <= 1.0) {
                local_count++;
            }
        }

        #pragma omp atomic
        inside_circle += local_count;
    }

    double pi = 4.0 * (double)inside_circle / (double)num_points;

    double end = omp_get_wtime();

    printf("Approximate value of PI = %.15f\n", pi);
    printf("Time taken with %d threads = %f seconds\n", num_threads, end - start);

    return 0;
}
