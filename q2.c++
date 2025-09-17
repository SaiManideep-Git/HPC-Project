#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Example function f(i)
int f(int i) {
    return i * i;   // square of i
}

int main(int argc, char *argv[]) {
    int n = 5000;  // size of loop
    int *v = (int *)calloc(n, sizeof(int));
    int *indices = (int *)malloc(n * sizeof(int));

    // Initialize indices (identity mapping)
    for (int i = 0; i < n; i++) {
        indices[i] = i;
    }

    // Use 12 threads
    omp_set_num_threads(12);

    // Parallel loop with critical section
    #pragma omp parallel for default(none) shared(v, indices, n)
    for (int i = 0; i < n; i++) {
        #pragma omp critical
        v[indices[i]] += f(i);
    }

    // Print results
    printf("Result vector (using critical, 12 threads):\n");
    for (int i = 0; i < n; i++) {
        printf("v[%d] = %d\n", i, v[i]);
    }

    free(v);
    free(indices);
    return 0;
}
