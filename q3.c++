#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Example function f(i)
int f(int i) {
    return i * i;   // square of i
}

int main() {
    int n = 1000000;  // bigger size to see timing differences
    int *v1 = (int *)calloc(n, sizeof(int));  // for critical version
    int *v2 = (int *)calloc(n, sizeof(int));  // for atomic version
    int *indices = (int *)malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) {
        indices[i] = i % 1000;  // force collisions on indices
    }

    omp_set_num_threads(12);

    double start, end;

    // ---------------- Critical Version ----------------
    start = omp_get_wtime();
    #pragma omp parallel for default(none) shared(v1, indices, n)
    for (int i = 0; i < n; i++) {
        #pragma omp critical
        v1[indices[i]] += f(i);
    }
    end = omp_get_wtime();
    printf("Time with critical: %f seconds\n", end - start);

    // ---------------- Atomic Version ----------------
    start = omp_get_wtime();
    #pragma omp parallel for default(none) shared(v2, indices, n)
    for (int i = 0; i < n; i++) {
        #pragma omp atomic
        v2[indices[i]] += f(i);
    }
    end = omp_get_wtime();
    printf("Time with atomic:   %f seconds\n", end - start);

    // Check correctness (first 10 values)
    printf("\nSample results (index : v1 vs v2):\n");
    for (int i = 0; i < 10; i++) {
        printf("%d : %d vs %d\n", i, v1[i], v2[i]);
    }

    free(v1);
    free(v2);
    free(indices);
    return 0;
}
