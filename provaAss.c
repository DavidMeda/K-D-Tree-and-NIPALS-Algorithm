#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

// extern void euc_dist(float *, int, float *, int, int, float *);
// extern void calcolaTAss(float *vect, int numEle, float *result);
// extern void prodottoMatriceAss(float *ds, float *V, float *U, int i, int k);
// extern void aggiornaDatasetAss(float *, float *, float *, int, int);
// extern void dividiAss(float *, int, float);
// extern void calcolaSqrt(float, float *);

int main(int argc, char const *argv[])
{
    int k = 34;
    int n = 2;
    float *a = _mm_malloc(sizeof(float) * n * k, 16);
    float *u = _mm_malloc(sizeof(float) * n, 16);

    float *v = _mm_malloc(sizeof(float) * k, 16);

    float res = 0, res2 = 0;
    for (int i = 0; i < n * k; i++)
    {
        a[i] = i;
        if (i < 34)
            res += i;
        else
            res2 += i;
    }

    for (int i = 0; i < k; i++)
    {
        v[i] = 1;
    }

    u[0] = 2;
    u[1] = 3;

    aggiornaDatasetAss(a, u, v, 0, k);
    aggiornaDatasetAss(a, u, v, 1, k);

    printf("\n");

    for (int i = 0; i < n * k; i++)
    {
        printf(" %f, ", a[i]);
        if (i  == k-1)
            printf("\n");
    }
    printf("\n");

    return 0;
}
