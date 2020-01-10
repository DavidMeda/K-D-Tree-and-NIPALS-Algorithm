#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

extern void multi(float *ds, float *V, float *U, int cut, int i, int k, int h);

int main(int argc, char const *argv[])
{

    int n = 1;
    int k = 4;
    int h = 2;
    int cut = 0;
    float *ds = _mm_malloc(sizeof(float) * n * k, 16);
    float *V = _mm_malloc(sizeof(float) * k * h, 16);
    float *U = _mm_malloc(sizeof(float) * n * h, 16);

    printf("\nds: [");
    for (int i = 0; i < n * k; i++)
    {
        ds[i] = i + 1;
        printf(" %f ", ds[i]);
    }
    printf("]\n V: [");
    for (int i = 0; i < h * k; i++)
    {
        V[i] = i + 10;
        printf(" %f ", V[i]);
    }
    printf(" ] ");

    multi(ds, V, U, cut, 0, k, h);

    for (int i = 0; i < n * h; i++)
    {
        printf("\nU[%d]: %f", i, U[i]);
    }
    printf("\n");
    /* code */
    return 0;
}
