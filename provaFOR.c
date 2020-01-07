#include <stdlib.h>
#include <stdio.h>
#include <xmmintrin.h>

#define BLOCKSIZE 32

extern void forCount(float *A, int i, int n);

void provaFOR(float *A, int n)
{
    int i;
    for (i = 0; i < n; i += BLOCKSIZE)
        forCount(A, i, n);
}

int main(int argc, char const *argv[])
{
    int n = 100;
    float *A = _mm_malloc(n * sizeof(float), 16);
    provaFOR(A, n);

    printf("Stampo vettore");
    for (int i = 0; i < n; i++)
    {
        printf("%f ", A[i]);
    }

    return 0;
}
