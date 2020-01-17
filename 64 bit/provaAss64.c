#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

extern void dividiAss64(float *a, int n, float *value);

int main(int argc, char const *argv[])
{
    int n = 34;
    float *a = malloc(sizeof(float) * n);
    float res = 2;
    for (int i = 0; i < n; i++)
    {
        a[i] = 4;
    }

    dividiAss64(a, n, &res);

    printf("\n");
    for (int i = 0; i < n; i++)
    {
        printf(" %f ", a[i]);
    }
    printf("\n");
    return 0;
}
