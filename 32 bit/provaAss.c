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

// int main(int argc, char const *argv[])
// {
//     int k = 34;
//     int n = 2;
//     float *a = _mm_malloc(sizeof(float) * n * k, 16);
//     float *u = _mm_malloc(sizeof(float) * n, 16);

//     float *v = _mm_malloc(sizeof(float) * k, 16);

//     float res = 0, res2 = 0;
//     for (int i = 0; i < n * k; i++)
//     {
//         a[i] = i;
//         if (i < 34)
//             res += i;
//         else
//             res2 += i;
//     }

//     for (int i = 0; i < k; i++)
//     {
//         v[i] = 1;
//     }

//     u[0] = 2;
//     u[1] = 3;

//     aggiornaDatasetAss(a, u, v, 0, k);
//     aggiornaDatasetAss(a, u, v, 1, k);

//     printf("\n");

//     for (int i = 0; i < n * k; i++)
//     {
//         printf(" %f, ", a[i]);
//         if (i  == k-1)
//             printf("\n");
//     }
//     printf("\n");

//     return 0;
// }
// int main(int argc, char const *argv[])
// {
//     int **a = (int *)malloc(sizeof(int *) * 3);

//     // for (int i = 0; i < 3; i++)
//     // // {
//     // a[0] = malloc(sizeof(int) * 5);
//     // a[1] = malloc(sizeof(int) * 3);
//     // a[2] = malloc(sizeof(int) * 6);
//     // for (int j = 0; j < 5; j++)
//     // {
//     //     a[0][j] = j + 1;
//     // }

//     // for (int j = 0; j < 3; j++)
//     // {
//     //     a[1][j] = j + 1;
//     // }

//     // for (int j = 0; j < 6; j++)
//     // {
//     //     a[2][j] = j + 1;
//     // }
//     // // }

//     // printf("\n");

//     // for (int j = 0; j < 5; j++)
//     // {
//     //     printf(" %d ", a[0][j]);
//     // }
//     // printf("\n");

//     // for (int j = 0; j < 3; j++)
//     // {
//     //     printf(" %d ", a[1][j]);
//     // }
//     // printf("\n");
//     // for (int j = 0; j < 6; j++)
//     // {
//     //     printf(" %d ", a[2][j]);
//     // }

//     printf(" \n%ld \n", sizeof(int));

//     return 0;
// }

void prodottoMatriceTrasp(float *v, float *ds, float *u, int numEleU, int k)
{
    int i, j;
    float sum = 0;
    for (i = 0; i < k; i++)
    {
        sum = 0;
        for (j = 0; j < numEleU; j++)
        {
            sum += ds[j * k + i] * u[j];
        }
        v[i] = sum;
    }
}

extern void prodMatriceTrasAss(float *, float *, float *, int, int);
int main(int argc, char const *argv[])
{

    int k = 64;
    int n = 53;
    float *a = malloc(sizeof(float) * n * k);
    float *u = malloc(sizeof(float) * n);
    float *v = malloc(sizeof(float) * k);

      float *v1 = malloc(sizeof(float) * k);
    // printf("dataset\n");

    for (int i = 0; i < n * k; i++)
    {
        a[i] = i + 1;
        if (i == k - 1)
        {
            // printf("\n");
        }
        // printf(" %f, ", a[i]);
    }
    // printf("matrice u\n");

    for (int i = 0; i < n; i++)
    {
        u[i] = i + 10;
        // printf(" %f, ", u[i]);
    }

    prodottoMatriceTrasp(v, a, u, n, k);

    printf("\nmatrice risultato C u\n");

    for (int i = 0; i < k; i++)
    {
        printf(" %f, ", v[i]);
    }

    // memset(v,0, sizeof(float)*k);

    prodMatriceTrasAss(a, v1, u, n, k);

    printf("\nmatrice risultato assembly u\n");

    for (int i = 0; i < k; i++)
    {
        printf(" %f, ", v1[i]);
    }

    printf("\n\n\n");


    prodottoMatriceTrasp(v, a, u, n, k);

    printf("\nmatrice risultato C u\n");

    for (int i = 0; i < k; i++)
    {
        printf(" %f, ", v[i]);
    }

    memset(v1,0, sizeof(float)*k);

    prodMatriceTrasAss(a, v1, u, n, k);

    printf("\nmatrice risultato assembly u\n");

    for (int i = 0; i < k; i++)
    {
        printf(" %f, ", v1[i]);
    }

    return 0;
}
