#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

extern void sommatoria(float *ds, int n, int k, int col);

int main(int argc, char const *argv[])
{
    float *ds = _mm_malloc(sizeof(float) * 16, 16);

    for(int i =0; i<16; i++){
        ds[i] = i+10;
    }
    int n = 8;
    int k = 2;
    int cut = 0;

    sommatoria(ds,n,k,cut);

        for(int i =0; i<16; i++){
            printf("\nvec[%d]: %f", i, ds[i]);
        }


    /* code */
    return 0;
}
