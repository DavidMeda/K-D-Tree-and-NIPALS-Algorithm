// Compilare:
//		gcc -O0 -fopenmp -m32 -msse gemm5.o gemm5c.c -o gemm5c
//
// Test di velocità:
//		time ./gemm5c
//


#include <stdlib.h>
#include <stdio.h>
#include <xmmintrin.h>

#define TYPE float
#define BLOCKSIZE 32


extern void gemm5(float* A, float* B, float* C, int si, int sj, int sk, int n);


void gemm2(float* A, float* B, float* C, int i, int n)  
{
	int j, k;
	for (j = 0; j < n; j += BLOCKSIZE)
		for (k = 0; k < n; k += BLOCKSIZE){
			gemm5(A,B,C,i,j,k,n);
		}
}


void gemm1(float* A, float* B, float* C, int n) {
	int i;
	#pragma omp parallel for
	for (i = 0; i < n; i += BLOCKSIZE)
		gemm2(A,B,C,i,n);
}


int main(int argc, char* argv[]) {
	int n = 4000;
	float* A;
	float* B;
	float* C;
	int i, j, k;

	A = (float*) _mm_malloc(n*n*sizeof(TYPE), 16);
	B = (float*) _mm_malloc(n*n*sizeof(TYPE), 16);
	C = (float*) _mm_malloc(n*n*sizeof(TYPE), 16);

	float x = 0.0f;
	for (i = 0; i < n*n; i++) {
		C[i] = 0.0f;
		if (i % 2 == 0)
			x = 0.0f;
		else
			x = 1.0f;
		A[i] = x;
		B[i] = x;
	}

	gemm1(A,B,C,n);

	printf("%f %f\n", C[n*n-2], C[n*n-1]);
}
