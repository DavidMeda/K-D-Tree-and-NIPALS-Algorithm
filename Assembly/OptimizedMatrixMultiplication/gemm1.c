// Compilare:
//		gcc -O0 -m32 gemm1.c -o gemm1
//
// Test di velocità:
//		time ./gemm1
//


#include <stdlib.h>
#include <stdio.h>


#define TYPE float


void gemm1(TYPE* A, TYPE* B, TYPE* C, int n) {
	int i, j, k;
	for (i = 0; i < n; i ++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				C[i+j*n] = C[i+j*n] + A[i+k*n] * B[k+j*n];
}


int main(int argc, char* argv[]) {
	const int n = 4000;
	TYPE* A;
	TYPE* B;
	TYPE* C;
	int i, j, k;
	TYPE x;
	
	A = (TYPE*) calloc(n*n, sizeof(TYPE));
	B = (TYPE*) calloc(n*n, sizeof(TYPE));
	C = (TYPE*) calloc(n*n, sizeof(TYPE));
	
	x = 0;
	for (i = 0; i < n*n; i++) {
		C[i] = 0;
		A[i] = x;
		B[i] = x;
		x = 1 - x;
	}

	gemm1(A, B, C, n);

	printf("%16.16f\n", C[n*n-1]);
}