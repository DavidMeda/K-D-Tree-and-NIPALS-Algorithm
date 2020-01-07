// Compilare:
//		gcc -O0 -m32 gemm0.c -o gemm0
//
// Test di velocità:
//		time ./gemm0
//


#include <stdlib.h>
#include <stdio.h>


#define TYPE float


void gemm0(TYPE** A, TYPE** B, TYPE** C, int n) {
	int i, j, k;
	for (i = 0; i < n; i ++)
		for (j = 0; j < n; j++)
			for (k = 0; k < n; k++)
				C[i][j] = C[i][j] + A[i][k] * B[k][j];
}


int main(int argc, char* argv[]) {
	const int n = 4000;
	TYPE** A;
	TYPE** B;
	TYPE** C;
	int i, j, k;
	TYPE x;
	
	A = (TYPE**) malloc(n*sizeof(TYPE*));
	B = (TYPE**) malloc(n*sizeof(TYPE*));
	C = (TYPE**) malloc(n*sizeof(TYPE*));
	
	x = 0;
	for (i = 0; i < n; i++) {
		A[i] = malloc(n*sizeof(TYPE));
		B[i] = malloc(n*sizeof(TYPE));
		C[i] = malloc(n*sizeof(TYPE));
		for (j = 0; j < n; j++) {
			C[i][j] = 0;
			A[i][j] = x;
			B[i][j] = x;
			x = 1 - x;
		}
	}
	
	gemm0(A, B, C, n);

	printf("%16.16f\n", C[n-1][n-1]);
}