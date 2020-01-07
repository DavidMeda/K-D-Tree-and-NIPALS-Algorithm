; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;

;
; Per eseguire:
;
;     ./run gemm4
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

n		equ		4000
dim		equ		4
p		equ		4
UNROLL		equ		4
BLOCKSIZE	equ		32

; 
; Le matrici si assumono memorizzate per colonne
;

align 16
inizio:		dd		0.0, 1.0, 0.0, 1.0
align 16
quattro:	dd		0.0, 0.0, 0.0, 0.0

section .bss			; Sezione contenente dati non inizializzati

alignb 16
A:		resd	n*n
alignb 16
B:		resd	n*n
alignb 16
C:		resd	n*n

section .text			; Sezione contenente il codice macchina

global	main

main:	start

		; ----------------------------------------
		; carica le matrici A e B e azzera la matrice C
		;
		movaps		xmm0, [inizio]
		movaps		xmm1, [quattro]
		xorps		xmm2, xmm2
		mov		ebx, 0
		mov		ecx, n*n/4
ciclof:		movaps		[A+ebx], xmm0
		movaps		[B+ebx], xmm0
		movaps		[C+ebx], xmm2
		addps		xmm0, xmm1
		add		ebx, 16
		dec		ecx
		jnz		ciclof

		; ----------------------------------------
		; Prodotto di matrici single-precision

		mov		eax, 0			; i = 0
.fori:		mov		ebx, 0			; j = 0
.forj:		mov		ecx, 0			; k = 0
.fork:		
		push		eax
		push		ebx
		push		ecx
		call		gemm_block
		pop		ecx
		pop		ebx
		pop		eax
		
		add		ecx, dim*BLOCKSIZE
		cmp		ecx, dim*n
		jb		.fork

		add		ebx, dim*BLOCKSIZE
		cmp		ebx, dim*n
		jb		.forj
		
		add		eax, dim*BLOCKSIZE
		cmp		eax, dim*n
		jb		.fori
		
		printps		C+n*n-16, 1
		
		stop
		
		

starti		equ		16
startj		equ		12
startk		equ		8


gemm_block:	push		ebp
		mov		ebp, esp

		mov		eax, [ebp+starti]	; i = starti
.fori:		mov		ebx, [ebp+startj]	; j = startj
.forj:		imul		edi, ebx, n		; 4*j*n

		;prefetchnta	[C+eax+edi+4*n]
		;prefetchnta	[C+eax+edi+4*n+dim*p]
		;prefetchnta	[C+eax+edi+4*n+2*dim*p]
		;prefetchnta	[C+eax+edi+4*n+3*dim*p]

		movaps		xmm4, [C+eax+edi]	; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
		movaps		xmm5, [C+eax+edi+p*dim]
		movaps		xmm6, [C+eax+edi+2*p*dim]
		movaps		xmm7, [C+eax+edi+3*p*dim]
		
		mov		ecx, [ebp+startk]	; k = 0

		mov		edx, ecx
		add		edx, dim*BLOCKSIZE
		
.fork:		imul		esi, ecx, n		; 4*k*n

		;prefetcht0	[A+eax+esi+4*n]

		;prefetcht0	[B+ecx+edi+dim*n]
		
		movss		xmm1, [B+ecx+edi]	; B[k][j] = B[4*k+4*j*n]
		shufps		xmm1, xmm1, 0

		movaps		xmm0, [A+eax+esi]	; A[i][k]
		mulps		xmm0, xmm1
		addps		xmm4, xmm0		; 
		movaps		xmm0, [A+eax+esi+p*dim]	; A[i+1][k]
		mulps		xmm0, xmm1
		addps		xmm5, xmm0		;
		movaps		xmm0, [A+eax+esi+2*p*dim]	; A[i+2][k] 
		mulps		xmm0, xmm1
		addps		xmm6, xmm0		; 
		movaps		xmm0, [A+eax+esi+3*p*dim]	; A[i+3][k]
		mulps		xmm0, xmm1
		addps		xmm7, xmm0		; 

;		mov		edx, [ebp+startk]
;		add		edx, dim*BLOCKSIZE
		add		ecx, dim
		cmp		ecx, edx
		jb		.fork

		movaps		[C+eax+edi], xmm4	;
		movaps		[C+eax+edi+p*dim], xmm5	;
		movaps		[C+eax+edi+2*p*dim], xmm6	;
		movaps		[C+eax+edi+3*p*dim], xmm7	;
		
		mov		edx, [ebp+startj]
		add		edx, dim*BLOCKSIZE
		add		ebx, dim
		cmp		ebx, edx
		jb		.forj
		
		mov		edx, [ebp+starti]
		add		edx, dim*BLOCKSIZE
		add		eax, p*UNROLL*dim		; UNROLL = 4
		cmp		eax, edx
		jb		.fori

		pop		ebp
		ret

