; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;

;
; Per eseguire:
;
;     ./run gemm3
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

n		equ		4000
dim		equ		4
p		equ		4
UNROLL		equ		4

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
.forj:		imul		edi, ebx, n		; 4*j*n

		movaps		xmm4, [C+eax+edi]	; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
		movaps		xmm5, [C+eax+edi+16]	; + p*dim
		movaps		xmm6, [C+eax+edi+32]	; + 2*p*dim
		movaps		xmm7, [C+eax+edi+48]	; + 3*p*dim
		
		mov		ecx, 0			; k = 0
.fork:		imul		esi, ecx, n		; 4*k*n

		movss		xmm1, [B+ecx+edi]	; B[k][j] = B[4*k+4*j*n]
		shufps		xmm1, xmm1, 0

		movaps		xmm0, [A+eax+esi]	; 
		mulps		xmm0, xmm1
		addps		xmm4, xmm0		; 
		movaps		xmm0, [A+eax+esi+16]	; 
		mulps		xmm0, xmm1
		addps		xmm5, xmm0		; 
		movaps		xmm0, [A+eax+esi+32]	; 
		mulps		xmm0, xmm1
		addps		xmm6, xmm0		; 
		movaps		xmm0, [A+eax+esi+48]	; 
		mulps		xmm0, xmm1
		addps		xmm7, xmm0		; 

		add		ecx, dim
		cmp		ecx, dim*n
		jb		.fork

		movaps		[C+eax+edi], xmm4	;
		movaps		[C+eax+edi+16], xmm5	;
		movaps		[C+eax+edi+32], xmm6	;
		movaps		[C+eax+edi+48], xmm7	;
		
		add		ebx, dim
		cmp		ebx, dim*n
		jb		.forj
		
		add		eax, p*dim*UNROLL	; UNROLL = 4
		cmp		eax, dim*n
		jb		.fori

		printps		C+n*n-16, 1

		stop
		
