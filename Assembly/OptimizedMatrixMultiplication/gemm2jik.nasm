; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;

;
; Per eseguire:
;
;     ./run gemm2
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

n		equ		4000
dim		equ		4
p		equ		4

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

		mov		ebx, 0			; j = 0
forj:		imul		edi, ebx, n		; 4*j*n
		mov		eax, 0			; i = 0
fori:		movaps		xmm0, [C+eax+edi]	; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
		mov		ecx, 0			; k = 0
fork:		imul		esi, ecx, n		; 4*k*n
		movaps		xmm1, [A+eax+esi]	; A[i..i+p-1][k] = A[4*i+4*k*n..4*(i+*k*n+p-1)]
		movss		xmm2, [B+ecx+edi]	; B[k][j] = B[4*k+4*j*n]
		shufps		xmm2, xmm2, 0
		mulps		xmm2, xmm1
		addps		xmm0, xmm2		; c0 += A[i..i+p-1][k]*B[k+j*n]

		add		ecx, dim		; k++
		cmp		ecx, dim*n		; (k < n) ?
		jb		fork
		
		movaps		[C+eax+edi], xmm0	; C[i..i+p-1][j]  = C[i+j*n..i+j*n+p-1] = c0
		
		add		eax, dim*p		; i += p
		cmp		eax, dim*n		; (i < n) ?
		jb		fori
		
		add		ebx, dim		; j++
		cmp		ebx, dim*n		; (j < n) ?
		jb		forj

		printps		C+n*n-16, 1
				
		stop
