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

		mov		eax, 0			; i = 0
fori:		mov		ebx, 0			; j = 0
forj:		imul		edi, ebx, n		; 4*j*n
		movss		xmm0, [C+eax+edi]	; c0 = C[i][j]		
		mov		ecx, 0			; k = 0
fork:		imul		esi, ecx, n		; 4*k*n
		movss		xmm1, [A+eax+esi]	; xmm1 = A[i][k]
		movss		xmm2, [B+ecx+edi]	; xmm2 = B[k][j]
		mulss		xmm2, xmm1
		addss		xmm0, xmm2		; c0 += A[i][k]*B[k][j]

		add		ecx, dim		; k++
		cmp		ecx, dim*n		; (k < n) ?
		jb		fork
		
		movss		[C+eax+edi], xmm0	; C[i][j] = c0
		
		add		ebx, dim		; j++
		cmp		ebx, dim*n		; (j < n) ?
		jb		forj
		
		add		eax, dim		; i ++
		cmp		eax, dim*n		; (i < n) ?
		jb		fori

		printps		C+n*n-16, 1
				
		stop
