; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;


;%include "sseutils.nasm"

A		equ		8
B		equ		12
C		equ		16
starti		equ		20
startj		equ		24
startk		equ		28
n		equ		32

dim		equ		4
p		equ		4
UNROLL		equ		4
BLOCKSIZE	equ		32

section .data			; Sezione contenente dati inizializzati

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	gemm5

gemm5:		push		ebp
		mov		ebp, esp
		push		ebx
		push		esi
		push		edi
		;
		;
		;
		mov		eax, [ebp+starti]	; i = starti
		imul		eax, dim
		
.fori:		mov		ebx, [ebp+startj]	; j = startj
		imul		ebx, dim
.forj:		
		mov		edi, [ebp+n]		; edi = n
		imul		edi, ebx		; edi = 4*j*n

		mov		esi, [ebp+C]		; esi = C
		add		esi, edi		; esi = C + 4*j*n
		movaps		xmm4, [eax+esi]		; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
		movaps		xmm5, [eax+esi+p*dim]
		movaps		xmm6, [eax+esi+2*p*dim]
		movaps		xmm7, [eax+esi+3*p*dim]
		
		mov		ecx, [ebp+startk]	; k = 0
		imul		ecx, dim

		mov		edx, ecx
		add		edx, dim*BLOCKSIZE
		
.fork:		
		mov		esi, [ebp+B]		; esi = B
		add		esi, edi		; esi = B + 4*j*n
		movss		xmm1, [ecx+esi]		; B[k][j] <--- Mem[B + 4*k + 4*j*n]
		shufps		xmm1, xmm1, 0

		mov		esi, [ebp+n]		; esi = n
		imul		esi, ecx		; esi = 4*k*n
		add		esi, [ebp+A]		; esi = A + 4*k*n

		movaps		xmm0, [eax+esi]		; A[i][k] <--- Mem[A + 4*i + 4*k*n]
		mulps		xmm0, xmm1
		addps		xmm4, xmm0		; 
		movaps		xmm0, [eax+esi+p*dim]	; A[i+1][k]
		mulps		xmm0, xmm1
		addps		xmm5, xmm0		;
		movaps		xmm0, [eax+esi+2*p*dim]	; A[i+2][k] 
		mulps		xmm0, xmm1
		addps		xmm6, xmm0		; 
		movaps		xmm0, [eax+esi+3*p*dim]	; A[i+3][k]
		mulps		xmm0, xmm1
		addps		xmm7, xmm0		; 

		add		ecx, dim
		cmp		ecx, edx
		jb		.fork

		mov		esi, [ebp+C]		; esi = C
		add		esi, edi		; esi = C + 4*j*n
		movaps		[eax+esi], xmm4		; C[i][k] <--- Mem[C + 4*i + 4*k*n]
		movaps		[eax+esi+p*dim], xmm5	;
		movaps		[eax+esi+2*p*dim], xmm6	;
		movaps		[eax+esi+3*p*dim], xmm7	;
		
		mov		edx, [ebp+startj]
		imul		edx, dim
		add		edx, dim*BLOCKSIZE
		add		ebx, dim
		cmp		ebx, edx
		jb		.forj
		
		mov		edx, [ebp+starti]
		imul		edx, dim
		add		edx, dim*BLOCKSIZE
		add		eax, p*dim*UNROLL	; 
		cmp		eax, edx
		jb		.fori
		;
		;
		;
		pop		edi
		pop		esi
		pop		ebx
		mov		esp, ebp
		pop		ebp
		ret

