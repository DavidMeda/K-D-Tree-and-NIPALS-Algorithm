; %include "sseutils.nasm"

A		equ		8
n       equ     12        ; variabile iterativa
size     equ     16

dim		equ		4
p		equ		4
UNROLL		equ		4
BLOCKSIZE	equ		32

section .data             ; Sezione contenente dati inizializzati

section .bss              ; Sezione contenente dati non inizializzati

section .text             ; Sezione contenente il codice macchina

global prova

prova:
        push		ebp
        mov		ebp, esp
        push		ebx
        push		esi
        push		edi


        mov		eax, [ebp+n] ; i = n varibile itarativa
        imul		eax, dim

    fori:
        

