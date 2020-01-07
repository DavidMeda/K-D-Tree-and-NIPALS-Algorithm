; --------------------------------------------------
; Esempio di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per eseguire:
;
;     nasm -f elf32 primo.nasm 
;     gcc -m32 primo.o -o primo
;     ./primo
;
; oppure:
;
;     ./run primo
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

svec1:	db		'vec1 = ',0
svec2:	db		'vec2 = ',0
svec3:	db		'vec3 = ',0

align 16 
vec1:	dd		1.0, 2.0, 3.0, 4.0
align 16
vec2:	dd		1.0, 3.0, 5.0, 7.0

section .bss			; Sezione contenente dati non inizializzati

alignb 16
vec3:	resd	4

section .text			; Sezione contenente il codice macchina

global	main

main:	start

		movaps		xmm0, [vec1]
		movaps		xmm1, [vec2]
		addps		xmm0, xmm1
		movaps		[vec3], xmm0

		prints		svec1				; prints <string-address>
		printps		vec1, 1				; printps <ps-address>, <ps-length>
		prints		svec2
		printps		vec2, 1
		prints		svec3
		printps		vec3, 1

		stop
