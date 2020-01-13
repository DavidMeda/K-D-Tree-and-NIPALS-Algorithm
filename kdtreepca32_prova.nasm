; ---------------------------------------------------------
; PQNN con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017



; Software necessario per l'esecuzione:

; NASM (www.nasm.us)
; GCC (gcc.gnu.org)

; entrambi sono disponibili come pacchetti software
; installabili mediante il packaging tool del sistema
; operativo; per esempio, su Ubuntu, mediante i comandi:

; sudo apt-get install nasm
; sudo apt-get install gcc

; potrebbe essere necessario installare le seguenti librerie:

; sudo apt-get install lib32gcc-4.8-dev (o altra versione)
; sudo apt-get install libc6-dev-i386

; Per generare file oggetto:

; nasm -f elf32 pqnn32.nasm

%include "sseutils.nasm"

section .data                             ; Sezione contenente dati inizializzati

uno:		dd		1.0

; align 16
; inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss                              ; Sezione contenente dati non inizializzati

; alignb 16
; vec2:		resq	4

section .text                             ; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria

; getmem	<size>,<elements>

; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

; fremem	<address>

; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
    mov	eax, %1
    push	eax
    mov	eax, %2
    push	eax
    call	get_block
    add	esp, 8
%endmacro

%macro	fremem	1
    push	%1
    call	free_block
    add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global sommatoria

data		equ		8
n 		equ 	12
k 		equ 	16
col 	equ 	20

; msg	db	'n:',0
; nl	db	10,0

sommatoria:
        ; ------------------------------------------------------------
        ; Sequenza di ingresso nella funzione
        ; ------------------------------------------------------------
        push		ebp                         ; salva il Base Pointer
        mov			ebp, esp                    ; il Base Pointer punta al Record di Attivazione corrente
        push		ebx                         ; salva i registri da preservare
        push		esi
        push		edi
        ; ------------------------------------------------------------
        ; legge i parametri dal Record di Attivazione corrente
        ; ------------------------------------------------------------

        ; [EAX] contiene l'indirizzo di memoria del dataset
        ; [EAX+4] contiene il numero di righe del dataset (n)
        ; [EAX+8] contiene il numero di colonne del dataset (k)
        ; [EAX+12] contiene il numero della colonna di iterazione (cut)
        mov     eax,[ebp+data]      ;eax= primo valore dataseet
        mov 	esi, 0              ; esi = 0
        ; mov 	ebx, [ebp+n]        ; in ebx=n
        mov 	edx, [ebp+k]        ; in edx ci sarà il valore k
        imul	edx, 4              ;edx = 4*k
        mov		edi, [ebp+col]      ; in edi ci sarà la varibile di cut
        imul    edi, 4              ;edi = 4*cut

        mov     ebx,4           ;cont for4
    .for4:
        mov     ecx, 1
        imul    ecx, edx        ;ecx= i
        imul    ecx, esi        ;ecx= i*4*k
        add     ecx, edi        ;ecx= i*4*k+4*cut
        movss   xmm0,[eax+ecx]  ;[base + i*4*k + 4*cut]
        ; shufps	xmm0, xmm0, 147 ;
        xorps   xmm1, xmm0
        shufps	xmm1, xmm1, 147
        add     esi,1           ;i++
        sub     ebx,1           ; cont--

        cmp     ebx, 0
        jnz     .for4
        
        
        movups 	[eax],	xmm1


        ; ------------------------------------------------------------
        ; Sequenza di uscita dalla funzione
        ; ------------------------------------------------------------

        pop	edi                           ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                      ; ripristina lo Stack Pointer
        pop	ebp                           ; ripristina il Base Pointer
        ret                               ; torna alla funzione C chiamante
