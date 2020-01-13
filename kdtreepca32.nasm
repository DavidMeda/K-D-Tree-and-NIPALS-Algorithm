; ---------------------------------------------------------
; PQNN con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
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
; Per generare file oggetto:
;
;     nasm -f elf32 pqnn32.nasm 
;
%include "sseutils.nasm"



section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
;vec2:		resq	4

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
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

;   multi3: moltiplica matrice ds per il vettore colonna di V 
;   il risultato viene scritto in una colle della matrice U

global multi3

dataset      equ     8
V       equ     12
U       equ     16
cut		equ		20              ; parametro colonna fissato per colonna di V e U
riga	equ		24              ; parametro riga i che scorre per le righe di ds
k		equ		28
h		equ		32

dim		equ		4
UNROLL  equ     4
multi3:

        ; ------------------------------------------------------------
        ; Sequenza di ingresso nella funzione
        ; ------------------------------------------------------------
       
        push	ebp                 ; salva il Base Pointer
        mov		ebp, esp            ; il Base Pointer punta al Record di Attivazione corrente
        push	ebx                 ; salva i registri da preservare
        push	esi
        push	edi
        
        ; ------------------------------------------------------------
        ; legge i parametri dal Record di Attivazione corrente
        ; ------------------------------------------------------------
        
        mov     eax, [ebp+k]        ; eax = k
        mov		ecx, [ebp+riga]	    ; i = var iterativa da 0 a n
        imul    ecx, eax            ; ecx = i*k
        sub     eax, dim*UNROLL     ; eax= k-16
        imul    ecx, dim            ; ecx= i*k*4
        add     ecx, [ebp+dataset]  ; ecx= ds + i*k*4
        mov     edx, [ebp+cut]      ; edx= cut

        xor     esi, esi            ; esi è j variabile iterativa da 0 a k-4
        xorps   xmm3, xmm3          ;azzero xmm3 per le somme parziali


        ;non modificare ecx, esi, eax, edx
        ;ecx= ds + i*k*4
        ;esi = j
        ;eax= k-4
        ;edx= cut

    loop_q:
        cmp     esi, eax            ; if(j>k-16)
        jg      h_add               ; se eax =k-16 salta alle istruzioni h_add altrimenti continua

        ; devo riempire il registro xmm0 con 4 valori non consecutivi di V
        movups  xmm0, [ecx+4*esi]   ; [ds + 4*i*k + j*4]
        movups  xmm4, [ecx+4*esi+1*dim*UNROLL]   ; [ds + 4*i*k + j*4 + 16]
        movups  xmm5, [ecx+4*esi+2*dim*UNROLL]   ; [ds + 4*i*k + j*4 + 32]
        movups  xmm6, [ecx+4*esi+3*dim*UNROLL]   ; [ds + 4*i*k + j*4 + 48]

        ;1° iterazione su V
        mov     edi, 4              ; edi: cont è il contatore fino a 4
        xorps   xmm2,xmm2           ; azzero xmm2

        shf1:
            shufps	xmm2, xmm2, 57  ; shift di una posizione a sinistra

        loop_4:
            mov     ebx, [ebp+h]    ; ebx = h
            imul    ebx, dim        ; ebx= h*4
            imul    ebx, esi        ; ebx= h*4*j
            add     ebx, [ebp+V]    ; ebx =V+ h*4*j
            
            ;prendo un valore di V
            movss   xmm1, [ebx+edx*4]; [V + 4*j*h + 4*cut] 
            xorps   xmm2, xmm1      ; sposto il valore in xmm2
            
            dec     edi             ; cont--
            inc     esi             ; j++ (0 < j < k-4)
            cmp     edi, 0          ; devo fare 4 iterazioni
            jg      shf1      
            jne     loop_4

        shufps	xmm2, xmm2, 57      ; shift di una posizione a destra
        mulps   xmm0, xmm2          ; xmm0 moltiplico i 4 valori di ds con quelli di V
        addps   xmm3, xmm0          ; xmm3 registro per somma parziale dei valori moltiplicati

        ;2° iterazione di V
        mov     edi, 4              ; edi: cont è il contatore fino a 4
        xorps   xmm2,xmm2           ; azzero xmm2

        shf2:
            shufps	xmm2, xmm2, 57  ; shift di una posizione a sinistra

        loop_4_2:
            mov     ebx, [ebp+h]    ; ebx = h
            imul    ebx, dim        ; ebx= h*4
            imul    ebx, esi        ; ebx= h*4*j
            add     ebx, [ebp+V]    ; ebx =V+ h*4*j
            
            ;prendo un valore di V
            movss   xmm1, [ebx+edx*4+1*dim*UNROLL]; [V + 4*j*h + 4*cut] 
            xorps   xmm2, xmm1      ; sposto il valore in xmm2
            
            dec     edi             ; cont--
            inc     esi             ; j++ (0 < j < k-4)
            cmp     edi, 0          ; devo fare 4 iterazioni
            jg      shf2       
            jne     loop_4_2

        shufps	xmm2, xmm2, 57      ; shift di una posizione a destra
        mulps   xmm4, xmm2          ; xmm0 moltiplico i 4 valori di ds con quelli di V
        addps   xmm3, xmm4          ; xmm3 registro per somma parziale dei valori moltiplicati

        ;3° iterazione V
        mov     edi, 4              ; edi: cont è il contatore fino a 4
        xorps   xmm2,xmm2           ; azzero xmm2

        shf3:
            shufps	xmm2, xmm2, 57  ; shift di una posizione a sinistra

        loop_4_3:
            mov     ebx, [ebp+h]    ; ebx = h
            imul    ebx, dim        ; ebx= h*4
            imul    ebx, esi        ; ebx= h*4*j
            add     ebx, [ebp+V]    ; ebx =V+ h*4*j
            
            ;prendo un valore di V
            movss   xmm1, [ebx+edx*4+3*dim*UNROLL]; [V + 4*j*h + 4*cut] 
            xorps   xmm2, xmm1      ; sposto il valore in xmm2
            
            dec     edi             ; cont--
            inc     esi             ; j++ (0 < j < k-4)
            cmp     edi, 0          ; devo fare 4 iterazioni
            jg      shf3       
            jne     loop_4_3

        shufps	xmm2, xmm2, 57      ; shift di una posizione a destra
        mulps   xmm5, xmm2          ; xmm0 moltiplico i 4 valori di ds con quelli di V
        addps   xmm3, xmm5          ; xmm3 registro per somma parziale dei valori moltiplicati

        ; 4° iterazione di V
        mov     edi, 4              ; edi: cont è il contatore fino a 4
        xorps   xmm2,xmm2           ; azzero xmm2

        shf4:
            shufps	xmm2, xmm2, 57  ; shift di una posizione a sinistra

        loop_4_4:
            mov     ebx, [ebp+h]    ; ebx = h
            imul    ebx, dim        ; ebx= h*4
            imul    ebx, esi        ; ebx= h*4*j
            add     ebx, [ebp+V]    ; ebx =V+ h*4*j
            
            ;prendo un valore di V
            movss   xmm1, [ebx+edx*4+4*dim*UNROLL]; [V + 4*j*h + 4*cut] 
            xorps   xmm2, xmm1      ; sposto il valore in xmm2
            
            dec     edi             ; cont--
            inc     esi             ; j++ (0 < j < k-4)
            cmp     edi, 0          ; devo fare 4 iterazioni
            jg      shf4       
            jne     loop_4_4

        shufps	xmm2, xmm2, 57      ; shift di una posizione a destra
        mulps   xmm6, xmm2          ; xmm0 moltiplico i 4 valori di ds con quelli di V
        addps   xmm3, xmm6          ; xmm3 registro per somma parziale dei valori moltiplicati
        jmp     loop_q    

    h_add:       
        haddps xmm3, xmm3           ; riduco la somma a un valore solo
        haddps xmm3, xmm3           ; riduco la somma a un valore solo
   
   loop_r:
        mov     eax, [ebp+k]        ; eax = k
        cmp     esi, eax
        jge     end                 ; se j == k no loop resto vai a end altrimenti loop_r
        
        ; prendo il valore di ds
        movss   xmm0, [ecx+4*esi]   ; [ds + 4*i*k + j*4]
        
        ;prendo il valore di V
        mov     ebx, [ebp+h]        ; ebx = h
        imul    ebx, dim            ; ebx= h*4
        imul    ebx, esi            ; ebx= h*4*j
        add     ebx, [ebp+V]        ; ebx = h*4*j+V
        ; movss   xmm1, [ebx+edx*4]   ; [V + 4*j*h + 4*cut]
        mulss   xmm0,[ebx+edx*4]    ; moltiplico il valore di ds per quello di V [V + 4*j*h + 4*cut]
        addps   xmm3, xmm0          ; aggiungo alle somme parziali il risultato ottenuto

        inc esi
        jmp loop_r

    end:
        mov		ecx, [ebp+riga]	    ; ecx è i = var iterativa da 0 a n
        imul    ecx, [ebp+h]        ; ecx = i*h
        imul    ecx, dim            ; ecx = i*h*4
        add     ecx, [ebp+U]        ; ecx = i*h*4 + U
        movss   [ecx+ edx*4], xmm3  ; aggiungo a [U+ i*4*h + cut*4] il valore in xmm3
        
        ; ------------------------------------------------------------
        ; Sequenza di uscita dalla funzione
        ; ------------------------------------------------------------
        pop	edi                     ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                ; ripristina lo Stack Pointer
        pop	ebp                     ; ripristina il Base Pointer
        ret                         ; torna alla funzione C chiamante

global euc_dist


dataset      equ     8
i1       equ     12
queryset       equ     16
i2		equ		20
k		equ		24
output  equ  28
dim		equ		4
UNROLL		equ		16

section .data                             ; Sezione contenente dati inizializzati
section .bss                              ; Sezione contenente dati non inizializzati
section .text                             ; Sezione contenente il codice macchina

euc_dist:
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
        xorps xmm2,xmm2
        xorps xmm1,xmm1
        xorps xmm0,xmm0
        xor     esi, esi  ;esi è j variabile iterativa da 0 a k-4
    loop_q:
        mov		ecx, [ebp+i1]	; i1 = var iterativa su ds da 0 a n
        mov		edx, [ebp+i2]	; i2 = var iterativa su qs da 0 a n
        mov     eax, [ebp+k]    ;eax = k
        sub     eax, UNROLL       ;eax= k-unroll
        cmp     esi, eax            ;if(j>k-4)
        jg      h_add               ;se eax =k-4 salta alle istruzioni h_add altrimenti continua
        mov     edi, [ebp+k]        ;edi= k
        imul    ecx, edi            ;ecx = i1*k
        imul    ecx, dim            ;ecx= i1*k*4
        imul    edx, edi            ;ecx = i2*k
        imul    edx, dim            ;edx= i2*k*4
        mov     eax, [ebp+dataset]     ;eax= ds
        add     ecx, eax            ;ecx= ds + i1*k*4
        mov     ebx, [ebp+queryset]     ;ebx= qs
        add     edx, ebx            ;edx= qs + i2*k*4
        ;non modificcare ecx e edx ed esi
        movups  xmm0, [ecx+4*esi]  ;[ds + 4*i1*k + j*4]
        movups  xmm1, [edx+4*esi]  ;[qs + 4*i2*k + j*4]
        subps   xmm0, xmm1          ;xmm0 registro per somma parziale dei valori sotratti
        mulps   xmm0, xmm0         ;xmm0 per il quadrato
        addps   xmm2, xmm0          ;xmm3 registro per somma parziale dei valori
        movups  xmm0, [ecx+4*esi+4*dim]  ;[ds + 4*i1*k + j*4]
        movups  xmm1, [edx+4*esi+4*dim]  ;[qs + 4*i2*k + j*4]
        subps   xmm0, xmm1          ;xmm0 registro per somma parziale dei valori sotratti
        mulps   xmm0, xmm0         ;xmm0 per il quadrato
        addps   xmm2, xmm0          ;xmm3 registro per somma parziale dei valori
        movups  xmm0, [ecx+4*esi+8*dim]  ;[ds + 4*i1*k + j*4]
        movups  xmm1, [edx+4*esi+8*dim]  ;[qs + 4*i2*k + j*4]
        subps   xmm0, xmm1          ;xmm0 registro per somma parziale dei valori sotratti
        mulps   xmm0, xmm0         ;xmm0 per il quadrato
        addps   xmm2, xmm0          ;xmm3 registro per somma parziale dei valori
        movups  xmm0, [ecx+4*esi+12*dim]  ;[ds + 4*i1*k + j*4]
        movups  xmm1, [edx+4*esi+12*dim]  ;[qs + 4*i2*k + j*4]
        subps   xmm0, xmm1          ;xmm0 registro per somma parziale dei valori sotratti
        mulps   xmm0, xmm0         ;xmm0 per il quadrato
        addps   xmm2, xmm0          ;xmm3 registro per somma parziale dei valori
        add esi,UNROLL
        jmp     loop_q    

    h_add:       
        haddps xmm2, xmm2      ;riduco la somma a un valore solo
        haddps xmm2, xmm2       ;riduco la somma a un valore solo
     loop_r:
        mov     eax, [ebp+k]    ;eax = k
        cmp     esi, eax
        jge      end       ;se j == k no loop resto vai a end altrimenti loop_r
        mov		ecx, [ebp+i1]	; i1 = var iterativa su ds da 0 a n
        imul    ecx, eax            ;ecx = i1*k
        imul    ecx, dim            ;ecx= i1*k*4
        mov     eax, [ebp+dataset]     ;eax= ds
        add     ecx, eax            ;ecx= ds + i1*k*4

        mov		edx, [ebp+i2]	; i2 = var iterativa su qs da 0 a n
        mov     eax, [ebp+k]    ;eax = k
        imul    edx, eax            ;edx = i2*k
        imul    edx, dim            ;edx= i2*k*4
        mov     ebx, [ebp+queryset]     ;ebx= qs
        add     edx, ebx            ;edx= ds + i1*k*4
        movss   xmm0, [ecx+4*esi] ; [ds + 4*i1*k + j*4]
        movss   xmm1, [edx+4*esi] ; [qs + 4*i2*k + j*4]
        subss   xmm0, xmm1          ;xmm0 registro per sottrazione parziale dei valori sotratti
        mulss   xmm0, xmm0         ;xmm0 per il quadrato
        addss   xmm2, xmm0          ;xmm3 registro per somma parziale dei valori   
        inc esi
        jmp loop_r

    end:
        mov eax,[ebp+output]
        movss     [eax],xmm2      ; eax = risultato
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante

