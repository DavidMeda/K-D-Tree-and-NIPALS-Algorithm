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

global calcolaSqrt

valoreSqrt  equ     8

calcolaSqrt:
    push	ebp                 ; salva il Base Pointer
    mov		ebp, esp            ; il Base Pointer punta al Record di Attivazione corrente
    push	ebx                 ; salva i registri da preservare
    push	esi
    push	edi


    mov     eax, [ebp+valoreSqrt]
    movss   xmm0,[eax]
    sqrtss     xmm0, xmm0
    movss       [eax], xmm0

    pop	edi                     ; ripristina i registri da preservare
    pop	esi
    pop	ebx
    mov	esp, ebp                ; ripristina lo Stack Pointer
    pop	ebp                     ; ripristina il Base Pointer
    ret                         ; torna alla funzione C chiamante




global dividiAss

vett        equ     8
numEleVett  equ     12
value       equ     16

dividiAss:

    push	ebp                 ; salva il Base Pointer
    mov		ebp, esp            ; il Base Pointer punta al Record di Attivazione corrente
    push	ebx                 ; salva i registri da preservare
    push	esi
    push	edi

    mov     eax, [ebp+vett]         ; eax = vett
    mov     edi, [ebp+numEleVett]   ; edi = numEle
    sub     edi, 16                 ; edi = numele - 16

    xor esi, esi
    movss   xmm0, [ebp+value]       ; valore divisore
    shufps  xmm0, xmm0, 0h
    movups  xmm5,xmm0
    movups  xmm6,xmm0
    movups  xmm7,xmm0


    loop16_div:

        cmp esi, edi
        jg  Bloop4_div
        movups  xmm1, [eax+ 4*esi]   ; [vettoer + 4*i]
        movups  xmm2, [eax+ 4*esi+16]   ; [vettoer + 4*i]
        movups  xmm3, [eax+ 4*esi+32]   ; [vettoer + 4*i]
        movups  xmm4, [eax+ 4*esi+48]   ; [vettoer + 4*i]
        divps   xmm1, xmm0
        divps   xmm2, xmm5
        divps   xmm3, xmm6
        divps   xmm4, xmm7
        movups  [eax+4*esi], xmm1
        movups  [eax+4*esi+16], xmm2
        movups  [eax+4*esi+32], xmm3
        movups  [eax+4*esi+48], xmm4
        add     esi, 16
        jmp     loop16_div
    
    Bloop4_div:
    add     edi, 12                      ;edi = numEle - 4

    loop4_div:
        cmp     esi, edi
        jg      Bloop1_div

        movups  xmm1, [eax+ 4*esi]   ; [vettoer + 4*i]
        divps   xmm1, xmm0
        movups  [eax+4*esi], xmm1

        add     esi, 4
        jmp loop4_div
    Bloop1_div:
    add     edi, 4

    loop1_div:
        cmp     esi, edi
        jge     fine_div
        movss  xmm1, [eax+ 4*esi]   ; [vettoer + 4*i]
        divss   xmm1, xmm0
        movss  [eax+4*esi], xmm1
        inc     esi
        jmp     loop1_div
    
    fine_div:
        pop	edi                     ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                ; ripristina lo Stack Pointer
        pop	ebp                     ; ripristina il Base Pointer
        ret                         ; torna alla funzione C chiamante



global aggiornaDatasetAss

datasetA    equ     8
vettoreU    equ     12 
vettoreV    equ     16 
rigaA       equ     20
kcoll        equ     24

aggiornaDatasetAss:

    push	ebp                 ; salva il Base Pointer
    mov		ebp, esp            ; il Base Pointer punta al Record di Attivazione corrente
    push	ebx                 ; salva i registri da preservare
    push	esi
    push	edi
    
    mov     edi, [ebp+rigaA]    ; edi = i
    mov     eax, [ebp+vettoreU] ; eax = u base
    mov     ebx, [ebp+vettoreV] ; eax = v base
    mov     ecx, [ebp+kcoll]     ; ecx = k size

    xor     esi, esi            ; esi = j

    movss   xmm0, [eax+4*edi] ; xmm0= [u + i*4]
    shufps  xmm0, xmm0, 0h      ;riproduco 1° valore in tutto xmm0

    mov     eax, [ebp+datasetA] ; eax = dataset base
    imul    edi, 4            ; edi = i *4
    imul     edi, ecx         ;edi = i*4*k
    add     eax, edi                ; eax= ds + i*4*k
    sub     ecx, 16     ; ecx = k-16

    loop16:

        cmp    esi, ecx
        jg     Bloop4_r

        ; 1° ite
        movups  xmm1, [ebx+4*esi]  ; xmm1 = [v + j*4]
        mulps   xmm1, xmm0
        movups  xmm2, [eax+esi*4]
        subps  xmm2, xmm1
        movups  [eax+esi*4],xmm2    ; store in ds

        ; 2° ite
        movups  xmm1, [ebx+4*esi+16]  ; xmm1 = [v + j*4+16]
        mulps   xmm1, xmm0
        movups  xmm2, [eax+esi*4+16]
        subps   xmm2, xmm1
        movups  [eax+esi*4+16],xmm2    ; store in ds

        ; 3° ite
        movups  xmm1, [ebx+4*esi+32]  ; xmm1 = [v + j*4]
        mulps   xmm1, xmm0
        movups  xmm2, [eax+4*esi+32]
        subps   xmm2, xmm1
        movups  [eax+4*esi+32],xmm2    ; store in ds

        ; 4° ite
        movups  xmm1, [ebx+4*esi+48]  ; xmm1 = [v + j*4]
        mulps   xmm1, xmm0
        movups  xmm2, [eax+4*esi+48]
        subps   xmm2, xmm1
        movups  [eax+4*esi+48],xmm2    ; store in ds

        add     esi, 16

        jmp     loop16
    Bloop4_r:

        add     ecx, 12            ; ecx = k-4

    loop4_r:
        cmp     esi, ecx
        jg     Bloop1_r

        movups  xmm1, [ebx+4*esi]  ; xmm1 = [v + j*4]
        mulps   xmm1, xmm0
        movups  xmm2, [eax+esi*4]
        subps   xmm2, xmm1
        movups  [eax+esi*4],xmm2    ; store in ds
       
        add     esi, 4
        jmp     loop4_r
    Bloop1_r:
        add     ecx, 4     ; ecx = k size

    loop1_r:
        cmp     esi, ecx
        jge     fine

        movss  xmm1, [ebx+4*esi]  ; xmm1 = [v + j*4]
        mulss   xmm1, xmm0
        movss   xmm2, [eax+esi*4]
        subss   xmm2, xmm1
        movss  [eax+esi*4],xmm2    ; [ds + k*i*4 + j*4]

        inc     esi
        jmp     loop1_r


    fine:
        pop	edi                     ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                ; ripristina lo Stack Pointer
        pop	ebp                     ; ripristina il Base Pointer
        ret                         ; torna alla funzione C chiamante



global calcolaTAss

vect      equ     8
numEle  equ     12
result  equ     16

calcolaTAss:

    push	ebp                 ; salva il Base Pointer
    mov		ebp, esp            ; il Base Pointer punta al Record di Attivazione corrente
    push	ebx                 ; salva i registri da preservare
    push	esi
    push	edi

    mov     eax, [ebp+vect]       ; eax = index base vettore
    mov     ebx, [ebp+numEle]       ; ebx = size vettore
    sub     ebx, 16         ; ebx = size - 16

    xor     esi, esi                ; variabile i di iterazione
    xorps   xmm0, xmm0              ; xmm0 somme parziali

    loop_quo:
        cmp     esi, ebx
        jg      Bloop_4

        ;prendo 16 valori consecutivi
        movups  xmm1, [eax+4*esi]                ; [vect + 4*i ]
        movups  xmm2, [eax+4*esi+16]   ; [vect + 4*i  + 16]
        movups  xmm3, [eax+4*esi+32]   ; [vect + 4*i  + 32]
        movups  xmm4, [eax+4*esi+48]   ; [vect + 4*i  + 48]
        mulps   xmm1, xmm1
        mulps   xmm2, xmm2
        mulps   xmm3, xmm3
        mulps   xmm4, xmm4
        addps   xmm0, xmm1
        addps   xmm0, xmm2
        addps   xmm0, xmm3
        addps   xmm0, xmm4

        add esi, 16
        jmp loop_quo

    Bloop_4:
        add     ebx, 12               ; ebx = size elem

    loop_4:
        cmp     esi, ebx
        jg     h_add_2

        movups  xmm1, [eax+4*esi]                ; [vect + 4*i + j*4]
        mulps   xmm1, xmm1
        addps   xmm0, xmm1

        add     esi, 4
        jmp     loop_4


    h_add_2:
        haddps  xmm0, xmm0                      ; riduzione somma
        haddps  xmm0, xmm0                      ; riduzione somma

        add     ebx, 4               ; ebx = size elem
    
    loop_rest:
        cmp     esi, ebx
        jge     end_2


        movss   xmm1, [eax+esi*4]               ;[vect + i*4] 
        mulss   xmm1, xmm1
        addss   xmm0, xmm1

        inc     esi
        jmp loop_rest

    end_2:
        mov     eax, [ebp+result]
        movss   [eax], xmm0

    pop	edi                     ; ripristina i registri da preservare
    pop	esi
    pop	ebx
    mov	esp, ebp                ; ripristina lo Stack Pointer
    pop	ebp                     ; ripristina il Base Pointer
    ret                         ; torna alla funzione C chiamante




;   prodottoMatriceAss: moltiplica matrice ds per il vettore colonna di V 
;   il risultato viene scritto in una colle della matrice U

global prodottoMatriceAss

dataset      equ     8
V       equ     12
U       equ     16
rigaI	equ		20              ; parametro riga i che scorre per le righe di ds
k		equ		24




prodottoMatriceAss:

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
        mov		edx, [ebp+rigaI]	    ; edx = n     i = var iterativa da 0 a n
        mov     eax, [ebp+k]        ; eax = k
        xor     edi, edi                

    forI:
        cmp     edi, edx
        je      close
        mov     ecx, edi
        imul    ecx, eax            ; ecx = i*k
        sub     eax, 16     ; eax= k-16
        imul    ecx, 4            ; ecx= i*k*4
        add     ecx, [ebp+dataset]  ; ecx= ds + i*k*4
        mov     ebx, [ebp+V]        ; ebx = V        

        xor     esi, esi            ; esi è j variabile iterativa da 0 a k-4
        xorps   xmm3, xmm3          ;azzero xmm3 per le somme parziali

        ;non modificare ecx, esi, eax, ebx
        ;ecx= ds + i*k*4
        ;esi = j
        ;eax= k-4
        ;ebx = V 

        loop_q_1:
            cmp     esi, eax            ; if(j>k-16)
            jg      h_add_1               ; se eax =k-16 salta alle istruzioni h_add altrimenti continua

            ;prendo 16 valori consecutivi di ds
            movups  xmm0, [ecx+4*esi]                ; [ds + 4*i*k + j*4]
            movups  xmm4, [ecx+4*esi+16]   ; [ds + 4*i*k + j*4 + 16]
            movups  xmm5, [ecx+4*esi+32]   ; [ds + 4*i*k + j*4 + 32]
            movups  xmm6, [ecx+4*esi+48]   ; [ds + 4*i*k + j*4 + 48]

            ;prendo 16 valori consecutivi di v

            mulps  xmm0, [ebx+esi*4]                ; [V+j*4]
            addps   xmm3, xmm0                      ; xmm3 somme parziali

            mulps   xmm4, [ebx+esi*4+16]    ; [V+j*4+ 16]
            addps   xmm3, xmm4                      ; xmm3 somme parziali

            mulps   xmm5, [ebx+esi*4+32]    ; [V+j*4+32]
            addps   xmm3, xmm5                      ; xmm3 somme parziali

            mulps   xmm6, [ebx+esi*4+48]    ; [V+j*4+48]
            addps   xmm3, xmm6                      ; xmm3 somme parziali


            add     esi, 16    ; esi += 16
            jmp     loop_q_1    

        h_add_1:       
            haddps xmm3, xmm3           ; riduco la somma a un valore solo
            haddps xmm3, xmm3           ; riduco la somma a un valore solo
            add    eax, 16       ; eax = k
    
    loop_r_1:
            cmp     esi, eax
            jge     end_1                 ; se j == k no loop resto vai a end altrimenti loop_r
            
            ; prendo il valore di ds
            movss   xmm0, [ecx+4*esi]   ; [ds + 4*i*k + j*4]
            
            ;prendo il valore di V
            mulss   xmm0,[ebx+esi*4]    ; moltiplico il valore di ds per quello di V [V + 4*j*h + 4*cut]
            addps   xmm3, xmm0          ; aggiungo alle somme parziali il risultato ottenuto

            inc esi
            jmp loop_r_1

        end_1:
            mov     esi, [ebp+U]
            movss   [esi+edi*4], xmm3       ; scrivo su U[u+ i*4] con 0 < i < n

            inc     edi
            jmp forI

    close:   
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

dsecu      equ      8
i1ecu      equ     12
qsecu      equ     16
i2ecu      equ     20
kecu       equ     24
outputecu  equ     28

section .data                             ; Sezione contenente dati inizializzati
section .bss                              ; Sezione contenente dati non inizializzati
section .text                             ; Sezione contenente il codice macchina

euc_dist:
        ; ------------------------------------------------------------
        ; Sequenza di ingresso nella funzione
        ; ------------------------------------------------------------
        push    ebp                         ; salva il Base Pointer
        mov      ebp, esp                    ; il Base Pointer punta al Record di Attivazione corrente
        push    ebx                         ; salva i registri da preservare
        push    esi
        push    edi
        ; ------------------------------------------------------------
        ; legge i parametri dal Record di Attivazione corrente
        ; ------------------------------------------------------------
        xorps xmm4,xmm4
        xor     esi, esi               ;esi=variabile iterativa da 0 a k (fisso)
        mov     edi, [ebp+kecu]        ;edi= k (fisso)
        mov     ecx, [ebp+i1ecu]       ;ecx =i1 var iterativa su ds
        imul    ecx, edi               ;ecx = i1*k
        imul    ecx, 4                 ;ecx= i1*k*4
        mov     eax, [ebp+dsecu]       ;eax= ds
        add     ecx, eax               ;ecx= ds + i1*k*4(fisso)
        mov     edx, [ebp+i2ecu]       ;edx= i2 var iterativa su qs
        imul    edx, edi               ;ecx = i2*k
        imul    edx, 4                 ;edx= i2*k*4
        mov     ebx, [ebp+qsecu]       ;ebx= qs
        add     edx, ebx               ;edx= qs + i2*k*4(fisso)
        mov     eax, [ebp+kecu]        ;eax = k
        sub     eax, 16                ;eax= k-unroll (fisso)
        mov     ebx, [ebp+kecu]        ;ebx = k
        sub     ebx, 8                 ;ebx= k-8 (fisso)
        
    loop_q_ecu:     
            cmp     esi, eax              
            jg      loop_q4_ecu             
            movups  xmm0, [ecx+4*esi]  
            subps   xmm0, [edx+4*esi] 
            movups  xmm1, [ecx+4*esi+16]    
            subps   xmm1, [edx+4*esi+16]
            movups  xmm2, [ecx+4*esi+32]   
            subps   xmm2, [edx+4*esi+32]
            movups  xmm3, [ecx+4*esi+48]   
            subps   xmm3, [edx+4*esi+48]
            mulps   xmm0, xmm0 
            mulps   xmm1, xmm1 
            mulps   xmm2, xmm2
            mulps   xmm3, xmm3
            addps   xmm4, xmm0                 
            addps   xmm4, xmm1
            addps   xmm4, xmm2                 
            addps   xmm4, xmm3          ;xmm4 registro per somma parziale dei valori
            add esi,16
            jmp     loop_q_ecu 

    loop_q4_ecu:
            
            cmp     esi, ebx              
            jg      h_add_ecu           
            movups  xmm0, [ecx+4*esi]  ;[ds + 4*i1*k + j*4]
            subps   xmm0, [edx+4*esi]          
            mulps   xmm0, xmm0         
            addps   xmm4, xmm0          
            movups  xmm0, [ecx+4*esi+16]   
            subps   xmm0, [edx+4*esi+16]          
            mulps   xmm0, xmm0         
            addps   xmm4, xmm0          
            add esi,8
            ; jmp     loop_q4_ecu

    h_add_ecu:       
        haddps xmm4, xmm4      ;riduco la somma a un valore solo
        haddps xmm4, xmm4       ;riduco la somma a un valore solo

    loop_r_ecu:
            cmp     esi, edi
            jge      end_ecu       ;se j == k no loop resto vai a end altrimenti loop_r
            movss   xmm0, [ecx+4*esi] 
            subss   xmm0, [edx+4*esi]         
            mulss   xmm0, xmm0         
            addss   xmm4, xmm0          ;xmm3 registro per somma parziale dei valori   
            inc esi
            jmp loop_r_ecu

    end_ecu:
        mov eax,[ebp+outputecu]
        sqrtss      xmm4, xmm4
        movss     [eax],xmm4      ; eax = risultato
    ; ------------------------------------------------------------
    ; Sequenza di uscita dalla funzione
    ; ------------------------------------------------------------

    pop  edi                  ; ripristina i registri da preservare
    pop  esi
    pop  ebx
    mov  esp, ebp              ; ripristina lo Stack Pointer
    pop  ebp                  ; ripristina il Base Pointer
    ret                    ; torna alla funzione C chiamante