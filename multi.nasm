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

; %include "sseutils.nasm"


data      equ     8
V       equ     12
U       equ     16
cut		equ		20
i		equ		24
k		equ		28
h		equ		32

dim		equ		4
UNROLL		equ		4

section .data                             ; Sezione contenente dati inizializzati
section .bss                              ; Sezione contenente dati non inizializzati
section .text                             ; Sezione contenente il codice macchina

global multi

multi:
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
        mov		ecx, [ebp+i]	; i = var iterativa da 0 a n
		imul	ecx, dim        ;ecx= i*4 byte
        mov     esi, 0  ;esi è j variabile iterativa da 0 a k-4
    loop_q:
        mov     eax, [ebp+k]    ;eax = k
        sub     eax, UNROLL        ;eax= k-4
        cmp     esi, eax            
        jg      h_add               ;se eax =k-4 salta alle istruzioni h_add altrimenti continua

        mov     edi, [ebp+k]        ;edi= k
        imul    ecx, edi            ;ecx = i*k
        imul    ecx, dim            ;ecx= i*k*4
        mov     eax, [ebp+data]     ;eax= ds
        add     ecx, eax            ;ecx= ds + i*k*4

        ;non modificcare ecx ed esi
        movups  xmm0, [ecx+4*esi]   ;[ds + 4*i*k + j*4]


        ; devo riempire il registro xmm2 con 4 valori non consecutivi di V
        mov     edi, 0          ; edi è il contatore fino a 4

        loop_4:
            mov     eax, [ebp+V]    ; eax=V
            mov     ebx, [ebp+h]    ; ebx = h
            imul    ebx, dim        ; ebx= h*4
            imul     ebx, esi        ; ebx= h*4*j
            add     ebx, eax        ; ebx = h*4*j+V
            mov     eax, [ebp+cut]  ; eax= cut
            movss   xmm1, [ebx+eax*4]; [V + 4*j*h + 4*cut] 

            xorps   xmm2, xmm1      ;sposto il valore in xmm2
            shufps	xmm2, xmm2, 147 ;shift di una posizione 

            inc     edi             ;cont++
            inc     esi             ; j++ (0 < j < k-4)

            cmp     edi, 4          ;devo fare 4 iterazioni
            jne     loop_4       
        
        mulps   xmm0, xmm2          ;xmm0 moltiplico i 4 valori di V con quelli di ds
        addps   xmm3, xmm0          ;xmm3 registro per somma parziale dei valori moltiplicati

        jmp     loop_q    

    h_add:       
        
        haddps xmm3, xmm3       ;riduco la somma a un valore solo
        haddps xmm3, xmm3       ;riduco la somma a un valore solo

    loop_r:
        mov     eax, [ebp+k]    ;eax = k
        cmp     esi, eax
        jge      end             ;se j == k no loop resto vai a end altrimenti loop_r

        ; prendo il valore di ds
        movss   xmm0, [ecx+4*esi]   ;[ds + 4*i*k + j*4]

        ;prendo il valore di V
        mov     eax, [ebp+V]    ; eax=V
        mov     ebx, [ebp+h]    ; ebx = h
        imul    ebx, dim        ; ebx= h*4
        imul     ebx, esi        ; ebx= h*4*j
        add     ebx, eax        ; ebx = h*4*j+V
        mov     eax, [ebp+cut]  ; eax= cut
        movss   xmm1, [ebx+eax*4]; [V + 4*j*h + 4*cut]

        mulss   xmm0,xmm1       ; moltiplico il valore di ds per quello di V
        addps     xmm3, xmm0       ; aggiungo alle somme parziali il risultato ottenuto

        inc esi
        jmp loop_r

    end:
        mov     eax, [ebp+U]       ; eax = U
        mov     ebx, [ebp+h]    ; ebx = h
        mov		ecx, [ebp+i]	; ecx è i = var iterativa da 0 a n
        imul     ecx, ebx            ; ecx = i*h
        imul    ecx, dim            ; ecx = i*h*4
        add     ecx, eax            ; ecx = i*h*4 + U
        mov     edx, [ebp+cut]      ; edx = cut
        movss   [ecx+ edx*4], xmm3        ; aggiungo a [U+ i*4*h + cut*4] il valore in xmm3

        ; ------------------------------------------------------------
        ; Sequenza di uscita dalla funzione
        ; ------------------------------------------------------------

        pop	edi                           ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                      ; ripristina lo Stack Pointer
        pop	ebp                           ; ripristina il Base Pointer
        ret                               ; torna alla funzione C chiamante
