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


dataset      equ     8
V       equ     12
U       equ     16
cut		equ		20
riga		equ		24
k		equ		28
h		equ		32

dim		equ		4

section .data                             ; Sezione contenente dati inizializzati
var     dd      4
section .bss                              ; Sezione contenente dati non inizializzati
section .text                             ; Sezione contenente il codice macchina

global multi2

multi2:
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
        
        mov     eax, [ebp+k]        ; eax = k
        mov		ecx, [ebp+riga]	    ; i = var iterativa da 0 a n
        imul    ecx, eax            ; ecx = i*k
        sub     eax, dim            ; eax= k-4
        imul    ecx, dim            ; ecx= i*k*4
        add     ecx, [ebp+dataset]  ; ecx= ds + i*k*4

        xor     esi, esi            ; esi è j variabile iterativa da 0 a k-4
        xorps   xmm1, xmm1          ; azzero xmm1 per le somme parziali

        ;non modificare ecx, esi, eax 
        ;esi = j
        ;eax= k-4
        ;ecx= ds + i*k*4

    loop_q:
        cmp     esi, eax            ; if(j>k-4)
        jg      h_add               ; se eax =k-4 salta alle istruzioni h_add altrimenti continua

        ; devo riempire il registro xmm0 con 4 valori di ds
        movups  xmm0, [ecx+4*esi]   ; [ds + 4*i*k + j*4]


        ; mov     ebx, [ebp+h]    ; ebx = h
        ; imul    ebx, dim        ; ebx= h*4
        ; imul    ebx, esi       ; ebx= h*4*j
        ; add     ebx, eax        ; ebx =V+ h*4*j

        ; prendo valore V
        mov     ebx, [ebp+V]    ; ebx = V
        movups  xmm2, [ebx+4*esi]       ;[V+4*j]


        ;moltiplico xmm0 per 4 valori di V
        ; mulps   xmm0, [ebx+4*esi]          ;xmm0 moltiplico i 4 valori di V[V + 4*j] con quelli di ds
        mulps   xmm0, xmm2
        addps   xmm1, xmm0          ;xmm1 registro per somma parziale dei valori moltiplicati
        
        
        add     esi, dim              ; esi (j+=4)
        jmp     loop_q    

    h_add:       
        haddps xmm1, xmm1       ;riduco la somma a un valore solo
        haddps xmm1, xmm1       ;riduco la somma a un valore solo
   
    loop_r:
        mov     eax, [ebp+k]        ; eax = k
        cmp     esi, eax
        jge     end                 ; se j == k no loop resto vai a end altrimenti loop_r
        
        ; prendo il valore di ds
        movss   xmm0, [ecx+4*esi]   ; [ds + 4*i*k + j*4]

        ;moltiplico per il valore di V
        ; mov     ebx, [ebp+h]    ; ebx = h
        ; imul    ebx, dim        ; ebx= h*4
        ; imul    ebx, esi       ; ebx= h*4*j
        ; ; mov     eax, [ebp+V]    ; eax=V
        ; add     ebx, eax        ; ebx =V+ h*4*j
        mov     ebx, [ebp+V]
        mulss  xmm0, [ebx+4*esi];[V + *4*j]
        addps   xmm1, xmm0          ; aggiungo alle somme parziali il risultato ottenuto


        inc esi
        jmp loop_r

    end:
        mov		ecx, [ebp+riga]	    ; ecx è i = var iterativa da 0 a n
        imul     ecx, [ebp+h]           ; ecx = i*h
        imul    ecx, dim            ; ecx = i*h*4
        add     ecx, [ebp+U]            ; ecx = i*h*4 + U
        mov     edx, [ebp+cut]      ; edx= cut

        movss   [ecx+ edx*4], xmm1  ; aggiungo a [U+ i*4*h + cut*4] il valore in xmm3m3
        ; ------------------------------------------------------------
        ; Sequenza di uscita dalla funzione
        ; ------------------------------------------------------------

        pop	edi                           ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                      ; ripristina lo Stack Pointer
        pop	ebp                           ; ripristina il Base Pointer
        ret                               ; torna alla funzione C chiamante
