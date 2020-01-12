; ---------------------------------------------------------
; PQNN con istruzioni SSE a 32 bit
; ---------------------------------------------------------




dataset      equ     8
V       equ     12
U       equ     16
cut		equ		20              ; parametro colonna fissato per colonna di V e U
riga	equ		24              ; parametro riga i che scorre per le righe di ds
k		equ		28
h		equ		32

dim		equ		4

section .data                       ; Sezione contenente dati inizializzati
section .bss                        ; Sezione contenente dati non inizializzati
section .text                       ; Sezione contenente il codice macchina

;   multi3: moltiplica matrice ds per il vettore colonna di V 
;   il risultato viene scritto in una colle della matrice U

global multi3

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
        
        xor     esi, esi            ; esi è j variabile iterativa da 0 a k-4
        mov     eax, [ebp+k]        ; eax = k
        mov		ecx, [ebp+riga]	    ; i = var iterativa da 0 a n
        imul    ecx, eax            ; ecx = i*k
        sub     eax, dim            ; eax= k-4
        imul    ecx, dim            ; ecx= i*k*4
        add     ecx, [ebp+dataset]  ; ecx= ds + i*k*4
        mov     edx, [ebp+cut]      ; edx= cut

        ;non modificcare ecx, esi, eax, edx
        ;ecx= ds + i*k*4
        ;esi = j
        ;eax= k-4
        ;edx= cut

    loop_q:
        cmp     esi, eax            ; if(j>k-4)
        jg      h_add               ; se eax =k-4 salta alle istruzioni h_add altrimenti continua

        ; devo riempire il registro xmm2 con 4 valori non consecutivi di V
        movups  xmm0, [ecx+4*esi]   ; [ds + 4*i*k + j*4]

        mov     edi, 4              ; edi: cont è il contatore fino a 4
        xorps   xmm2,xmm2 ;per azzerarlo xmm2

        shf:
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
            jg      shf       
            jne     loop_4

        shufps	xmm2, xmm2, 57      ; shift di una posizione a destra
        mulps   xmm0, xmm2          ; xmm0 moltiplico i 4 valori di ds con quelli di V
        addps   xmm3, xmm0          ; xmm3 registro per somma parziale dei valori moltiplicati
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
        imul     ebx, esi           ; ebx= h*4*j
        add     ebx, [ebp+V]        ; ebx = h*4*j+V
        movss   xmm1, [ebx+edx*4]   ; [V + 4*j*h + 4*cut]
        mulss   xmm0,xmm1           ; moltiplico il valore di ds per quello di V
        addps     xmm3, xmm0        ; aggiungo alle somme parziali il risultato ottenuto

        inc esi
        jmp loop_r

    end:
        mov     eax, [ebp+U]        ; eax = U
        mov     ebx, [ebp+h]        ; ebx = h ()tanto un solo valore devo scrivere
        mov		ecx, [ebp+riga]	    ; ecx è i = var iterativa da 0 a n
        imul     ecx, ebx           ; ecx = i*h
        imul    ecx, dim            ; ecx = i*h*4
        add     ecx, eax            ; ecx = i*h*4 + U
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
