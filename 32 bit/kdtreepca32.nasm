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

;   questa funzione divide ogni elemento di un vettore per un numero dato
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
    mov     ebx, edi                ; ebx = numEle
    sub     edi, 16                 ; edi = numele - 16
    
    xor     esi, esi                    ; esi = i
    movss   xmm0, [ebp+value]       ; valore divisore
    shufps  xmm0, xmm0, 0h

    loop16_div:

        cmp     esi, edi
        jg      loop1_div

        movups  xmm1, [eax+ 4*esi]   ; [vettore + 4*i]
        divps   xmm1, xmm0
        movups  [eax+4*esi], xmm1

        movups  xmm2, [eax+ 4*esi+16]   ; [vettore + 4*i]
        divps   xmm2, xmm0
        movups  [eax+4*esi+16], xmm2

        movups  xmm3, [eax+ 4*esi+32]   ; [vettore + 4*i]
        divps   xmm3, xmm0
        movups  [eax+4*esi+32], xmm3

        movups  xmm4, [eax+ 4*esi+48]   ; [vettore + 4*i]
        divps   xmm4, xmm0
        movups  [eax+4*esi+48], xmm4

        add     esi, 16
        jmp     loop16_div
    
    loop1_div:
        cmp     esi, ebx
        jge     fine_div

        movss  xmm6, [eax+ 4*esi]   ; [vettoer + 4*i]
        divss   xmm6, xmm0
        movss  [eax+4*esi], xmm6

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
nRighe       equ     20
kcoll        equ     24

aggiornaDatasetAss:

    push	ebp                 ; salva il Base Pointer
    mov		ebp, esp            ; il Base Pointer punta al Record di Attivazione corrente
    push	ebx                 ; salva i registri da preservare
    push	esi
    push	edi
    
    xor     edi, edi    ; edi = i var ite (0, n)

    forI_agg:
        mov     eax, [ebp+nRighe]   ; eax = n
        cmp     edi, eax
        je      end_agg

        forJ_agg:
            mov     eax, [ebp+vettoreU] ; eax = u base
            mov     ebx, [ebp+vettoreV] ; eax = v base
            xor     esi, esi            ; esi = j
            movss   xmm0, [eax+4*edi]   ; xmm0= [u + i*4]
            shufps  xmm0, xmm0, 0h      ; riproduco 1° valore in tutto xmm0
            mov     eax, edi            ; eax = i
            imul    eax, 4              ; eax = i*4
            imul    eax, [ebp+kcoll]    ; eax = i*4*k
            add     eax, [ebp+datasetA]    ; eax= ds + i*4*k
            mov     ecx, [ebp+kcoll]    ; ecx = k size
            mov     edx, ecx            ; edx = k
            sub     ecx, 16             ; ecx = k-16

        loop16_agg:

            cmp    esi, ecx
            jg     loopResto_agg

            ; 1° ite
            movups  xmm1, [ebx+4*esi]  ; xmm1 = [v + j*4]
            mulps   xmm1, xmm0
            movups  xmm2, [eax+esi*4]   ; xmm2 = [ds ]
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

            jmp     loop16_agg

        loopResto_agg:
            cmp     esi, edx            ;
            jge     endForI_agg

            movss   xmm1, [ebx+4*esi]   ; xmm1 = [v + j*4]
            mulss   xmm1, xmm0
            movss   xmm2, [eax+esi*4]
            subss   xmm2, xmm1
            movss   [eax+esi*4],xmm2    ; [ds + k*i*4 + j*4]

            inc     esi
            jmp     loopResto_agg
    
    endForI_agg:
        inc     edi
        jmp     forI_agg

    end_agg:
        pop	edi                     ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                ; ripristina lo Stack Pointer
        pop	ebp                     ; ripristina il Base Pointer
        ret                         ; torna alla funzione C chiamante


;   prodMatriceAss: moltiplica matrice ds per il vettore colonna di V 
;   il risultato viene scritto in una colle della matrice U

global prodMatriceAss

datasetProd      equ     8
VProd      equ     12
UProd       equ     16
nProd	equ		20              ; parametro riga i che scorre per le righe di ds
kProd		equ		24

prodMatriceAss:

        push	ebp                 ; salva il Base Pointer
        mov		ebp, esp            ; il Base Pointer punta al Record di Attivazione corrente
        push	ebx                 ; salva i registri da preservare
        push	esi
        push	edi
        
        mov		edx, [ebp+nProd]	    ; edx = n     
        xor     edi, edi             ; i = var iterativa da 0 a n
        mov     ebx, [ebp+VProd]        ; ebx = V        

    forI_prod:
        cmp     edi, edx            ; i >= n
        jge      close

        mov     eax, [ebp+kProd]        ; eax = k
        mov     ecx, edi            ; ecx = i
        imul    ecx, eax            ; ecx = i*k
        sub     eax, 16             ; eax= k-16
        imul    ecx, 4            ; ecx= i*k*4
        add     ecx, [ebp+datasetProd]  ; ecx= ds + i*k*4

        xor     esi, esi            ; esi è j variabile iterativa da 0 a k-4
        xorps   xmm3, xmm3          ; azzero xmm3 per le somme parziali

        ;non modificare ecx, esi, eax, ebx
        ;ecx= ds + i*k*4
        ;esi = j
        ;eax= k-16
        ;ebx = V 

        forJ_prod:
            cmp     esi, eax            ; if(j>k-16)
            jg      hadd_prod               ; se eax =k-16 salta alle istruzioni h_add altrimenti continua

            ;prendo 16 valori consecutivi di ds
            movups  xmm0, [ecx+4*esi]                ; [ds + 4*i*k + j*4]
            movups  xmm4, [ecx+4*esi+16]                ; [ds + 4*i*k + j*4 + 16]
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


            add     esi, 16     ; esi += 16
            jmp     forJ_prod    

        hadd_prod:       
            haddps xmm3, xmm3           ; riduco la somma a un valore solo
            haddps xmm3, xmm3           ; riduco la somma a un valore solo
            mov     eax, [ebp+kProd]        ; eax = k
    
        forJResto_prod:
            cmp     esi, eax
            je     endForI_prod                 ; se j == k no loop resto vai a end altrimenti loop_r
            
            ; prendo il valore di ds
            movss   xmm0, [ecx+4*esi]   ; [ds + 4*i*k + j*4]
            ;prendo il valore di V
            mulss   xmm0,[ebx+esi*4]    ; moltiplico il valore di ds per quello di V [V + 4*j*h + 4*cut]
            addss   xmm3, xmm0          ; aggiungo alle somme parziali il risultato ottenuto

            inc esi
            jmp forJResto_prod

    endForI_prod:
        mov     eax, [ebp+UProd]
        movss   [eax+edi*4], xmm3       ; scrivo su U[u+ i*4] con 0 < i < n

        inc     edi
        jmp     forI_prod

    close:   
        pop	edi                     ; ripristina i registri da preservare
        pop	esi
        pop	ebx
        mov	esp, ebp                ; ripristina lo Stack Pointer
        pop	ebp                     ; ripristina il Base Pointer
        ret                         ; torna alla funzione C chiamante


global euclideanDistanceAss

dsEucl      equ      8
i1eucl      equ     12
qsEucl      equ     16
i2Eucl      equ     20
kEucl       equ     24
outputEucl  equ     28

section .data                             ; Sezione contenente dati inizializzati
section .bss                              ; Sezione contenente dati non inizializzati
section .text                             ; Sezione contenente il codice macchina

euclideanDistanceAss:
        
        push    ebp                         ; salva il Base Pointer
        mov      ebp, esp                    ; il Base Pointer punta al Record di Attivazione corrente
        push    ebx                         ; salva i registri da preservare
        push    esi
        push    edi

        xorps   xmm4,xmm4
        xor     esi, esi               ;esi=variabile iterativa da 0 a k (fisso)
        mov     edi, [ebp+kEucl]        ;edi= k (fisso)
        mov     ecx, [ebp+i1eucl]       ;ecx =i1 var iterativa su ds
        imul    ecx, edi               ;ecx = i1*k
        imul    ecx, 4                 ;ecx= i1*k*4
        mov     eax, [ebp+dsEucl]       ;eax= ds
        add     ecx, eax               ;ecx= ds + i1*k*4(fisso)
        mov     edx, [ebp+i2Eucl]       ;edx= i2 var iterativa su qs
        imul    edx, edi               ;ecx = i2*k
        imul    edx, 4                 ;edx= i2*k*4
        mov     ebx, [ebp+qsEucl]       ;ebx= qs
        add     edx, ebx               ;edx= qs + i2*k*4(fisso)
        mov     eax, edi        ;eax = k
        sub     eax, 16                ;eax= k-16 (fisso)
        
    loop16_eucl:     
            cmp     esi, eax              
            jg       hadd_eucl             
            
            movups  xmm0, [ecx+4*esi] 
            movups  xmm1, [edx+4*esi]
            subps   xmm0, xmm1
            mulps   xmm0, xmm0
            addps   xmm4, xmm0 
            
            movups  xmm0, [ecx+4*esi+16]    
            movups  xmm1, [edx+4*esi+16]
            subps   xmm0, xmm1
            mulps   xmm0, xmm0
            addps   xmm4, xmm0
            
            movups  xmm0, [ecx+4*esi+32]   
            movups  xmm1, [edx+4*esi+32]
            subps   xmm0, xmm1
            mulps   xmm0, xmm0
            addps   xmm4, xmm0
           
            movups  xmm0, [ecx+4*esi+48]   
            movups  xmm1, [edx+4*esi+48]
            subps   xmm0, xmm1
            mulps   xmm0, xmm0
            addps   xmm4, xmm0
            
            add     esi,16
            jmp     loop16_eucl 

    hadd_eucl:       
        haddps xmm4, xmm4      ;riduco la somma a un valore solo
        haddps xmm4, xmm4       ;riduco la somma a un valore solo

    loopResto_eucl:
            cmp     esi, edi
            jge     end_eucl       ;se j == k no loop resto vai a end altrimenti loop_r
            
            movss   xmm0, [ecx+4*esi] 
            movss   xmm1, [edx+4*esi]
            subss   xmm0,  xmm1       
            mulss   xmm0, xmm0         
            addss   xmm4, xmm0          ;xmm4 registro per somma parziale dei valori   
            
            inc     esi
            jmp     loopResto_eucl

    end_eucl:
        mov     eax,[ebp+outputEucl]
        sqrtss  xmm4, xmm4
        movss   [eax],xmm4      ; eax = risultato

    pop  edi                  ; ripristina i registri da preservare
    pop  esi
    pop  ebx
    mov  esp, ebp              ; ripristina lo Stack Pointer
    pop  ebp                  ; ripristina il Base Pointer
    ret                    ; torna alla funzione C chiamante


global prodMatriceTrasAss

datasetTras     equ     8
vettoreVTras    equ     12
vettoreUTras    equ     16
nTras           equ     20
kTras           equ     24


prodMatriceTrasAss:
    push    ebp                         ; salva il Base Pointer
    mov      ebp, esp                    ; il Base Pointer punta al Record di Attivazione corrente
    push    ebx                         ; salva i registri da preservare
    push    esi
    push    edi


    xor     esi, esi                ; esi = i

    fori_tras:
        mov     edx, [ebp+nTras]        ; edx = n
        cmp     esi, edx
        jge      endTras                ; se i >= n

        xor     edi, edi                ; edi = j     
        mov     ecx, [ebp+kTras]        ; ecx= k
        mov     eax, [ebp+datasetTras]  ; eax= ds
        imul    edx, esi, 4            ; edx = i*4
        mov     ebx, [ebp+vettoreUTras] ; ebx = U
        imul     ecx, edx                ; ecx = i*4*K
        add     ecx, eax                ; ecx= i*4*k + ds
        mov     eax, [ebp+vettoreVTras] ; eax = V
        movss   xmm0, [ebx +esi*4]      ; xmm0 c'è il valore di U
        shufps  xmm0, xmm0, 0              ; xmm0 con propagazione valore
        mov     ebx, [ebp+kTras] 
        ; non cambiare esi, edi, ecx (per ds), eax (per V)

        forj_tras:
            mov     edx, [ebp+kTras]    ; edx = k
            sub     edx, 16              ; edx = k-4
            cmp     edi, edx                      
            jg      forjResto_tras           ; se j >= k-4
        
            movups  xmm1, [eax+4*edi]   ;  xmm1 4 valori di V per scrittura
            movups  xmm2, [ecx +4*edi]  ; xmm2 4 valori di ds
            mulps   xmm2, xmm0          ; moltiplico 4 val ds per 1 val di U
            addps   xmm1, xmm2
            
            movups  xmm3, [eax+4*edi+16]   ;  xmm1 4 valori di V per scrittura
            movups  xmm4, [ecx +4*edi+16]  ; xmm2 4 valori di ds
            mulps   xmm4, xmm0          ; moltiplico 4 val ds per 1 val di U
            addps   xmm3, xmm4

            movups  xmm5, [eax+4*edi+32]   ;  xmm1 4 valori di V per scrittura
            movups  xmm6, [ecx +4*edi+32]  ; xmm2 4 valori di ds
            mulps   xmm6, xmm0          ; moltiplico 4 val ds per 1 val di U
            addps   xmm5, xmm6

            movups  xmm2, [eax+4*edi+48]   ;  xmm1 4 valori di V per scrittura
            movups  xmm4, [ecx +4*edi+48]  ; xmm2 4 valori di ds
            mulps   xmm4, xmm0          ; moltiplico 4 val ds per 1 val di U
            addps   xmm2, xmm4

            movups  [eax+4*edi], xmm1
            movups  [eax+4*edi+16], xmm3
            movups  [eax+4*edi+32], xmm5
            movups  [eax+4*edi+48], xmm2

            add     edi, 16
            jmp     forj_tras
    
        forjResto_tras:
            cmp     edi, ebx
            je      Bfori_tras
            movss  xmm1, [eax+4*edi]   ;  xmm1 1 valori di V per scrittura
            movss  xmm2, [ecx +4*edi]  ; xmm2 1 valori di ds
            mulss   xmm2, xmm0          ; moltiplico 4 val ds per 1 val di U
            addss   xmm1, xmm2
            
            movss   [eax+4*edi], xmm1
            
            inc     edi
            jmp     forjResto_tras

    Bfori_tras:
        inc     esi                     ; i++
        jmp     fori_tras


endTras:
    pop  edi                  ; ripristina i registri da preservare
    pop  esi
    pop  ebx
    mov  esp, ebp              ; ripristina lo Stack Pointer
    pop  ebp                  ; ripristina il Base Pointer
    ret                    ; torna alla funzione C chiamante