
%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 32
;vec1:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 32
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
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

global euc_dist_64

;euc_dist_64(rdi=ds, rsi=qs, rdx=k);
;euc_dist_64(xmm0=ds, xmm1=qs, rdi=k);

euc_dist_64:

	push rbp
    mov rbp, rsp
    push rbx							

    xorps   xmm2,xmm2				; store valore finale
    xor     rax, rax               ;rax=variabile iterativa da 0 a k (fisso)

    mov     r10, rdx		; r10 = i      
    sub     r10, 32         ;r10= k-32 (fisso)

    mov     r11, rdx       ;r11 = k
    sub     r11, 16         ;r11= k-16 (fisso)

    vxorps ymm2,ymm2

    loop_q_ecu:     
        cmp     rax, r10             
        jg      loop_q16_ecu 

        vmovups  ymm0, [rdi+rax*4]  ;[ds+i*4]
        vmovups  ymm1, [rsi+rax*4]  ;[qs+i*4]
        vsubps   ymm0, ymm1          
        vmulps   ymm0, ymm0         
        vaddps   ymm2, ymm0 
        vmovups  ymm0, [rdi+rax*4+32] 
        vmovups  ymm1, [rsi+rax*4+32] 
        vsubps   ymm0, ymm1         
        vmulps   ymm0, ymm0         
        vaddps   ymm2, ymm0          
        vmovups  ymm0, [rdi+rax*4+64]  
        vmovups  ymm1, [rsi+rax*4+64] 
        vsubps   ymm0, ymm1           
        vmulps   ymm0, ymm0         
        vaddps   ymm2, ymm0          
        vmovups  ymm0, [rdi+rax*4+96]  
        vmovups  ymm1, [rsi+rax*4+96] 
        vsubps   ymm0, ymm1           
        vmulps   ymm0, ymm0         
        vaddps   ymm2, ymm0          ;ymm2 registro per somma parziale dei valori
        add rax,32
        jmp     loop_q_ecu 

    loop_q16_ecu:
            
        cmp     rax, r11              
        jg      h_add_ecu          ; se i > k-16

        vmovups  ymm0, [rdi+rax*4]  
        vsubps   ymm0, [rsi+rax*4]          
        vmulps   ymm0, ymm0         
        vaddps   ymm2, ymm0          
        vmovups  ymm0, [rdi+rax*4+32]    
        vsubps   ymm0, [rsi+rax*4+32]          
        vmulps   ymm0, ymm0         
        vaddps   ymm2, ymm0          
        add rax,16

    h_add_ecu:       
        vhaddps ymm2, ymm2      ;riduco la somma a un valore solo
        vhaddps ymm2, ymm2       ;riduco la somma a un valore solo
		vhaddps ymm2, ymm2       ;riduco la somma a un valore solo
        vextractf128 xmm2,ymm2,0

    loop_r_ecu:
        cmp     rax, rdx
        jge      end_ecu       ;se j == k no loop resto vai a end altrimenti loop_r
        movss   xmm0, [rdi+rax*4] 
        subss   xmm0, [rsi+rax*4]          ;xmm0 registro per sottrazione parziale dei valori sotratti
        mulss   xmm0, xmm0         ;xmm0 per il quadrato
        addss   xmm2, xmm0          ;xmm3 registro per somma parziale dei valori   
        inc rax
        jmp loop_r_ecu

    end_ecu:
        movss       xmm0,xmm2
        sqrtss      xmm0, xmm0
      
      

	pop	rbx
    mov rsp, rbp
    pop rbp
    ret


global dividiAss64

;dividiAss64( rdi=vett, rsi= numEleVett, rdx = value )
dividiAss64:

	push rbp
    mov rbp, rsp
    push rbx	

	xor	rax, rax	; rax = i
	vbroadcastss ymm0, [rdx] 		; ymm0 contiene il valore che deve dividere
	mov		rbx, rsi				; rbx= size vettore
	sub		rbx, 32					; rbx = size -8

	loop32_div:
		cmp		rax, rbx
		jg		BloopResto_div

		vmovups	ymm1, [rdi+rax*4]		;[vect + i*4]
		vdivps	ymm1, ymm0				; [vect + i*4] / value
		vmovups	[rdi+rax*4], ymm1		; store

		vmovups	ymm1, [rdi+rax*4+32]
		vdivps	ymm1, ymm0
		vmovups	[rdi+rax*4+32], ymm1

		vmovups	ymm1, [rdi+rax*4+64]
		vdivps	ymm1, ymm0
		vmovups	[rdi+rax*4+64], ymm1

		vmovups	ymm1, [rdi+rax*4+96]
		vdivps	ymm1, ymm0
		vmovups	[rdi+rax*4+96], ymm1

		add 	rax, 32
		jmp		loop32_div

	BloopResto_div:
		vextractf128 xmm0,ymm0,0

	loopResto_div:
		cmp		rax, rsi		; compare numElevett
		je 		end_div

		vmovss	xmm1, [rdi+rax*4]
		vdivss 	xmm1,xmm0
		vmovss	[rdi+rax*4], xmm1

		inc		rax
		jmp		loopResto_div

	end_div:
		pop	rbx
		mov rsp, rbp
		pop rbp
		ret
	
global calcolaTAss_64
;calcolaTAss_64( rdi= vect, rsi= numEle) return sommatoria dei quadrati

calcolaTAss_64:
    push rbp
    mov rbp, rsp
    push rbx

    mov     rbx, rsi        ; rbx = numEle
    sub     rbx, 32         ; rbx = numEle - 32
    xor     rax, rax        ; rax = i var iterativa
    vxorps  ymm0, ymm0      ; azzero registro per le somme parziali

    loop32_cal:
        cmp     rax, rbx
        jg      Bloop8_cal

        vmovups ymm1, [rdi+4*rax]   ; prendo [vect + 4*i]
        vmovups ymm2, [rdi+4*rax+32]   ; prendo [vect + 4*i]
        vmovups ymm3, [rdi+4*rax+64]   ; prendo [vect + 4*i]
        vmovups ymm4, [rdi+4*rax+96]   ; prendo [vect + 4*i]

        vmulps  ymm1, ymm1              ; quadrato dei valori
        vmulps  ymm2, ymm2
        vmulps  ymm3, ymm3
        vmulps  ymm4, ymm4

        vaddps  ymm0, ymm1
        vaddps  ymm0, ymm2        ; somme parziali
        vaddps  ymm0, ymm3
        vaddps  ymm0, ymm4

        add     rax, 32
        jmp     loop32_cal

    Bloop8_cal:
        mov     rbx, rsi        ; rbx = numEle
        sub     rbx, 8         ; rbx = numEle - 8

    loop8_cal:
        cmp     rax, rbx
        jg      BloopResto_cal

        vmovups ymm1, [rdi+4*rax]   ; prendo [vect + 4*i]
        vmulps  ymm1, ymm1              ; quadrato dei valori
        vaddps  ymm0, ymm1        ; somme parziali

        add     rax, 8
        jmp     loop8_cal

    BloopResto_cal:
        vhaddps ymm0, ymm0      ;riduco la somma a un valore solo
        vhaddps ymm0, ymm0      
        vhaddps ymm0, ymm0      
        vextractf128 xmm0, ymm0, 0  ; copio somma parziale in xmm0

    loopResto_cal:
        cmp     rax, rsi
        je      end_cal

        vmovss  xmm1, [rdi+4*rax]   ; prendo 1 valore [vect + 4*i]
        vmulss  xmm1, xmm1           ; quadrato dei valore
        vaddss  xmm0, xmm1          ; somma parziale

        inc     rax
        jmp     loopResto_cal

    end_cal:
        pop	rbx
		mov rsp, rbp
		pop rbp
		ret


global aggiornaDatasetAss_64
;aggiornaDatasetAss_64(RDI = ds, RSI = u, RDX = v, RCX = rigaDs, R8 = k)
aggiornaDatasetAss_64:
    
    push rbp
    mov rbp, rsp
    push rbx

    mov     rax, r8                 ; rax = k
    sub     rax, 32                 ; rax = k - 32
    xor     r9, r9                  ; rdi = j
    vbroadcastss    ymm0, [rsi+4*rcx]   ; riproduco [u + 4*i] in tutto ymm0
    mov      rbx, rcx                   ; rbx = i
    imul     rbx, r8                    ; rbx = i*k
    imul     rbx, 4                     ; rbx = i*k*4
    add      rbx, rdi                      ; rbx = ds + i*k*4

    loop32_agg:
        cmp     r9, rax            
        jg      Bloop8_agg

        vmovups ymm1, [rdx+4*r9]       ; [v + 4*j]
        vmulps  ymm1, ymm0              
        vmovups ymm2, [rbx+4*r9]       ; [ds + i*4*k + 4*j]
        vsubps  ymm2, ymm1
        vmovups [rbx+4*r9], ymm2       ; store in ds

        vmovups ymm1, [rdx+4*r9+32]       ; [v + 4*j]
        vmulps  ymm1, ymm0              
        vmovups ymm2, [rbx+4*r9+32]       ; [ds + i*4*k + 4*j]
        vsubps  ymm2, ymm1
        vmovups [rbx+4*r9+32], ymm2       ; store in ds

        vmovups ymm1, [rdx+4*r9+64]       ; [v + 4*j]
        vmulps  ymm1, ymm0              
        vmovups ymm2, [rbx+4*r9+64]       ; [ds + i*4*k + 4*j]
        vsubps  ymm2, ymm1
        vmovups [rbx+4*r9+64], ymm2       ; store in ds

        vmovups ymm1, [rdx+4*r9+96]       ; [v + 4*j]
        vmulps  ymm1, ymm0              
        vmovups ymm2, [rbx+4*r9+96]       ; [ds + i*4*k + 4*j]
        vsubps  ymm2, ymm1
        vmovups [rbx+4*r9+96], ymm2       ; store in ds

        add     r9, 32
        jmp     loop32_agg

    Bloop8_agg:
        mov     rax, r8
        sub     rax, 8                  ; rax = k - 8

    loop8_agg:
        cmp     r9, rax
        jg      BloopResto_agg

        vmovups ymm1, [rdx+4*r9]       ; [v + 4*j]
        vmulps  ymm1, ymm0              
        vmovups ymm2, [rbx+4*r9]       ; [ds + i*4*k + 4*j]
        vsubps  ymm2, ymm1
        vmovups [rbx+4*r9], ymm2       ; store in ds

        add     r9, 8
        jmp     loop8_agg

    BloopResto_agg:
        vmovss  xmm0, [rsi+4*rcx]       ; prendo un valore di U [u + 4*i]
        vshufps xmm0, xmm0, 0           ; riproduco 1° valore in tutto   
    
    loopResto_agg:
        cmp     r9, r8
        jge     end_agg

        vmovss  xmm1, [rdx+4*r9]       ; [v + 4*j]
        vmulss  xmm1, xmm0
        vmovss  xmm2, [rbx+4*r9]       ; [ds + i*4*k + 4*j]
        vsubss  xmm2, xmm1
        vmovss  [rbx+4*r9], xmm2       ; store in ds

        inc     r9
        jmp     loopResto_agg

    end_agg:
        pop	rbx
		mov rsp, rbp
		pop rbp
		ret



global prodMatriceAss_64
;prodMatriceAss_64(RDI = ds, RSI = v, RDX = v, RCX = n, R8 = k )
prodMatriceAss_64:
    push rbp
    mov rbp, rsp
    push rbx

    xor     r9, r9        ; rdi = i variabile iterativa (0, n)

    forI_prod:
        cmp     r9, rcx
        je      end_prod

        mov     rax, r8
        sub     rax, 32         ; rax = k - 32

        mov     rcx, r8         ; rcx = k
        imul    rcx, rdi        ; rcx = k*i
        imul    rcx, 4          ; rcx = k*i*4
        add     rcx, rdi        ; rcx = k*i*4 + ds

        xor     r10, r10        ; r10 = k varibile ite (0,k)
        vxorps  ymm0, ymm0      ; ymm0 = azzero reg somme parziali

        loopJ32_prod:
            cmp     r10, rax
            jmp     BloopJ8_prod

            vmovups ymm1, [rcx+4*r10]       ; ymm1 [ds + 4*i*k + j*4]
            vmulps  ymm1, [rdx+4*r10]      ; ymm1 moltiplicato [v + j*4]
            vaddps  ymm0, ymm1              ; somme parziali




        BloopJ8_prod:





    end_prod:
        pop	rbx
		mov rsp, rbp
		pop rbp
		ret
    

