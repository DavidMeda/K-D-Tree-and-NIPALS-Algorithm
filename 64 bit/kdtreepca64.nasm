; ---------------------------------------------------------
; PageRank con istruzioni AVX a 64 bit
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
;     nasm -f elf64 pagerank64.nasm 
;

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
; Funzione prova
; ------------------------------------------------------------
; global prova

; msg	db 'n:',0
; nl	db 10,0

; prova:
; 		; ------------------------------------------------------------
; 		; Sequenza di ingresso nella funzione
; 		; ------------------------------------------------------------
; 		push		rbp				; salva il Base Pointer
; 		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
; 		pushaq						; salva i registri generali

; 		; ------------------------------------------------------------
; 		; I parametri sono passati nei registri
; 		; ------------------------------------------------------------
; 		; rdi = indirizzo della struct input

; 		; esempio: stampa input->n e di input->k
; 		; rdi contiente l'indirizzo della struttura contenente i parametri
; 		; [rdi] contiene l'indirizzo della stringa con il nome del file
; 		; [rdi+8] contiene l'indirizzo di partenza del data set
; 		; [rdi+16] contiene l'indirizzo di partenza del query set
; 		movsx rax, dword[rdi+24]		; [rdi+16] contiene n
; 		prints msg
; 		printi rax
; 		prints nl
; 		;movsx rax, dword[rdi+28]		; a 4 byte da n si trova k
; 		;printi rax
; 		; ------------------------------------------------------------
; 		; Sequenza di uscita dalla funzione
; 		; ------------------------------------------------------------
		
; 		popaq						; ripristina i registri generali
; 		mov		rsp, rbp			; ripristina lo Stack Pointer
; 		pop		rbp					; ripristina il Base Pointer
; 		ret							; torna alla funzione C chiamante


global euc_dist_64

;euc_dist_64(rdi=ds, rsi=qs, rdx=k);
;euc_dist_64(xmm0=ds, xmm1=qs, rdi=k);

euc_dist_64:
    ; ------------------------------------------------------------
    ; Sequenza di ingresso nella funzione
    ; ------------------------------------------------------------

	push rbp
    mov rbp, rsp
    push rbx							
    ; ------------------------------------------------------------
    ; I parametri sono passati nei registri
    ; ------------------------------------------------------------

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
      
      
    ; ------------------------------------------------------------
    ; Sequenza di uscita dalla funzione
    ; ------------------------------------------------------------

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
	


