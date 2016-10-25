; docformat = 'rst'

;+
; Display progress in teminal window.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Params:
; 
;    i : in, type="scalar number"
; 
;       The progress so far.
;
;    N : in, type="scalar number"
;    
;       Normalize i with this number.
; 
; :Keywords:
; 
;    message : in, optional, type=string, default="'Progress'"
;   
;       A message identifying the task for which progress is being
;       made. 
; 
;    finished : in, optional, type=boolean
; 
;       If set, will ignore i and N and write out 100% progress.
;
;    nobar : in, optional, type=boolean
;
;       If set, do not show a status bar.
;
;    barlength : in, optional, type=integer, default=20.
;
;       The length of the status bar in number of characters.
; 
; :History:
; 
;     2016-05-25 : MGL. First version, code taken from
;                  red_sumfiles.pro. 
; 
;     2016-10-10 : MGL. Implemented a real bar. 
; 
;     2016-10-25 : MGL. Deal with N=0 and N=1.
; 
; 
;-
pro red_progressbar, i, N, message = message, finished = finished, nobar = nobar, barlength = barlength

  if n_elements(message) eq 0 then message = 'Progress'
  bb = string(13B)

  if n_elements(barlength) eq 0 then barlength = 20
  
  if keyword_set(finished) then begin
     if keyword_set(nobar) then begin
        print, bb, message + ' -> ' $
               , 100., '%', FORMAT = '(A,A,F5.1,A,$)'
     endif else begin
        print, bb, '['+string(replicate(61B, barlength))+'] ' $
               , 100., '% '+message, FORMAT = '(A,A,F5.1,A,$)'
     endelse
     print, ' '
  endif else begin
     case N of
        0: return
        1: norm = 100.
        else: norm = 100. / (N - 1.0)  
     endcase
     if keyword_set(nobar) then begin
        print, bb, message + ' -> ' $
               , norm * i, '%', FORMAT = '(A,A,F5.1,A,$)'
     endif else begin
        elength = floor(norm*i/100.*barlength)
        mlength = barlength-elength
        bar = ''
        if elength gt 0 then bar += string(replicate(61B, elength))   ; Replicated '='
        if mlength gt 0 then bar += string(replicate(45B, mlength))   ; Replicated '-'
        print, bb, '[' + bar + '] ' $
               , norm*i, '% '+message, FORMAT = '(A,A,F5.1,A,$)'
     endelse
  
  endelse
  
end
