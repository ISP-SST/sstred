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
; :History:
; 
;     2016-05-25 : MGL. First version, code taken from
;                  red_sumfiles.pro. 
; 
; 
;-
pro red_progressbar, i, N, message = message, finished = finished

  if n_elements(message) eq 0 then message = 'Progress'
  bb = string(13B)

  if keyword_set(finished) then begin
     print, bb, message + ' -> ', 100., '%', FORMAT = '(A,A,F5.1,A,$)'
  endif else begin
     norm = 100. / (N - 1.0)     
     print, bb, message + ' -> ', norm * i, '%', FORMAT = '(A,A,F5.1,A,$)'
  endelse
  
end
