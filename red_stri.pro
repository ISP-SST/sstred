; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    var : 
;   
;   
;   
; 
; :Keywords:
; 
;    ni  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_stri, var, ni = ni
  if(n_elements(ni) eq 0) then begin
     res = strcompress(string(var), /remove_all)
  endif else begin
     res = string(var, FORMAT = ni)
  endelse
  return, res
end
