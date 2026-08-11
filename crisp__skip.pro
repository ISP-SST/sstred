; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    states : 
;   
;   
;   
;    nremove : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2016-05-19 : Copied (and slightly modified) to the crisp class.
; 
; 
;-
pro crisp::skip, states, nremove

  if( n_elements(nremove) eq 0 ) then nremove=0
  if( nremove eq 0 || n_elements(states) lt 1 ) then return

  os = states[0].tuning
  np = n_elements(states)-1
  
  if( np lt 1 ) then return
                                ;
  for ii = 1L, np  do begin
     if( states[ii].tuning ne os ) then begin
        os = states[ii].tuning
        states[ii:(ii+nremove-1)<np].skip = 1B
     endif
  endfor
                                ;
  return
end
