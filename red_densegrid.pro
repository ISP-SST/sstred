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
;    xl : 
;   
;   
;   
; 
; :Keywords:
; 
;    thres  : 
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
function red_densegrid, xl, thres = thres
                                ;
  if(n_elements(thres) eq 0) then thres = 0.2
                                ;
  nt = n_elements(xl)
  xld = fltarr(2L*nt-1) + 1.e11
  xld[0:*:2] = xl
                                ;
  for ii = 0L, nt-2 do begin
     sl = (xl[ii] + xl[ii+1])
     dl = abs(xl[ii] - xl[ii+1])
     if(dl ge thres) then xld[2L*ii+1] = 0.5 * sl ; Threshold (only place points if the dlambda > thres)
  endfor
                                ;
  idx = where(xld lt 1.e10, count)
  if(count ge 1) then xld = xld[idx]
                                ;
  return,xld
end
