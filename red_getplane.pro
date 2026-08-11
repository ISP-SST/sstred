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
;    kx : 
;   
;   
;   
;    xx : 
;   
;   
;   
;    yy : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_getplane, kx, xx, yy
  dim = size(kx, /dim)
  res = xx * 0.0
                                ;
  for ii = 0L, dim[1] - 1 do for jj = 0L, dim[1] - 1 do begin
     res+=kx[jj,ii] * yy^float(jj) * xx^float(ii)
  endfor
                                ;
  return, res
end
