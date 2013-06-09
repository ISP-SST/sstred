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
;    dx : 
;   
;   
;   
;    dy : 
;   
;   
;   
; 
; :Keywords:
; 
;    cubic  : 
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
function red_shift_im, var, dx, dy, cubic = cubic 
  if(n_elements(cubic) eq 0) then cubic = -0.5 

  dim = size(var, /dimension)
                                ;
                                ; get the index of each matrix element
                                ;
  xgrid = findgen(dim[0]) # (fltarr(dim[1]) + 1.0) - total(dx)
  ygrid = (fltarr(dim[0]) + 1.0) # findgen(dim[1]) - total(dy)

                                ;
                                ; Interpolate
                                ; 
  return, interpolate(var, xgrid, ygrid, missi = median(var), cubic = cubic)
end
