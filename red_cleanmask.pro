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
;    gain : 
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
function red_cleanmask, gain

  dims = size(gain, /dim)

  mask = gain ne 0

  oc = bytarr(dims[0]+2, dims[1]+2)
  oc[1, 1] = morph_open(morph_close(mask,replicate(1,5,5)),replicate(1,5,5))
  labels = (label_region(oc, /all_neighbors))[1:dims[0], 1:dims[1]]

  for i = 0, 100 do begin &$
     aa = where(labels eq i, cnt) &$
     if cnt eq 0 then break &$
     print, i
  endfor 
  mg = fltarr(i)
  for ii = 0, i-1 do mg[ii] = mean(gain(where(labels eq ii)))

  mmm = min(mg, minlabel)

  return, mask or (labels eq minlabel)

end
