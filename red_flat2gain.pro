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
;    flat : 
;   
;   
;   
; 
; :Keywords:
; 
;    badthreshold  : 
;   
;   
;   
;    mingain  : 
;   
;   
;   
;    maxgain  : 
;   
;   
;   
;    smoothsize  : 
;   
;   
;   
;    preserve  : 
;   
;   
;   
;    gain_nozero  : 
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
function red_flat2gain, flat, badthreshold = bad, mingain = min, maxgain = max, smoothsize = smoothparameter, preserve = preserve, gain_nozero = gain_nozero
  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(smoothparameter) eq 0) then smoothparameter = 7

  g = 0.0*flat
  indx = where(flat ne 0)
  g[indx] = median(flat[indx]) / flat[indx]
  gain_nozero = g

  ;; dgain = g - smooth(g, smoothparameter, /edge_truncate)
  psf = red_get_psf(round(3*smoothparameter), round(3*smoothparameter), smoothparameter, smoothparameter)
  dgain = g - red_convolve(g, psf / total(psf))
  mask = g GE min and g LE max and dgain LT bad AND finite(g)

  ker = replicate(1B, [5, 5])
  mask = morph_open(mask, ker)

  idx = where(mask AND finite(g), count, complement= idx1)
  if(count NE n_elements(flat)) then g[idx1] = 0.0
  if(count gt 0) then g[idx] = median(flat[idx])/flat[idx]

  if(~keyword_set(preserve)) then for ii = 1,7 do g[ii*128,*] = 0.0                   
  return, g
end
