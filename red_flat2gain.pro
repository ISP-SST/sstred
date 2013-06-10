; docformat = 'rst'

;+
; Make (inverse) gain table from flat field (or sum thereof). 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    Mats Löfdahl (MGL), 2008
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    flat : in, type="2D array"
;   
;      A flat field image.
;   
; 
; :Keywords:
; 
;    badthreshold : in, type=float
;   
;      Unsharp masking threshold for bad pixels .
;   
;    mingain : in, type=float
;   
;      Thresholds on the gain itself
;   
;    maxgain : in, type=float
;   
;      Threshold on the gain itself
;   
;    smoothsize : in, type=float
;   
;      Unsharp masking smoothing kernel width.
;   
;    preserve : in, optional, boolean
;   
;      If set, don't zero borders between Sarnoff taps.
;   
;    gain_nozero : out, optional
;   
;      The gain before zeroing the bad pixels.
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-06-09 : Added documentation. MGL.
; 
; 
;-
function red_flat2gain, flat $
                        , badthreshold = badthreshold $
                        , mingain = mingain $
                        , maxgain = maxgain $
                        , smoothsize = smoothparameter $
                        , preserve = preserve $
                        , gain_nozero = gain_nozero

  if(n_elements(badthreshold) eq 0) then badthreshold = 1.0
  if(n_elements(mingain) eq 0) then mingain = 0.1
  if(n_elements(maxgain) eq 0) then maxgain = 4.0
  if(n_elements(smoothparameter) eq 0) then smoothparameter = 7

  g = 0.0*flat
  indx = where(flat ne 0)
  g[indx] = median(flat[indx]) / flat[indx]
  gain_nozero = g

  ;; dgain = g - smooth(g, smoothparameter, /edge_truncate)
  psf = red_get_psf(round(3*smoothparameter), round(3*smoothparameter), smoothparameter, smoothparameter)
  dgain = g - red_convolve(g, psf / total(psf))
  mask = g GE mingain and g LE maxgain and dgain LT badthreshold AND finite(g)

  ker = replicate(1B, [5, 5])
  mask = morph_open(mask, ker)

  idx = where(mask AND finite(g), count, complement= idx1)
  if(count NE n_elements(flat)) then g[idx1] = 0.0
  if(count gt 0) then g[idx] = median(flat[idx])/flat[idx]

  if(~keyword_set(preserve)) then for ii = 1,7 do g[ii*128,*] = 0.0                   
  return, g
end
