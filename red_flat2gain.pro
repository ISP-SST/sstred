; docformat = 'rst'

;+
; Make (inverse) gain table from flat field (or sum thereof). 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl (MGL), 2008
; 
; 
; :Returns:
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
;    gain_nozero : out, optional
;   
;      The gain before zeroing the bad pixels.
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-06-09 : MGL. Added documentation.
; 
;   2013-07-02 : JdlCR. Fixed a bug that would allow for NaNs in the
;                resulting gaintable
; 
;   2016-05-31 : THI. Remove instrument specific keyword /preserve and
;                use class-methods instead.
; 
;-
function red_flat2gain, flat, $
                        badthreshold = bad, $
                        mingain = min, $
                        maxgain = max, $
                        smoothsize = smoothparameter, $
                        gain_nozero = gain_nozero

  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(smoothparameter) eq 0) then smoothparameter = 7.0d0

  g = median(flat) / flat
  mask1 = ~finite(g)
  pos = where(mask1, count, complement=pos1)
  gain_nozero = red_fillnan(g)
  if(count gt 0) then begin
    g[pos]=0.0
  endif

 ;; dgain = g - smooth(g, smoothparameter, /edge_truncate)
  psf = red_get_psf(round(3*smoothparameter), round(3*smoothparameter) $
                    , double(smoothparameter), double(smoothparameter))
  dgain = g - red_convolve(g, psf / total(psf))
  mask = g GE min and g LE max and dgain LT bad AND finite(g)

  ker = replicate(1B, [5, 5])
  mask = morph_open(mask, ker)
  
  idx = where(mask AND finite(g), count, complement= idx1)
  if(count NE n_elements(flat)) then g[idx1] = 0.0
  if(count gt 0) then g[idx] = median(flat[idx])/flat[idx]

  ;; Jdlcr, Recheck nans
  idx = where(~finite(g), count)
  if(count gt 0) then g[idx] = 0.0

  return, g

end
