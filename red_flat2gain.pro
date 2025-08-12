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
;    minflat : in, type=float, default=0.05
;
;      Lower cutoff for flat.
;
;    mingain : in, type=float, default=0.1
;   
;      Thresholds on the gain itself
;   
;    maxgain : in, type=float, default=4.0
;   
;      Threshold on the gain itself
;   
;    smoothsize : in, type=float, default=7.0
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
;   2019-05-14 : MGL. Protect against small values in flat.
; 
;   2025-08-10 : MGL. New keyword minflat.
; 
;-
function red_flat2gain, flat, $
                        badthreshold = bad, $
                        minflat = minflat, $
                        mingain = min, $
                        maxgain = max, $
                        smoothsize = smoothparameter, $
                        gain_nozero = gain_nozero

  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(minflat) eq 0) then minflat = 0.05
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(smoothparameter) eq 0) then smoothparameter = 7.0d0

  indx = where(flat gt median(flat)*minflat, complement = cindx, ncompl = Nc) ; CRISP w/ new cameras
  
  ;;med = median(flat)  
  med = median(flat[indx]) 
  g = med / (flat > (med*1e-5))
  mask1 = ~finite(g)
  pos = where(mask1, count, complement=pos1)
  if arg_present(gain_nozero) then gain_nozero = red_fillnan(g)
  if(count gt 0) then g[pos]=0.0
  if Nc gt 0 then g[cindx] = 0.0

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

;; Data with a bump
fname = '/scratch_local/mats/2025-08-03/CRISP2-test/flats/camXXXI_10.00ms_G00.00_6563.flat.fits'
fname = '/scratch_local/mats/2025-08-03/CRISP2-test/flats/camXXXI_12.00ms_G00.00_8542.flat.fits'

f = readfits(fname)

g = red_flat2gain(f, max = 50., minflat = 0.01)


end

