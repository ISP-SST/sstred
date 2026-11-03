; docformat = 'rst'

;+
; Generates the effective transmission profile of a dual-etalon FPI
; system
;
; Given a set of reflectivities, cavity separations and F-ratio,
; generates the effective transmission profile of a dual-etalon FPI
; system. This routine implements the "conv" approximation described
; in de la Cruz Rodriguez et al. (2026), ignoring phase-errors between
; different rays in the integral of the pupil.
;
; :Categories:
;
;    SST pipeline
;
; Author:
;
;    Jaime de la Cruz Rodriguez, Institute for solar physics
;
;
; :Returns:
;
;    The fpi transmission profile
;
; :Params:
;
;     tw : in, type=dblarr
;  
;          Array with wavelength offsets where the profile will be computed in Angstroms
;
;     w0 : in, type=double
;
;          Central wavelength where the profile will be computed
; 
;     hc : in, type=float
;
;          Nominal high resolution etalon cavity separation in Angstroms
;
;     hr : in, type=float
;
;          High resolution etalon reflectivity 
;
;     lc : in, type=float
;
;          Nominal low resolution etalon cavity separation in Anstroms
;
;     lr : in, type=float
;
;          Low resolution etalon reflectivity 
;
;     Fr : in, type=float
;
;          F-ratio of the FPI system
;
; :History:
;
;    2026-11-03 : JdlCR. Created! 
;
;-
function red_fpi_profile, tw, w0, hc, hr, lc, lr, Fr

  
  ;; --- default refractive index --- 

  n_hr = 1.0d
  n_lr = 1.0d

  
  ;; --- Init cavity separations ---
  
  nhc = long(0.5d0 + hc*n_hr / (w0 * 0.5d0)) 
  hc = nhc * w0 * 0.5d0
        
  nlc = long(0.5d0 + lc*n_lr / (w0 * 0.5d0)) 
  lc = nlc * w0 * 0.5d0
  
  
  ;; --- Integration angles and weights ---

  PI2 = 2.0d0 * 3.1415926535897932384626433832d0
  calp = cos(sqrt([0.04691007703067d,0.23076534494716d $
                   , 0.50d,0.76923465505284d $
                   , 0.95308992296933d]*(0.5d0/Fr)^2)) * PI2
  wng = [0.11846344252809d $
         , 0.23931433524968d $
         , 0.28444444444444d $
         , 0.23931433524968d $
         , 0.11846344252809d]
  wng /= total(wng) ;; make sure they add up to 1.0

  
  ;; --- Finesse ---
  
  fhr = 4.0d * hr / (1.0d - hr)^2
  flr = 4.0d * lr / (1.0d - lr)^2

  
  ;; --- Phase  ----

  ntw = n_elements(tw)
  tr = dblarr(ntw)

  phr = calp * hc
  plr = calp * lc


  ;; --- calculate blue-shift due to inclined rays: lam = w*cos(theta)*lam0 --- 
  
  dtw = total(wng*calp)
  dnorm = total(wng)

  dtw /= PI2 * dnorm
  dtw = (dtw - 1.0d0)*w0
  
  tw1 = tw + w0 + dtw ;; remove the blueshift

  
  ;; --- Integrate over the pupil assuming symmetry --- 

  for n=0, n_elements(calp)-1 do begin
    tr += (1.0d / (1.0d + fhr * (sin(phr[n] / (tw1)))^2)) $
          * (1.0d / (1.0d + flr * (sin(plr[n] / (tw1)))^2)) * wng[n]
  endfor
  
  return, tr
  
end

ntw = 1001
dw = 0.001d0
w0 = 6173.0d0

tw = (dindgen(ntw) - ntw/2)*dw

;;
instrument = 'CRISP' ;; can be also CRISP2 or CHROMIS
p = red_get_FPI_params(instrument, w0)
tr = red_fpi_profile(tw, w0,  p[0], p[1], p[2] , p[3], p[4])

plot, tw, tr, yrange=[0.0,1.0]

end
