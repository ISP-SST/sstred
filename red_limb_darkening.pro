; docformat = 'rst'

;+
; A normalized limb darkening function of wavelength and mu=cos(theta).
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
;    The limb darkening function values.
; 
; :Params:
; 
;     lambda : in, type=float
; 
;       The wavelength in units of meters.
;
;     mu : in, type=fltarr
;
;       The mu values. 
; 
; :Keywords:
; 
;     neckel_coeffs : in, out, optional, type=fltarr
;   
;        The coefficients from Neckel & Labs (1994).
; 
;     analytical : in, optional, type=boolean
; 
;        Get coefficients from Neckel (2005) instead.
; 
; 
; :History:
; 
;    2021-09-06 : MGL. First version, based on private IDL library
;                 function neckel_p5().
; 
;-
function red_limb_darkening, lambda, mu, neckel_coeffs = A, analytical = analytical

  if n_elements(A) ne 6 then begin
    ;; We do analytical or tabulated only if we have not received the
    ;; correct number of coefficients.
    
    A = red_neckel_coefficients(lambda, analytical = analytical)
    
  endif

  P5 = A[5]
  P5 = A[4] + P5*mu
  P5 = A[3] + P5*mu
  P5 = A[2] + P5*mu
  P5 = A[1] + P5*mu
  P5 = A[0] + P5*mu

  return, P5
  
end

lambda = 617e-9
mu = (findgen(100)+1)/100.
n = n_elements(mu)


ld = red_neckel_coefficients(lambda, mu)
ld_anal = red_neckel_coefficients(lambda, mu, /analytical)


cgplot, mu, ld, xtitle = '$\mu$', ytitle = 'Intensity' $
        , xrange = [0.8, 1.01],  yrange = [0.8, 1.01] $
        , color = 'red', title = string(lambda*1e9, format = '(f5.1)')+' nm'
cgplot, /over, mu, ld_anal, color = 'blue'

end
