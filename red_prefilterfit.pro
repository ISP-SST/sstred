; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz, ISP
; 
; 
; :Returns:
; 
;    Prefilter model profile as function of wavelength.
; 
; 
; :Params:
; 
;    par : in, type=fltarr
;
;       Profile model parameters.
;       par[0] = 
;       par[1] = 
;       par[2] = 
;       par[3] = 
;       par[4] = 
;       par[5] = 
;       par[6] = 
;
;
;    iwav : in, type=fltarr
;
;       Wavelengths for which the profile is to be returned.
;
;    pref : in, type=float
;
;       Prefilter wavelength.
;   
; 
; 
; :History:
;
;      2016-12-05 : MGL. Added documentation header.
; 
; 
;-
function chromis_prefilter, par, iwav, pref
    iwav1 = iwav - pref
    res = par[0] / (1.d0+((2.d0*(iwav1 - par[2]) / par[3])^2)^par[4]) * (1.0d0 + par[5]*iwav1 + par[6]*iwav1^3)
  return, res
end
