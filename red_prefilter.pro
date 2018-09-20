; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
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
;    lambda : in, type=fltarr
;
;       Wavelengths in Ångström, for which the profile is to be
;       returned.
;
;    cwl : in, type=float
;
;       Prefilter central wavelength in Ångström.
;   
; 
; 
; :History:
;
;      2016-12-05 : MGL. Added documentation header.
;
;      2018-09-20 : MGL. Renamed..
; 
;-
function red_prefilter, par, lambda, cwl

  lambda1 = lambda - cwl
  res = par[0] / (1.d0+((2.d0*(lambda1 - par[2]) / par[3])^2)^par[4]) * (1.0d0 + par[5]*lambda1 + par[6]*lambda1^3)

  return, res

end
