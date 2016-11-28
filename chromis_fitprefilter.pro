; docformat = 'rst'

;+
; Prefilter model for CHROMIS.
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
; 
; :Params:
; 
;     par : in
; 
;        The parameters of the prefilter model.
;   
;   
;     iwav : in
;   
;        The wavelengths in [Å].
;
;
;     pref : in
; 
;        The prefilter central wavelength in [Å].
; 
; 
; :History:
; 
;   2016-11-28 : MGL. Split from chromis__fitprefilter.pro and added
;                header. 
; 
; 
;-
function chromis_prefilter, par, iwav, pref

  iwav1 = iwav - pref
  res = par[0] / (1.d0+((2.d0*(iwav1 - par[2]) / par[3])^2)^par[4]) * (1.0d0 + par[5]*iwav1 + par[6]*iwav1^3)

  return, res

end
