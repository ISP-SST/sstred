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
; 
; :Params:
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;   2016-11-28 : MGL. Split from chromis__fitprefilter.pro and added
;                header.
;
;   2016-12-04 : JdlCR. Added weights to remove masked sections.
;
;   2018-09-20 : MGL. Renamed.
; 
;-
function red_prefilterfit, par, xl = xl, yl = yl, spectrum = spectrum, lambda = lambda, pref=pref, w=w
  
  lambda1 = lambda - pref
  xl1 = xl - pref
  
  y1 = interpol(yl, xl1+par[1], lambda1) 
  prefc = chromis_prefilter(par, lambda, pref)
  
  return, (spectrum - (y1 * prefc)) * w

end
