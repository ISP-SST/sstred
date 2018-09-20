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
; 
;-
function chromis_prefilterfit, par, xl = xl, yl = yl, ispec = ispec, iwav = iwav, pref=pref, w=w
  
  iwav1 = iwav - pref
  xl1 = xl - pref
  
  y1 = interpol(yl, xl1+par[1], iwav1) 
  prefc = chromis_prefilter(par, iwav, pref)
  
  
  return, (ispec - (y1 * prefc)) * w

end
