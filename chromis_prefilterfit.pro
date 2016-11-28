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
; 
;-
function chromis_prefilterfit, par, xl = xl, yl = yl, ispec = ispec, iwav = iwav, pref=pref
  
  iwav1 = iwav - pref
  xl1 = xl - pref
  
  y1 = interpol(yl, xl1+par[1], iwav1) 
  prefc = chromis_prefilter(par, iwav, pref)
  
  res = ispec - (y1 * prefc)

  plot, iwav, ispec, /line
  oplot, iwav, y1*prefc
  oplot, iwav, prefc/par[0] * max(ispec), line=2
  wait, 0.01
  
  return, res

end
