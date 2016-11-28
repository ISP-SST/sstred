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
function chromis_profile, wav, ang = ang, erh = erh

  if(n_elements(ang) eq 0) then ang = 0.0d0
  if(n_elements(erh) eq 0) then erh = 0.0d0

  w0 = median(wav)
  
  ;; Reflectivities
  if(w0 lt 4010) then begin
    thr = 0.90d0 + erh
    tlr = 0.80d0
  endif else begin
    thr = 0.91d0 + erh
    tlr = 0.8d0
  endelse
  
  ;; Fix cavity separation
  shr = 357.8d4
  nhr = long(0.5d0 + shr / (w0 * 0.5d0))
  hc = nhr * w0 * 0.5d0
  ;;
  slr = 136.9d4
  nlr = long(0.5d0 + slr / (w0 * 0.5d0)) 
  lc = nlr * w0 * 0.5d0

  ;; Finesse
  fhr = 4.d0 * thr / (1.d0 - thr)^2
  flr = 4.d0 * tlr / (1.d0 - tlr)^2

  ;; Phase
  ca = 6.28318530779d0  * cos(ang)
  phr = hc * ca
  plr = lc * ca
  
  ;; Transmission profiles
  hre = 1.d0 / (1.d0 + fhr * (sin(phr / (wav)))^2)
  lre = 1.d0 / (1.d0 + flr * (sin(plr / (wav)))^2)
  
  return, lre*hre

end
