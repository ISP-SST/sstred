; docformat = 'rst'

;+
; Computes the az/el of the sun in radians given ha in hours and dec in degrees
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    Dick Shine, LMSAL (ANA routine lapalma_azel.ana)
;
;    Jaime de la Cruz Rodriguez (Ported to IDL)
; 
; 
; :returns:
; 
; 
; :Params:
; 
;   ha : 
;   
;   
;   
;   dec : 
;   
;   
;   
;   az : 
;   
;   
;   
;   el : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red_get_azel,ha,dec,az,el

;  latitude of obs.
  drrat=!pi/180.d0              ;deg to rad rate
  la=28.758d0*drrat
  dr=dec*drrat                  ;radian form of declination
  hr=ha*(360.D0/24.D0)*drrat	;and of hour angle
  
;use Ken's formula from gdr.ana
  s1=sin(la)*sin(dr)
  c1=cos(la)*cos(dr)*cos(hr)
  xq=s1+c1
;clamp xq
  xq=xq<1.0
  xq=xq>(-1.0)
  el=asin(xq)
  s1=sin(dr)-sin(la)*xq
  c1=cos(la)*sqrt(1.-xq*xq)
  xq=s1/c1
  xq=xq<1.0
  xq=xq>(-1.0)
  az=acos(xq)
;Ken's formula loses the az sign, we can restore it using the ha sign, at
;least for La Palma
  az=az+(ha gt 0.)*(2.d0*!pi-2.d0*az)
  az=az-!pi
  return
end
