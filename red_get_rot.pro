; docformat = 'rst'

;+
; Gives the image rotation as a function of telescope coordinates and
; observation table at the Swedish Solar Observatory, La Palma.
;
; Looking at the projected primary image the rotation is clockwise
; during the morning up the meridian passage and then counterclockwise
; during the afernoon.
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;   Göran Hosinsky (ANA version)
;
;   Jaime de la Cruz Rodriguez (ported to IDL)
; 
; :returns:
; 
; 
; :Params:
; 
;   az : in, type=float
;   
;     Azimuth of Sun in radians
;   
;   el : in, type=float
;   
;     Elevation of Sun in radians
;   
;   dec : 
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
function red_get_rot,az,el,dec

  drrat=!pi/180.d0              ;deg to rad rate
  LAT=28.758d0*drrat            ;observatory latitude (La Palma)
  TC=318.0d0*drrat              ;table constant

; TC is the table constant in radians a constant depending on which
; observation table is used. TC is about 48 to give the angle between
; the table surface and the N-S direction at the first observation
; table.

; According to spherical astronomy the angle between the N-S
; great circle and the vertical great circle in an AZ-EL 
; telescope varies as:

  ra1=asin(cos(LAT)*sin(AZ)/cos(dec*drrat))
; In the image plane the angle of the movement in Elevation
; varies as:

  ra=az+(atan(cos(EL),sin(EL)))+TC-ra1
  return,ra
end
