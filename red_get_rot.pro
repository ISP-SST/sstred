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
; :Author:
; 
;   Göran Hosinsky (ANA version)
;
;   Jaime de la Cruz Rodriguez (ported to IDL)
; 
; :Returns:
; 
; 
; :Params:
; 
;   az : in, type=float
;   
;     Azimuth angle of Sun in radians
;   
;   el : in, type=float
;   
;     Elevation angle of Sun in radians
;   
;   dec : 
;   
;     Declination angle of Sun in radians.
;   
; 
; :Keywords:
; 
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-05-07 : MGL. Use !dpi since everything else is in double
;                precision. Some cleanup in the comments.
; 
; 
;-
function red_get_rot,az,el,dec

  drrat = !dpi/180.d0           ;deg to rad rate
  lat = 28.758d0*drrat          ;observatory latitude (La Palma)
  tc = 318.0d0*drrat            ;table constant
  
  ;; TC is the table constant in radians a constant depending on which
  ;; observation table is used. TC is about 48 to give the angle
  ;; between the table surface and the N-S direction at the first
  ;; observation table.
  
  ;; According to spherical astronomy the angle between the N-S great
  ;; circle and the vertical great circle in an AZ-EL telescope varies
  ;; as:
  ra1=asin(cos(lat)*sin(az)/cos(dec*drrat))

  ;; In the image plane the angle of the movement in Elevation varies
  ;; as:
  ra=az+(atan(cos(el),sin(el)))+tc-ra1

  return,ra

end
