; docformat = 'rst'

;+
; Calculate WCS coordinates OBSGEO-{X,Y,Z} from lat, long, height.
;
; Based on Rots et al. (2015), 2015A&A...574A..36R, WCS paper IV,
; section 4.1.3.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;     Mats Löfdahl, ISP
; 
; 
; :Returns:
;
;     The array [OBSGEO-X,OBSGEO-Y,OBSGEO-Z]
; 
; :Params:
; 
;   B : in, type=float
; 
;      The latitude in degrees, north positive.
;   
;   L : in, type=float
; 
;      The longitude in degrees, east positive.
;   
;   H : in, type=float
; 
;      The altitude in meters.
; 
; 
; :History:
; 
;    2017-03-06 : MGL. First version.
; 
; 
;-
function red_obsgeo, B, L, H

  a = 6378140.d                 ; [m] semi-major axis
  f = 1d/298.2577d              ; "inverse of the inverse flattening"

  d2r = !dpi/180d

  e2 = 2d*f - f*f
  sinB = sin(B*d2r) & cosB = cos(B*d2r)
  sinL = sin(L*d2r) & cosL = cos(L*d2r)
  
  N = a/sqrt(1d - e2*sinB^2)

  X = (N+H) * cosL * cosB
  Y = (N+H) * sinL * cosB
  Z = (N * (1d - e2) + H) * sinB
  
;  print, 'OBSGEO-{B,L,H} =', B, L, H
;  print, 'OBSGEO-{X,Y,Z} =', X, Y, Z

  return, [X, Y, Z]
  
end

tmp = red_obsgeo(28.7305802d, -17.896868d, 2426d)

end

