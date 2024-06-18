; docformat = 'rst'

;+
; The minimum distance from a grid for each point. Used for pinhole
; grid fitting.
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; :Returns:
; 
;    An array with the distances between all points and the closest
;    grid point.
; 
; 
; :Params:
; 
;    P : in
; 
; :Keywords:
; 
;    x : in, required, type=array
;   
;       X coordinates of data points.
;   
;    y :  in, required, type=array
;   
;       Y coordinates of data points.
;   
; 
;    Nx : in, required, type=integer
;   
;       The number of grid points in the X direction.
;   
;    Ny : in, required, type=integer
;   
;       The number of grid points in the Y direction. 
;   
;   
;   
; 
; 
; :History:
;  
;    2013-09-16 : MGL. First version.
; 
;    2024-06-18 : MGL. Re-implemented.
;-
function red_fitgrid_deviates, p, X = x, Y = y

  x0    = p[0]
  y0    = p[1]
  dx    = p[2]
  dy    = p[3]
  theta = p[4]

  ;; 1. Subtract (x0,y0) and rotate the data coordinates by -theta

  xr = (x-x0)*cos(-theta) - (y-y0)*sin(-theta)
  yr = (x-x0)*sin(-theta) + (y-y0)*cos(-theta)
  
  ;; 2. Divide by dx and dy
  xr /= dx
  yr /= dy
  
  ;; 3. The fit error is now the distance from integer coordinates
  ;;dev = sqrt( (xr - round(xr))^2 + (yr - round(yr))^2 )
  dev = [ xr - round(xr), yr - round(yr) ]

  return, dev

end

