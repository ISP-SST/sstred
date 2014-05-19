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
; 
;-
function red_fitgrid_deviates, p, X = x, Y = y, Nx = Nx, Ny = Ny

  grid = red_gridpoints(p[0], p[1], p[2], p[3], p[4], Nx, Ny)

  Npoints = n_elements(x)
  if Npoints ne n_elements(y) then begin
     print, 'Different numbers of x and y coordinates'
     help, x, y
     retall
  endif

  dev = fltarr(Npoints)
  
  for i = 0, Npoints-1 do dev[i] = sqrt(min( (grid.x-x[i])^2 + (grid.y-y[i])^2 ))

  return, dev

end
