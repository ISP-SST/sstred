; docformat = 'rst'

;+
; Model grid point coordinates.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
;    The coordinates of the grid points in a structure with the fields x and y.
; 
; :Params:
; 
;   x0 : in, type=float
;
;     Grid offset in X direction [pixels]. 
; 
;   y0 : in, type=float 
;   
;     Grid offset in Y direction [pixels]. 
; 
;   dx : in, type=float 
; 
;     Grid spacing in X direction [pixels].
;  
;   dy : in, type=float
; 
;     Grid spacing in Y direction [pixels].
; 
;   theta : in, type=float
; 
;     Grid rotation angle [radians].
; 
;   Nx : in, type=integer
;
;     Number of grid points in the X direction. 
; 
;   Ny : in, type=integer
;
;     Number of grid points in the Y direction.
; 
; :History:
;  
;    2013-09-16 : MGL. First version.
;    
;    
;    
;-
function red_gridpoints, x0, y0, dx, dy, theta, Nx, Ny

  x = (findgen(Nx)*dx+x0) # replicate(1., Ny)
  y = (findgen(Ny)*dy+y0) ## replicate(1., Nx)
  
  xr = x*cos(theta) - y*sin(theta)
  yr = x*sin(theta) + y*cos(theta)

  return, {x:xr, y:yr}

end
