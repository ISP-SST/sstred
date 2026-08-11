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
;   x0 : in, type="float or struct"
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
;
; :Keywords:
;
;   grid : in, optional, type="struct or array"
;
;     If a struct, the members are the parameter values. Otherwise,
;     the elements of the array are the parameter values. The
;     parameters themselves are ignored if this keyword is used.
;
; :History:
;  
;    2013-09-16 : MGL. First version.
;  
;    2019-12-09 : MGL. Allow the first parameter to be a struct or an
;                 array, replacing the parameters.
; 
;    2024-06-18 : MGL. Changes due to re-implementation of red_fitgrid.
;    
;-
function red_gridpoints, x0, y0, dx, dy, theta, Nx, Ny, grid = grid

  if n_elements(grid) gt 0 then begin
    if size(grid, /tname) eq 'STRUCT' then begin
      x0    = grid.x0
      y0    = grid.y0
      dx    = grid.dx
      dy    = grid.dy
      theta = grid.theta
      Nx    = grid.nx
      Ny    = grid.ny
    endif else begin
      x0    = grid[0]
      y0    = grid[1]
      dx    = grid[2]
      dy    = grid[3]
      theta = grid[4]
      Nx    = grid[5]
      Ny    = grid[6]
    endelse
  endif
 
  x = ((findgen(2L*Nx)-Nx)*dx) # replicate(1., 2L*Ny)
  y = ((findgen(2L*Ny)-Ny)*dy) ## replicate(1., 2L*Nx)
  
  xr = x*cos(theta) - y*sin(theta) + x0
  yr = x*sin(theta) + y*cos(theta) + y0

  return, {x:xr, y:yr}

end
