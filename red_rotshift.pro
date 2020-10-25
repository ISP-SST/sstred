; docformat = 'rst'

;+
; Rotate and shift an image using linear or cubic interpolation.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
;
;    The rotated and shifted image.
; 
; :Params:
; 
;    arr : in, type=array
;
;       The image array to be rotated and shifted.
;
;    angle : in, type=float
;
;       The rotation angle in radians.
;
;    dx : in, type=float
;
;       The X shift in pixels.
;
;    dy : in, type=float
;
;       The Y shift in pixels.       
; 
; 
; :Keywords:
; 
;    background : in, optional, type=float, default="median(arr)"
;   
;       Padding for non-data pixels.
;   
;    cubic : in, optional, type=float
;   
;       Dummy keyword, kept for backwards compatibility
;   
;    linear : in, optional, type=boolean
;   
;       Use linear interpolation, takes precedence over using the
;       cubic keyword.
;   
;     stretch_grid : in, optional, type=array
;
;       stretch grid computed with red_gridmatch (2, nTiles_x, nTiles_y)
;
;    nthreads : in, optional, type=int
;   
;       Number of threads to use during the interpolation.
;
; :History:
; 
;    2020-06-22 : MGL. First version.
;
;    2020-10-01 : JdlCR. Added the stretch_grid correction and
;                        homebrewed multithreaded interpolation.
; 
;-
function red_rotshift, arr, angle, dx, dy $
                       , background = background $
                       , stretch_grid = stretch_grid $
                       , nthreads = nthreads $
                       , nearest = nearest $
                       , stretch_matrix = stretch_matrix

  if n_elements(background) eq 0 then background = median(arr)

  dim = size(arr, /dim)
  xsi = dim[0] * 0.5
  ysi = dim[1] * 0.5

  xg =  dindgen(dim[0])
  yg =  dindgen(dim[1])
  
  xgrid = xg # (dblarr(dim[1]) + 1.0d)
  ygrid = (dblarr(dim[0]) + 1.0d) # yg

  ;; Not needed anymore as we use our own interpolation routines
  ;;xgrid *= float(dim[0]) / max(xgrid)
  ;;ygrid *= float(dim[1]) / max(ygrid)

  if(n_elements(stretch_matrix) ne 0) then begin
    
    dim1 = size(stretch_matrix, /dim)
    
    if((dim1[0] ne dim[0]) and (dim1[1] ne dim[1])) then begin
      print, 'red_rotshift: ERROR, dimensions of the stretch_matrix do not match the image dimensions, STOP!'
      stop
    endif

    smatx = double(stretch_matrix[*,*,0]) - xgrid
    smaty = double(stretch_matrix[*,*,1]) - ygrid

    
  endif else if(n_elements(stretch_grid) ne 0) then begin

    smat = red_get_full_stretch_matrix(dim[0], dim[1], stretch_grid)
    smat[*,*,0] -= xgrid
    smat[*,*,1] -= ygrid 
     
    smatx = smat[*,*,0] 
    smaty = smat[*,*,1] 
  endif else begin
    smatx = xgrid*0
    smaty = ygrid*0
  endelse


  ;; Save some operations by precomputing terms
  
  co = cos(angle[0]) 
  si = sin(angle[0])
  
  xterm = (xgrid - xsi - dx + smatx)
  yterm = (ygrid - ysi - dy + smaty)

  xgrid1 = co * xterm - si * yterm + xsi
  ygrid1 = si * xterm + co * yterm + ysi 
  
  return, red_interpolate2D(xg, yg, arr, xgrid1, ygrid1, nthreads = nthreads, nearest = nearest)
end
