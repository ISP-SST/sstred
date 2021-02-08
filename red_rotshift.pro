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
;    Mats Löfdahl & Jaime de la Cruz Rodriguez, Institute for Solar Physics
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
;    unrotated_shifts : in, optional, type=array
;
;        Array of 2 elements or distorsion grid measured before
;        derrotating the images.
;
; :History:
; 
;    2020-06-22 : MGL. First version.
;
;    2020-10-01 : JdlCR. Added the stretch_grid correction and
;                        homebrewed multithreaded interpolation.
;
;    2020-12-09 : JdlCR. When original_dimensions are provided,
;                        The image must be rotated relative to
;                        the center of those dimensions.
;
;    2021-02-07 : JdlCR. Introduced unrotated_shifts, in order to
;                        apply shifts that were measured before derotation.
;-
function red_rotshift, arr, angle, dx, dy $
                       , background = background $
                       , stretch_grid = stretch_grid $
                       , nthreads = nthreads $
                       , nearest = nearest $
                       , stretch_matrix = stretch_matrix $
                       , original_dimensions = original_dimensions $
                       , unrotated_shifts = unrotated_shifts

  if n_elements(background) eq 0 then background = median(arr)
  if n_elements(unrotated_shifts) eq 0 then unrotated_shifts = [0,0]

  
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

    smat = red_get_full_stretch_matrix(dim[0], dim[1], stretch_grid,  /only_shifts)

    smatx = smat[*,*,0] 
    smaty = smat[*,*,1] 
  endif else begin
    smatx = xgrid*0
    smaty = ygrid*0
  endelse

  ;; Are there shifts that were measured before derotation?
  
  if(n_elements(unrotated_shifts) eq 2) then begin
    ushifts_x = unrotated_shifts[0]
    ushifts_y = unrotated_shifts[1]
  endif else if(n_elements(size(unrotated_shifts,/dim)) eq 3) then begin
    ;; Do not use, it does not seem to work!!!
    print, 'red_rotshift: ERROR, providing a distorsion grid as unrotated_shift is untested and it does not seem to work.'
    stop
    ;; I leave the code in case we figure out how to do it
    dum = red_get_full_stretch_matrix(dim[0], dim[1], unrotated_shifts,  /only_shifts)
    ushifts_x = dum[*,*,0]
    ushifts_y = dum[*,*,1]
  endif else begin
    ushifts_x = 0
    ushifts_y = 0
  endelse
  
  
  ;; Save some operations by precomputing terms
  
  co = cos(angle[0]) 
  si = sin(angle[0])
  
  xterm = (xgrid - xsi - dx + smatx)
  yterm = (ygrid - ysi - dy + smaty)

  xgrid1 = co * xterm - si * yterm + xsi + ushifts_x
  ygrid1 = si * xterm + co * yterm + ysi + ushifts_y

  
  return, red_interpolate2D(xg, yg, arr, xgrid1, ygrid1, nthreads = nthreads $
                            , nearest = nearest, missing = background)
end
