; docformat = 'rst'

;+
; Apply a distorsion matrix to an image using bilinear or nearest
; neighbor interpolation.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    J. de la Cruz Rodriguez, Institute for Solar Physics
; 
; 
; :Returns:
;
;    Interpolated version of the input array
; 
; :Params:
;    im : in, type=array
;
;      image/2D array to be corrected for distorsions.
;
;    gr : in, type=array
;
;      Distortion matrix from red_gridmatch. It usually has
;      dimensions (2, nTiles_x, nTiles_y)
; 
; :Keywords:
;
;    nearest : in, optional, type=boolean
;   
;       Use nearest neighbor interpolation
;   
;    nthreads: in, optional, type=int
;
;       Number of threads to use 
;
; :History:
; 
;    2020-10-10 : JdlCR. First version.
; 
;-
function red_stretch_linear, im, gr, nthreads = nthreads, nearest=nearest, missing = missing

  d = size(im, /dim)
  nx = long(d[0])
  ny = long(d[1])
  
  dmat = red_get_full_stretch_matrix(nx, ny, double(gr))

  x = dindgen(nx)
  y = dindgen(ny)

  res = dblarr([nx, ny])
  if(n_elements(nthreads) eq 0) then nthreads=long(4) $
  else nthreads = max([long(nthreads),long(1)])

  libfile = red_libfile('creduc.so')


  if(n_elements(missing) eq 0) then missing = mean(im)
  
  
  if(keyword_set(nearest)) then begin
    dum = call_external(libfile, 'nearest2D_wrap', ny, nx, y, x, double(im), ny, nx, dmat[*,*,1], dmat[*,*,0], res, nthreads, missing)
  endif else begin
    dum = call_external(libfile, 'bilint2D_wrap', ny, nx, y, x, double(im), ny, nx, dmat[*,*,1], dmat[*,*,0], res, nthreads, missing)
  endelse
  
  return, res
end
