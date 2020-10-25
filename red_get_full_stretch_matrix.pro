; docformat = 'rst'

;+
;  Computes the full distortion matrix from the output of
;  red_gridmatch, so it can be added to the rotation matrix.
;
;  This routine basically interpolates the grid to the dimensions
;  of the image.
;
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
; 
;    nx : in, type=int
;
;       x-dimension of the output array.
;
;    y : in, type=int
;
;       y-dimension of the output array.
;
;    gr : in, type=array
;
;      Distortion matrix from red_gridmatch. It usually has
;      dimensions (2, nTiles_x, nTiles_y)
; 
; :Keywords:
;
;
; :History:
; 
;    2020-10-10 : JdlCR. First version.
; 
;-

function red_get_full_stretch_matrix, nx, ny, gr

  libfile = red_libfile('creduc.so')

  res = dblarr([nx, ny, 2])

  d = size(gr, /dim)
  nx1 = long(d[1])
  ny1 = long(d[2])
  
  dum = call_external(libfile, 'ana_stretch_full_matrix_wrap', long(ny), long(nx), ny1, nx1, double(gr), res)
                      
  return, res
end
