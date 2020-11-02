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

function red_get_full_stretch_matrix, nx, ny, gr, original_size = original_size, only_shifts=only_shifts

  libfile = red_libfile('creduc.so')
  
  d = size(gr, /dim)
  nx1 = long(d[1])
  ny1 = long(d[2])
  
  res_final = dblarr([nx, ny, 2])

  if(n_elements(container_size) eq 2) then begin
    res = dblarr([original_size[0], original_size[1], 2])
    dim = size(res,/dim)

    dum = call_external(libfile, 'ana_stretch_full_matrix_wrap', long(original_size[0]), long(original_size[1]), ny1, nx1, double(gr), res)

    if(keyword_set(only_shifts)) then begin   
      res[*,*,0] -= dindgen(dim[0])#replicate(1.0d0, dim[1])
      res[*,*,1] -= replicate(1.0d0, dim[0])#dindgen(dim[1])
    endif
    
    dim = size(res,/dim)
    dim[0] = min([dim[0], nx])-1
    dim[1] = min([dim[1], ny])-1

    res_final[0:dim[0],0:dim[1],*] = res[0:dim[0],0:dim[1],*]
    
  endif else begin
    dum = call_external(libfile, 'ana_stretch_full_matrix_wrap', long(ny), long(nx), ny1, nx1, double(gr), res_final)

    if(keyword_set(only_shifts)) then begin   
      res_final[*,*,0] -= dindgen(nx)#replicate(1.0d0, ny)
      res_final[*,*,1] -= replicate(1.0d0, nx)#dindgen(ny)
    endif
    
  endelse
 

  return, res_final
end
