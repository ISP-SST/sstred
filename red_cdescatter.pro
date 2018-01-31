; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    img : 
;   
;   
;   
;    fgain : 
;   
;   
;   
;    fpsf : 
;   
;   
;   
; 
; :Keywords:
; 
;    nthreads  : in, optional, type=integer, default=4
;   
;      Number of threads to use.
;   
;    verbose  : in, optional, type=boolean
;   
;      Set this to get verbose info about what's going on.
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_cdescatter, img, fgain, fpsf, nthreads = nthreads, verbose = verbose
  dim = size(img, /dimension)
  nx = dim[0]
  ny = dim[1]
                                ; Ensure input is a float array
  dimg = fltarr(nx, ny)
                                ;
  if(n_elements(verbose) eq 0) then verbose = 0L
  if(n_elements(nthreads) eq 0) then nthreads = 4L
                                ;
  libfile = red_libfile('creduc.so')
                                ;
  b = call_external(libfile, 'cdescatter', long(nx), $
                    long(ny), float(img), float(fgain), float(fpsf), $
                    dimg, long(nthreads), long(verbose))
                                ;
  return, dimg
end
