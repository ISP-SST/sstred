; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    img : 
;   
;   
;   
;    psf : 
;   
;   
;   
; 
; :Keywords:
; 
;    nthreads  : 
;   
;   
;   
;    verbose  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_cconvolve, img, psf, nthreads = nthreads, verbose = verbose
                                ;
  dim = size(img, /dimension)
  nx = dim[0]
  ny = dim[1]
                                ;
  dim = size(psf, /dimension)
  nx1 = dim[0]
  ny1 = dim[1]
                                ;
  co_img = fltarr(nx, ny)
                                ;
  if(n_elements(verbose) eq 0) then verbose = 0L
  if(n_elements(nthreads) eq 0) then nthreads = 1L
  assign                        ;
  dir=getenv('CREDUC')
                                ;
  dum = call_external(dir+'/creduc.so', 'cconvolve', long(nx), long(ny), $
                      long(nx1), long(ny1), float(img), float(psf), co_img, $
                      long(nthreads), long(verbose))
                                ;
  return, co_img
end
