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
; 
; :Keywords:
; 
;    val  : 
;   
;   
;   
;    mask  : 
;   
;   
;   
;    nthreads  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;   2013-08-23 : JdlCR : make sure the C binary gets long(nthreads),
;                        otherwise it can crash.
; 
;-
function red_fillpix, img, val = val, mask = mask, nthreads = nthreads
  inam = 'red_fillpix : '
  if(n_elements(val) eq 0) then val = 0.0001
  if(n_elements(nthreads) eq 0) then nthreads = min([!cpu.TPOOL_NTHREADS,6L]) else nthreads = round(nthreads) 

  libfile = red_libfile('creduc.so')

  dim = size(img, /dim)
  nx = dim[0]
  ny = dim[1]

  ;;
  ;; make mask 
  ;;
  if(n_elements(mask) eq 0) then begin
     mask = bytarr(nx,ny) + 1B
     idx = where(img LT val, count)

     if(count eq 0) then begin
        ;print, inam + 'nothing to do'
        return, img
     endif else mask[idx] = 0B
  endif else mask = byte(temporary(mask))

  res = float(img)
  b = call_external(libfile, 'cfillpix', long(nx), $
                    long(ny), res, mask, long(nthreads))
  
  return, res
end
