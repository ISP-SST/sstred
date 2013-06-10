; docformat = 'rst'

;+
; Convolve an image with a kernel (2D) using FFT.
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
;    img : in, type="2D array"
;   
;      The image to be convolved.
;   
;    psf : in, type="2D array"
;   
;      The PSF or kernel to convolve with.
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_convl2d, img, psf
   
  ;; image dimensions
   
  dim = size(img, /dim)
  dps = size(psf, /dim)
  me = median(img)
   
  ;; npad
   
  npadx = dim[0] + dps[0]
  npady = dim[1] + dps[1]
   
  tmp = fltarr(npadx, npady)
  ppsf = fltarr(npadx, npady)
   
  ;; pad
   
  tmp[0:dim[0]-1, 0:dim[1]-1] = img - me
  ppsf[0:dps[0]-1, 0:dps[1]-1] = psf / total(psf)
  ppsf = shift(temporary(ppsf), -dps/2)
   
  ;; Convolve, clip and add median again.
   
  return, float((fft(fft(ppsf, 1) * fft(tmp, 1),-1))[0:dim[0]-1, 0:dim[1]-1]) + me

end
