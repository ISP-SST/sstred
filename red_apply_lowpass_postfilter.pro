; docformat = 'rst'

;+
; Apply a lowpass filter to a cube of images.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for solar physics.
; 
; 
; :Returns:
; 
;    A cube of low-pass filtered images.
; 
; :Params:
; 
;    ims : in, type="fltarr(sz,sz,Nim) or complexarr(sz,sz,Nim)"
;   
;      The image on which the filter should be based. (But see keyword
;      fft_done)
;   
;    filter : in, type="fltarr(sz,sz)"
;   
;      The low-pass filter to be applied to all images in ims.
; 
; :Keywords:
; 
;    apodisation : in, optional, type=float
;
;       By default, the input image is assumed to be apodised in order
;       to reduce the fft "dreadful cross". But if this keyword is
;       given and non-zero, then apodisation will be done. If the
;       value is unity, a default 1/16 of the image dimension will be
;       a smooth transition from 0 to 1. If the value is smaller than
;       unity, then that fraction will be used instead of the default.
; 
;   fft_done : in, optional, type=boolean
;
;       If this keyword is set, im will be assumed to be the fourier
;       transform of the image rather than the image itself. Then the
;       apodisation keyword is ignored.
;
;   doplot : in, optional, type=boolean
;
;       Set this to plot results of various fits.
;
;   returnfourier : in, optional, type=boolean
;
;       Set this to return the Fourier transforms of the filtered
;       images. 
;
;
; :History:
; 
;    2013-11-20 : MGL. First version.
;
;-
function red_apply_lowpass_postfilter, ims, filter $
                                 , apodisation = apodisation $
                                 , fft_done = fft_done $
                                 , doplot = doplot $
                                 , returnfourier = returnfourier
  
  dims = size(ims, /dim)
  sz = dims[0]
  if n_elements(dims) lt 3 then Nims = 1 else Nims = dims[2]

  ;; Require that the image dimensions are equal.
  if sz ne dims[1] then begin
     print, 'red_make_lowpass_postfilter : images are not square.'
     help, ims
     stop
  endif

  ;; Require that the filter and image dimensions match
  if sz ne (size(filter, /dim))[0] or sz ne (size(filter, /dim))[1] then begin
     print, 'red_apply_lowpass_postfilter : filter and image dimensions do not match'
     stop
  endif

  if keyword_set(returnfourier) then begin
     imsf = complexarr(dims)
  endif else begin
     imsf = fltarr(dims)
  endelse

  ;; Is the filter centered on the corner or on the middle of the array?
  if filter(sz/2, sz/2) gt filter(0, 0) then filt = shift(filter, sz/2, sz/2) else filt = filter

  for iim = 0, Nims-1 do begin

     if keyword_set(fft_done) then begin

        ;; Fourier transform done already     
        imf = ims[*, *, iim]

     endif else begin

        ;; We have to calculate the Fourier transform here.
        if n_elements(apodisation) eq 0 then begin
           ;; Fourier transform the image
           
           imf = fft(ims[*, *, iim])

        endif else begin
           ;; Do apodisation first
           
           if apodisation eq 1 then apodisation =  1/16.
           w = red_taper(sz*[1., 0, apodisation])

           ;;mn = mean(ims[*, *, iim])
           mn = total((1-w)*ims[*, *, iim])/total((1-w))
           imf = fft((ims[*, *, iim]-mn)*w+mn)
           
        endelse                 ; Apodisation
        
     endelse                    ; FFT_done

     if keyword_set(returnfourier) then begin
        imsf[*, *, iim] = imf*filt
     endif else begin     
        imsf[*, *, iim] = fft(imf*filt, /inv)
     endelse                    

  endfor                        ; iim
  
  return, imsf
  
end
