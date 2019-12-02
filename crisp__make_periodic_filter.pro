; docformat = 'rst'

;+
; Make a Fourier filter that can be used to remove periodic artifacts
; caused by polcal.
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
; :Params:
; 
;    prefilter : in, type=string
; 
;      The prefilter for which to make the filter.
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2019-12-02 : MGL. First version.
; 
;-
pro crisp::make_periodic_filter, prefilter

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  pfiles = file_search('polcal/cam*_'+prefilter+'_polcal.fits', count = Nfiles)
  if Nfiles ne 2 then begin
     print, inam + ' : No polcal files for prefilter ' + prefilter
     retall
  endif
  
  ;; Make power spectrum based on sum of polcal components.
  
  polcal1 = readfits(pfiles[0], h = h1) 
  polcal2 = readfits(pfiles[1], h = h2) 
  pp = total(polcal2,1) + total(polcal1,1)

  ;; Make a mask for pixel filling
  ccc = red_histo_gaussfit(pp)
  sigma = ccc[2]
  med = ccc[1]
  mask = ~(pp gt med-sigma*5 and pp lt med+sigma*5)

  pp = rdx_fillpix(pp*(1-mask), mask = mask)

;  w=makewindow([1024,40,50])
;  w = fltarr(1024, 1024)+1.

;med = median(pp)

;  ff = shiftfft(fft((pp-med)*w))
  ff = shiftfft(fft(pp-med))

  ;; Display the power spectrum and let the user click on peak.

  ;; This could be extended to allow for removing multiple peaks. For
  ;; now only removing a single peak is implemented.
  
  tvscl, alog(abs(ff)^2)
  fac = 10.
  sz = 100
  ffc = centerpic(alog(abs(ff)^2), sz = sz)
  tvscl, rebin(ffc, sz*fac, sz*fac, /samp)

  print, inam+ ' : Click on a significant peak (but not the peak at the origin):'
  cursor, xx, yy, /dev

  xx /= fac
  yy /= fac
  xx -= sz/2
  yy -= sz/2
  print, xx, yy

  ;; Define kernel for suppressing the chosen spatial frequencies. 
  ksz = 15
  kernel = fltarr(ksz, ksz)
  for i = 0, ksz-1 do for j = 0, ksz-1 do kernel[i, j] = sqrt((i-ksz/2)^2+(j-ksz/2)^2)
  kernel /= max(kernel*1.2)

  ;; Define the the Fourier domain filter by placing the kernels at
  ;; the clicked position, as well as the position on the other side
  ;; of the origin. 
  filt = fltarr(1024, 1024)+1.
  for i = 0, ksz-1 do for j = 0, ksz-1 do filt[512+round(xx) + i-ksz/2, 512+round(yy) +j-ksz/2] = kernel[i, j]
  for i = 0, ksz-1 do for j = 0, ksz-1 do filt[512-round(xx) + i-ksz/2, 512-round(yy) +j-ksz/2] = kernel[i, j]

  ;; Apply the filter

  ff2 = ff * filt
  ffc2 = centerpic(alog(abs(ff2)^2 >1e-10), sz = sz)
;  tvscl, rebin(ffc2, sz*fac, sz*fac, /samp)

  ;; Display the polcal sum before and after filtering
  pp2 = float(fft(shiftfft(ff2), /inv)) + med
  red_show, [pp, pp2]

  ;; Display the filter
  red_show, filt, w = 2

  ;; Write the filter to the polcal/ directory
  writefits, red_strreplace(pfiles[0], '_polcal', '_fringefilter'), filt

end 

