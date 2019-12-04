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

  Nx = 1024L
  Ny = 1024L
  
  ;; Make a mask for pixel filling
  ccc = red_histo_gaussfit(pp)
  sigma = ccc[2]
  med = ccc[1]
  mask = ~(pp gt med-sigma*5 and pp lt med+sigma*5)

  ;; Clean the mask
  x = Nx/2 + [-3, -2, -1, 0, 1, 2, 3]
  y = Ny/2
  roiPixels = x + y * Nx
  newroi = region_grow(mask,roipixels,threshold=[0,0]) 
  newmask = bytarr(Nx, Ny) + 1
  newmask[newroi] = 0
  
  pp = rdx_fillpix(pp*(1-newmask), mask = newmask)

  ;; Inspect the sum image and decide whether to continue
  red_show, pp
  s = ''
  read, 'Do you want to continue [y/N]? ', s
  if strupcase(strmid(s, 0, 1)) ne 'Y' then return

  ff = shift(fft(pp-med), Nx/2, Ny/2)

  ;; Display the power spectrum and let the user click on peak.

  ;; This could be extended to allow for removing multiple peaks. For
  ;; now only removing a single peak is implemented.
  
  tvscl, alog(abs(ff)^2)
  fac = 10.
  sz = 100
  ffc = red_centerpic(alog(abs(ff)^2), sz = sz)
  red_show, rebin(ffc, sz*fac, sz*fac, /samp)

  ;; Draw axes
  cgarrow, 0, Ny/2-fac*0.75, Nx, Ny/2-fac*0.75, hsize = 0, /solid, color = 'yellow'
  cgarrow, Nx/2-fac*0.75, 0, Nx/2-fac*0.75, Ny, hsize = 0, /solid, color = 'yellow'

  print, inam+ ' : Click on a significant peak (but not at the origin):'
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
  filt = fltarr(Nx, Ny)+1.
  for i = 0, ksz-1 do for j = 0, ksz-1 do filt[Nx/2+round(xx) + i-ksz/2, Ny/2+round(yy) +j-ksz/2] = kernel[i, j]
  for i = 0, ksz-1 do for j = 0, ksz-1 do filt[Nx/2-round(xx) + i-ksz/2, Ny/2-round(yy) +j-ksz/2] = kernel[i, j]

  ;; Apply the filter

  ff2 = ff * filt
  ffc2 = red_centerpic(alog(abs(ff2)^2 >1e-10), sz = sz)

  ;; Display the polcal sum before and after filtering
  pp2 = float(fft(shift(ff2, Nx/2, Ny/2), /inv)) + med
  red_show, [pp, pp2]

  ;; Display the filter
  red_show, filt, w = 2

  ;; Write the filter to the polcal/ directory
  writefits, 'polcal/periodic_filter_'+prefilter+'.fits', filt

end 

