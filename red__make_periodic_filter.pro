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
;   hole_width : in, optional, type=integer, default=15
;   
;      The width of the tapered hole around the spatial frequency
;      to be blocked.
; 
; 
; :History:
; 
;    2019-12-02 : MGL. First version.
; 
;    2019-12-06 : MGL. New keyword hole_width. Write some header info
;                 to the filter file.
; 
;    2022-09-14 : MGL. CRISP --> RED. Take missing pixels into
;                 account. 
; 
;-
pro red::make_periodic_filter, prefilter $
                               , hole_width = hole_width

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(hole_width) eq 0 then hole_width = 15

  ;; Define kernel for suppressing the chosen spatial frequencies. 
  kernel = fltarr(hole_width, hole_width)
  for i = 0, hole_width-1 do $
     for j = 0, hole_width-1 do $
        kernel[i, j] = sqrt((i-hole_width/2)^2+(j-hole_width/2)^2)
  kernel /= max(kernel*1.2)

  pfiles = file_search('polcal/cam*_'+prefilter+'_polcal.fits', count = Nfiles)
  if Nfiles ne 2 then begin
     print, inam + ' : No polcal files for prefilter ' + prefilter
     retall
  endif

  detector1 = (strsplit(file_basename(pfiles[0]), '_', /extract))[0]  
  detector2 = (strsplit(file_basename(pfiles[1]), '_', /extract))[0]  

  ;; Geometrical distortion maps for the two NB cameras
  self -> getalignment, align=align, prefilters=prefilter
  indx = where(align.state2.detector eq detector1, Nalign)
  case Nalign of
    0    : stop                 ; Should not happen!
    1    : amap1 =      align[indx].map
    else : amap1 = mean(align[indx].map, dim = 3)
  endcase
  amap1 /= amap1[2, 2]          ; Normalize
  indx = where(align.state2.detector eq detector2, Nalign)
  case Nalign of
    0    : stop                 ; Should not happen!
    1    : amap2 =      align[indx].map
    else : amap2 = mean(align[indx].map, dim = 3)
  endcase
  amap2 /= amap2[2, 2]          ; Normalize

  amap1_inv = invert(amap1)
  amap2_inv = invert(amap2)
  
  ;; Polcal components.  
  polcal1 = readfits(pfiles[0], h = h1) 
  polcal2 = readfits(pfiles[1], h = h2)
  polcal1 = total(polcal1,1)
  polcal2 = total(polcal2,1)
  polcal1 -= median(polcal1)  
  polcal2 -= median(polcal2)  
  ;; WB orientation:
  polcal1 = rdx_img_project(amap1_inv, polcal1, /preserve)
  polcal2 = rdx_img_project(amap2_inv, polcal2, /preserve)  
  
  pp = polcal1 + polcal2
  
  
  dims = size(pp, /dim)
  
  Nx = dims[0]
  Ny = dims[1]
  
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
  red_show, pp, /scroll
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

  ;; Define the the Fourier domain filter by placing the kernels at
  ;; the clicked position, as well as the position on the other side
  ;; of the origin. 
  filt = fltarr(Nx, Ny)+1.
  for i = 0, hole_width-1 do $
     for j = 0, hole_width-1 do $
        filt[Nx/2+round(xx) + i-hole_width/2 $
             , Ny/2+round(yy) +j-hole_width/2] = kernel[i, j]
  for i = 0, hole_width-1 do $
     for j = 0, hole_width-1 do $
        filt[Nx/2-round(xx) + i-hole_width/2 $
             , Ny/2-round(yy) +j-hole_width/2] = kernel[i, j]

  ;; Apply the filter

  ff2 = ff * filt
  ffc2 = red_centerpic(alog(abs(ff2)^2 >1e-10), sz = sz)

  ;; Display the polcal sum before and after filtering
  pp2 = float(fft(shift(ff2, Nx/2, Ny/2), /inv)) + med
  red_show, [pp, pp2], /scroll

  ;; Display the filter
  red_show, filt, w = 2, /scroll
  
  ;; Write the filter to the polcal/ directory
  red_mkhdr, hdr, filt
  fxaddpar, hdr, 'HOLEWDTH', hole_width, 'Width of "hole" around filtered frequency component'
  fxaddpar, hdr, 'HOLE_X', round(xx), 'Center x frequency coordinate of filtered component'
  fxaddpar, hdr, 'HOLE_Y', round(xx), 'Center y frequency coordinate of filtered component'
  writefits, 'polcal/periodic_filter_remapped_'+prefilter+'.fits', filt

end 

