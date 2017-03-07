; docformat = 'rst'

;+
; Makes a combined post-momfbd lowpass filter for a cube of images.
; The filter is meant to preserve the signal-dominated part but remove
; high spatial frequency noise from mosaicking. The mosaicking noise
; appears primarily near the axis directions.
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
;    A low-pass filter.
; 
; :Params:
; 
;    ims : in, type="fltarr(sz,sz,Nim) or complexarr(sz,sz,Nim)"
;   
;      The images on which the filter should be based. (But see
;      keyword fft_done)
;   
;   limfreq : in, optional, type=scalar, default="sz/2"
;   
;      The radial coordinate of the diffraction limit in pixels.
;   
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
;   isotropic : in, optional, type=boolean
;
;       Assume SNR is the same in all directions, base cutoff on
;       directions +/- 45 degrees from the axes. Suitable for e.g.
;       granulation images. If this keyword is not set, then take a
;       two-dimensional approach, perhaps more appropriate for images
;       of elongated structures like sunspot penumbrae or
;       chromospheric fibrils.
;
;   doplot : in, optional, type=boolean
;
;       Set this to plot results of various fits.
;
;   filtercube : out, optional, type="fltarr(sz,sz,Nim)"
;
;       The individual filters.
;
;   fouriercube : out, optional, type="complexarr(sz,sz,Nim)"
;
;       The Fourier transforms of the input images. Will not be
;       returned if fft_done is set.
;
;   wiener : in, optional, type=boolean
;
;       If this keyword is set, the filter generated is the classical
;       Wiener filter, P_signal/(P_signal+P_noise). The signal and
;       noise parts will be identified using an assumption that the
;       power spectrum of the image consists of an exponential falloff
;       "signal" part at low to mid frequencies and a constant "noise"
;       level at high frequencies. If not set, the assumption is that
;       the original noise is all filtered away in image restoration
;       and the ramaining noise is a much weaker component. So the
;       filter cutoff will be set by looking for a steep drop in
;       power.
;
;
; :History:
; 
;    2013-11-18 : MGL. First version.
;
;    2013-11-19 : MGL. Added keyword isotropic and code to make a
;                 simpler, circularly symmetrical Wiener filter.
;
;    2013-11-20 : MGL. Now accepts an image cube and optionally
;                 returns the corresponding cube of filters in keyword
;                 filtercube. Fourier transforms of images can be
;                 returned in keyword fouriercube. Also added keyword
;                 doplot.
;
;    2013-11-21 : MGL. Added wiener keyword and implemented
;                 (isotropic) filters based on drop in power.
;
;-
function red_make_lowpass_postfilter, ims, limfreq $
                                      , apodisation = apodisation $
                                      , fft_done = fft_done $
                                      , isotropic = isotropic $
                                      , doplot = doplot $
                                      , filtercube = filtercube $
                                      , fouriercube = fouriercube $
                                      , wiener = wiener $
                                      , pause = pause
  
  dims = size(ims, /dim)
  sz = dims[0]
  if n_elements(dims) lt 3 then Nims = 1 else Nims = dims[2]

  if n_elements(limfreq) eq 0 then limfreq = sz/2

  xx = (findgen(sz/2))

  ;; Require that the image dimensions are equal.
  if sz ne dims[1] then begin
     print, 'red_make_lowpass_postfilter : images are not square.'
     help, ims
     stop
  endif

  if (sz/2)*2 ne sz then  begin
     print, 'red_make_lowpass_postfilter : Need even image dimensions.'
     help, ims
     stop
  endif

  filtercube = fltarr(dims)
  if ~keyword_set(fft_done) and arg_present(fouriercube) then begin
     fouriercube = complexarr(dims)
  endif

  for iim = 0, Nims-1 do begin

     if keyword_set(fft_done) then begin

        ;; Fourier transform done already     
        pow = abs(ims[*, *, iim])^2         ; Power spectrum

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
           
        endelse 

	if arg_present(fouriercube) then begin
        	fouriercube[*, *, iim] = imf
	endif

        pow = abs(imf)^2        ; Power spectrum

     endelse
     

     if keyword_set(isotropic) then begin
        
        ;; Make a circularly symmetric filter 

        ;; Average power in 20 degree wide sectors around +/-45 degrees
        ;; from the axes.
        pow1d = total(exp(red_sectormean(alog(shiftfft(pow)), sz/2 $
                                     , angles = [-45., 45.]*!pi/180. $
                                     , width = 20*!pi/180.) $
                         ), 2)/2.
        mx = max(pow1d)
        mn = min(pow1d)
        
        if keyword_set(doplot) then begin
           window, 0
           cgplot,pow1d, /ylog
        endif

        if keyword_set(wiener) then begin
           ;; Wiener filter, P_signal/(P_signal+P_noise).

           ;; Fit logarithm of 1D power to a model consisting of a
           ;; constant level and an exponential falloff. Exclude
           ;; lowest and highest spatial frequencies

           ;; First find the constant part between ilo_noise and
           ;; ihi_noise
           ihi_noise = limfreq*.99
           minval = min(pow1d[0:ihi_noise])
           ilo_noise = min(where(abs(pow1d[0:ihi_noise]-minval)/minval lt .04))
           Noise_fit = median(pow1d[ilo_noise:ihi_noise])

           if keyword_set(doplot) then begin
              cgplot, xx[ilo_noise]*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'blue', /over
              cgplot, ihi_noise*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'blue', /over
              cgplot, xx, Noise_fit+0*xx, color = 'blue', /over
           endif

           ;; Then the exponential part between ilo_signal and
           ;; ihi_signal
           ilo_signal = round(sz*0.150)
           ihi_signal = ilo_noise*0.7
           Pstart = [0, -1]
           errs = 1.+0*xx[ilo_signal:ihi_signal]
           P = mpfitexpr('P[0]+P[1]*x', xx[ilo_signal:ihi_signal] $
                         , alog(pow1d[ilo_signal:ihi_signal]), errs, Pstart)
           Signal_fit = exp(P[0]+P[1]*xx)

           if keyword_set(doplot) then begin
              cgplot, ilo_signal*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'red', /over
              cgplot, xx[ihi_signal]*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'red', /over
              cgplot, xx, Signal_fit, color = 'red', /over        
              cgplot, xx, Signal_fit+Noise_fit, color = 'green', /over
           endif

           filtercube[*, *, iim] = red_roundmatrix(signal_fit/(signal_fit+noise_fit), sz/2)

        endif else begin
           ;; Non-Wiener. Find sudden drop in power, corresponding to
           ;; a bump in the derivative of the log.

           hicut = round(limfreq*.96)
           locut = round(hicut/3.)
           ld = deriv(smooth(alog(pow1d),15))
           indx = (indgen(hicut))[locut:*]

           ;; Location of the bump
           r = gaussfit(float(indx), ld[indx], gcoeff, nterms = 6)
           rc = gcoeff[1] 

           ;; Make the filter
           filt = aperture(sz/2, rc)
           kern = red_get_psf(sz, sz, sz/20., sz/20.)
           filtercube[*, *, iim] = red_convolve(filt, kern)
           filtercube[*, *, iim] = filtercube[*, *, iim] / max(filtercube[*, *, iim])

           if keyword_set(doplot) then begin
              window, 1
              cgplot, indx, ld[indx], title = 'Image cube index '+strtrim(iim, 2)
              cgplot, indx, r, color = 'red',/over
              cgplot, rc*[1., 1.], !y.crange, color = 'blue',/over
              cgplot, indx, rc-indx, color = 'cyan', /over
              cgplot, (filtercube[sz/2:*, sz/2, iim])*(!y.crange[1]-!y.crange[0]) + !y.crange[0], /over

              window, 0
              cgplot,pow1d, /ylog, title = 'Image cube index '+strtrim(iim, 2)
              cgplot, rc*[1., 1.], 10^!y.crange, color = 'blue',/over
           endif

        endelse                 ; Wiener

     endif else begin

        ;; Make a filter that is not necessarily circularly symmetric
        if keyword_set(wiener) then begin
           ;; Wiener filter, P_signal/(P_signal+P_noise).

           Nangles = 16
           Nangles = 32

           pow1d = exp(red_sectormean(alog(shiftfft(pow)), sz/2, Nangles $
                                      , angles = angles))
           
           filt1d = 0.0*pow1d

           for ia = 0, Nangles-1 do begin
              
              mx = max(pow1d[*, ia])
              mn = min(pow1d[*, ia])
              
              if keyword_set(doplot) then begin
                 window, 4
                 cgplot,pow1d[*, ia], /ylog
              endif

              ;; Fit logarithm of 1D power to a model consisting of a
              ;; constant (?) level and an exponential falloff. Exclude
              ;; lowest and highest spatial frequencies

              ;; First find the constant part between ilo_noise and
              ;; ihi_noise
              ihi_noise = limfreq*.99
              minval = min(pow1d[0:ihi_noise, ia])
              ilo_noise = min(where(abs(pow1d[0:ihi_noise, ia]-minval)/minval lt .04))
              ;;Noise_fit = median(pow1d[ilo_noise:ihi_noise, ia])
              Noise_fit = min(exp((smooth(alog(pow1d[*, ia]), 5))[ilo_noise:ihi_noise]))

              if keyword_set(doplot) then begin
                 cgplot, exp(smooth(alog(pow1d[*, ia]), 5)), color = 'cyan', /over
                 cgplot, xx, minval+0*xx, linestyle = 2, color = 'blue', /over
                 cgplot, xx[ilo_noise]*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'blue', /over
                 cgplot, ihi_noise*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'blue', /over
                 cgplot, xx, Noise_fit+0*xx, color = 'blue', /over
              endif

              ;; Then the exponential part between ilo_signal and
              ;; ihi_signal
              ilo_signal = round(sz*0.150)
              ihi_signal = ilo_noise*0.7
              Pstart = [0, -1]
              errs = 1.+0*xx[ilo_signal:ihi_signal]
              P = mpfitexpr('P[0]+P[1]*x', xx[ilo_signal:ihi_signal] $
                            , alog(pow1d[ilo_signal:ihi_signal, ia]), errs, Pstart)
              Signal_fit = exp(P[0]+P[1]*xx)

              if keyword_set(doplot) then begin
                 cgplot, ilo_signal*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'red', /over
                 cgplot, xx[ihi_signal]*[1, 1], [mn/10, mx*2], linestyle = 3, color = 'red', /over
                 cgplot, xx, Signal_fit, color = 'red', /over        
                 cgplot, xx, Signal_fit+Noise_fit, color = 'green', /over
              endif

              filt1d[0, ia] = signal_fit/(signal_fit+noise_fit)

           endfor               ; ia

           ;;filtercube[*, *, iim] = red_roundmatrix(signal_fit/(signal_fit+noise_fit), sz/2)
           filtercube[*, *, iim] = red_sectormatrix(filt1d, sz/2)

          
        endif else begin
           ;; Non-Wiener. 

           print, '2D non-Wiener filter not implemented yet.'
           stop

           ;; Cutoff levels to try
           Nlevels = 100
           mx = max(pow)
           mn = min(pow)
           ;; Equidistant in log space
           levels = exp(findgen(Nlevels+1)*(alog(mx)-alog(mn))/Nlevels+alog(mn))

           ;; Collect thresholded areas based on the levels
           areas = fltarr(Nlevels)
           for i = 0, Nlevels-1 do areas[i] = total(pow gt levels[i])
           print, min(deriv(areas), iselect)


           if keyword_set(doplot) then begin
              cgplot,areas                        
              cgplot,-deriv(areas)*10,/over,color='blue'
              cgplot,deriv(deriv(areas))*100,/over,color='red'
              
              cgplot, roundmean(alog(pow))
           endif

           ;; ... select level
           
           ;; Make a psf for smoothing the filter 
           fwhm = sz/32.
           psf = red_get_psf(sz, sz, fwhm, fwhm)
           
           filter = pow gt levels[iselect]               ; Threshold
           filter = shift(filter, sz/2, sz/2)            ; Origin in the center
           filter = morph_close(filter,replicate(1,5,5)) ; Close holes and remove isolated pixels
           filter = red_convolve(filter, psf)            ; Smooth the filter
           filter = shift(filter, sz/2, sz/2)            ; Origin in the corner
           
        endelse                 ; Wiener
     endelse                    ; isotropic

     if keyword_set(pause) then stop


  endfor                        ; iim

  if Nims eq 1 then begin
  	return, filtercube
  endif else begin
  	return, min(filtercube, dim = 3)
  endelse

end
