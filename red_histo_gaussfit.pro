; docformat = 'rst'

;+
; Fit a Gaussian to the histogram of input data.
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
; :Returns:
; 
; 
; :Params:
; 
;   data : in, type=array
; 
;      The data, to the histogram of which the Gaussian should be
;      fitted.
; 
; :Keywords:
; 
;   fwlevel : in, optional, type=float
;   
;      Use only data points larger than fwlevel times max(histogram). 
;   
;   nterms : in, optional, type=integer
;   
;      Use this when calling mpfitpeak.
;   
;   show_plots : in, optional, type=boolean
;   
;      Show plot with the fitting results.
; 
; 
; :History:
; 
;   2019-04-02 : MGL. First version.
; 
;-
function red_histo_gaussfit, data $
                             , fwlevel = fwlevel $
                             , nterms = nterms $
                             , show_plots = show_plots

  mn = median(data)
  st = stddev(data)
  
  ;;mn = biweight_mean(images[*,*,ich])
  ;;st = robust_sigma(images[*,*,ich])

  hmax = mn + 3.*st
  hmin = mn - 3.*st
  
  Nbins = 2000
  hh = histogram(data, min = hmin, max = hmax $
                 , Nbins = Nbins, locations = locations)
  binsize = (hmax - hmin) / (Nbins - 1)
  intensities = locations + binsize/2.

  hh = median(hh, 3)            ; Median filter to get rid of single-bin peak

  if n_elements(fwlevel) gt 0 then begin
    ;; Limit the fit to points above fwlwvel * max
    mx = max(hh)
    indx = where(hh ge fwlevel * mx, Nmax)
    if Nmax gt 10 then begin
      intensities = intensities[indx]
      hh = hh[indx]
    endif
  endif

  yfit=mpfitpeak(float(intensities), float(hh), a, nterms = nterms)
  if keyword_set(show_plots) then begin 
    cgplot, intensities, hh, /xstyle, /ystyle
    cgplot, /over, intensities, yfit, color='red'
    cgplot, /over, [0,0]+a[1], [0,1]*max(hh), color='red'
    cgplot, /over, [0,0]+mn,   [0,1]*max(hh), color='cyan'
  endif
  
  return, a
  
end
