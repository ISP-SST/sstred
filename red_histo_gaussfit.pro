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
;   lorentzian  : in, optional, type=boolean
;
;      Fit a Lorentzian kernel to the peak, rather than the default Gaussian. 
;
;   moffat  : in, optional, type=boolean
;
;      Fit a Moffat kernel to the peak, rather than the default Gaussian. 
;
;   nan : in, optional, type=boolean
;
;      Protect against NaN values.
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
;   2021-04-10 : MGL. New keyword nan.
; 
;   2021-04-25 : MGL. New keywords lorentzian and moffat.
; 
;-
function red_histo_gaussfit, data $
                             , lorentzian = lorentzian $
                             , moffat = moffat $
                             , nan = nan $
                             , nterms = nterms $
                             , fwlevel = fwlevel $
                             , show_plots = show_plots
  
  mn = median(data)
  st = stddev(data, nan = nan)
  
  ;;mn = biweight_mean(images[*,*,ich])
  ;;st = robust_sigma(images[*,*,ich])

  hmax = mn + 3.*st
  hmin = mn - 3.*st
  
  Nbins = 2000
  hh = histogram(data, min = hmin, max = hmax $
                 , Nbins = Nbins, locations = locations, nan = nan)
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

  yfit=mpfitpeak(float(intensities), float(hh), a, nterms = nterms $
                 , lorentzian = lorentzian, moffat = moffat)

  if keyword_set(show_plots) then begin 
    cgplot, intensities, hh, /xstyle, /ystyle
    cgplot, /over, intensities, yfit, color='red'
    cgplot, /over, [0,0]+a[1], [0,1]*max(hh), color='red'
    cgplot, /over, [0,0]+mn,   [0,1]*max(hh), color='cyan'
    cgplot, /over, [0,0]+A[2], [0,1]*max(hh), color='green'
    cgplot, /over, [0,0]-A[2], [0,1]*max(hh), color='green'
  endif
  
  return, a
  
end
