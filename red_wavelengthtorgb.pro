; docformat = 'rst'

;+
; Translate wavelengths to RGB triplets.
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; :Returns:
;
;    One or several RGB triplets approximating wavelengths within the
;    visible spectrum.
; 
; :Params:
; 
;    lambda : Wavelength in nm.
; 
; :Keywords:
; 
;    hex : in, optional, type=boolean
;
;       Set this to get the triplet in hexadecimal string representation.  
; 
;    num : in, optional, type=boolean
;
;       Set this to get the triplet in numerical form, suitable for
;       use with cgplot.
; 
; 
; :History:
; 
;    2016-12-05 : MGL. Added documentation header. Vectorized.
; 
;    2016-12-13 : MGL. Bugfix in hex and num output.
; 
; 
;-
function red_wavelengthtorgb, lambda, hex = hex, num = num

  N = n_elements(lambda)

  rgb = fltarr(3, N)
  factor = fltarr(N)
  
  IntensityMax = 255
  gamma        = 0.80  

  ;; Intensity falls off near the vision limits:

  indx = where(lambda ge 380 and lambda lt 420, Nindx)
  if Nindx gt 0 then factor[indx] = 0.3 + 0.7*(lambda[indx] - 380.) / (420. - 380.)
  indx = where(lambda ge 420 and lambda lt 700, Nindx)
  if Nindx gt 0 then factor[indx] = 1.0
  indx = where(lambda ge 700 and lambda lt 780, Nindx)
  if Nindx gt 0 then factor[indx] = 0.3 + 0.7*(780. - lambda[indx]) / (780. - 700.) 
  factor = transpose([[factor], [factor], [factor]])

  ;; RGB triplets:

  indx = where(lambda lt 440, Nindx)
  if Nindx gt 0 then begin
    rgb[0, indx] = -(lambda[indx]-440.)/(440.-380.)
    rgb[2, indx] = 1.0
  endif

  indx = where(lambda ge 440 and lambda lt 490, Nindx)
  if Nindx gt 0 then begin
    rgb[1, indx] = (lambda[indx] - 440.)/(490. - 440.) 
    rgb[2, indx] = 1.0
  endif

  indx = where(lambda ge 490 and lambda lt 510, Nindx)
  if Nindx gt 0 then begin
    rgb[1, indx] = 1.0
    rgb[2, indx] = -(lambda[indx] - 510.) / (510. - 490.) 
  endif

  indx = where(lambda ge 510 and lambda lt 580, Nindx)
  if Nindx gt 0 then begin
    rgb[0, indx] = (lambda[indx] - 510.) / (580. - 510.) 
    rgb[1, indx] = 1.0
  endif

  indx = where(lambda ge 580 and lambda lt 645, Nindx)
  if Nindx gt 0 then begin
    rgb[0, indx] = 1.0
    rgb[1, indx] = -(lambda[indx] - 645.) / (645. - 580.)
  endif
  
  indx = where(lambda ge 645 and lambda lt 780, Nindx)
  if Nindx gt 0 then rgb[0, indx] = 1.0

  ;; Apply scaling and gamma
  rgb = round(IntensityMax * (rgb * factor)^gamma)

  ;; Make RGB hexadecimal strings.
  if keyword_set(hex) then $
     return, string(reform(rgb[0, *]*65536L + rgb[1, *]*256L + rgb[2, *]) $
                    , format='(Z06)')
  
  ;; Make long integers, suitable for plotting.
  if keyword_set(num) then $
     return, reform(rgb[0,*] + rgb[1, *]*256L + rgb[2, *]*65536L)
  
  ;; Return the triples
  return, rgb

end                             ; WavelengthToRGB
