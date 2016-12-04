;; Return an RGB triple representing the color associated with the
;; wavelength(s) lambda given in nm.
function red_WavelengthToRGB, lambda, hex = hex, num = num

  N = n_elements(lambda)

  rgb = fltarr(3, N)
  
  IntensityMax = 255
  Gamma        = 0.80  

  for i = 0, N-1 do begin

     ;; RGB values
     if lambda[i] lt 380 then rgb[*, i] = [0., 0., 0.] $
     else if lambda[i] lt 440 then rgb[*, i] = [-(lambda[i]-440.)/(440.-380.), 0.0, 1.0] $
     else if lambda[i] lt 490 then rgb[*, i] = [0.0, (lambda[i] - 440.)/(490. - 440.), 1.0] $
     else if lambda[i] lt 510 then rgb[*, i] = [0.0, 1.0, -(lambda[i] - 510.) / (510. - 490.)] $
     else if lambda[i] lt 580 then rgb[*, i] = [(lambda[i] - 510.) / (580. - 510.), 1., 0.] $
     else if lambda[i] lt 645 then rgb[*, i] = [1., -(lambda[i] - 645.) / (645. - 580.), 0.] $
     else if lambda[i] lt 780 then rgb[*, i] = [1., 0., 0.] $
     else rgb[*, i] = [0., 0., 0.]
     
     ;; Let the intensity fall off near the vision limits
     if lambda[i] lt 380 then factor = 0.0 $
     else if lambda[i] lt 420 then factor = 0.3 + 0.7*(lambda[i] - 380.) / (420. - 380.) $
     else if lambda[i] lt 700 then factor = 1.0 $
     else if lambda[i] lt 780 then factor = 0.3 + 0.7*(780. - lambda[i]) / (780. - 700.) $
     else factor = 0.0
     
     ;; "Adjust"
     rgb[*, i] = IntensityMax * (rgb[*, i] * Factor)^Gamma
     
  endfor

  rgb = round(rgb)
  
  ;; Make RGB hexadecimal strings.
  if keyword_set(hex) then begin  
     rgbhex = strarr(N)
     for i = 0, N-1 do rgbhex[i] = string(rgb[0, i]*65536L+rgb[1, i]*256L+rgb[2, i],format='(Z06)')
     return, rgbhex
  endif

  ;; Make long integer array, suitable for plotting.
  if keyword_set(num) then begin
     rgbnum = lonarr(N)
     for i = 0, N-1 do rgbnum[i] = rgb[0, i] + rgb[1, i]*256L + rgb[2, i]*65536L
     return, rgbnum
  endif

  ;; Return the triples
  return, rgb

END                             ; WavelengthToRGB
