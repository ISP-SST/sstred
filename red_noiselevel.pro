; docformat = 'rst'

;+
; Measure the standard deviation of white noise in an image.
;
; Measures the level of the power spectrum of pic, in the area where
; there should be no signal and no `dreadful cross`.
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
;   Noise stddev in counts.
; 
; :Params:
; 
;   pic : in, type=array
; 
;     The image in which to measure the noise level. If Pic is 2D,
;     then return a scalar. If Pic is 3D, assume array has dimensions
;     [sz,sz,Nim] and return array with Nim noise levels.
; 
; :Keywords:
; 
;   fourier : in, optional, type=boolean
;
;     Pic is the Fourier transform of an image.
;
;   limfreq : in, optional, type=float
;
;     The diffraction limit in pixels.
;
;   oversamplingrate : in, optional, type=float
; 
;     The oversampling factor, limfreq/Nyquist.
; 
; :History:
; 
;   1993-09-02 : MGL. First version (in ANA)
; 
;   2019-05-22 : MGL. Incorporated in SSTRED.
;
;-
function red_noiselevel, Pic $
                         , OversamplingRate = OversamplingRate $
                         , LimFreq = LimFreq $
                         , Fourier = Fourier 

  dims = size(Pic, /dim)
  xDimen = dims[0]
  Ndim = size(Pic, /n_dim)      ; Number of dimensions
  if Ndim le 2 then Nim = 1 else Niim = dims[2]

  if keyword_set(OversamplingRate) and keyword_set(LimFreq) then begin
    print, 'NoiseLevel: Please set only one of LimFreq and OversamplingRate.'
    retall
  endif
  if ~keyword_set(OversamplingRate) and ~keyword_set(LimFreq) then begin
    LimFreq = xDimen/2
  endif else if keyword_set(OversamplingRate) then begin
    LimFreq = xDimen/2 / OversamplingRate
  endif

;print, 'limfreq = ', limfreq

  if keyword_set(Fourier) then begin ; argument_set
    f = Pic
  endif else begin
    if Ndim gt 2 then begin
      f = complexarr(xDimen, xDimen, Nim)
      for i = 0, Nim-1 do begin
        f[0, 0, i] = fft(Pic[*, *, i])
      endfor
    endif else f = fft(Pic)
  endelse

  NoiseLevels = fltarr(Nim)

;  a = aperture(xDimen/2,LimFreq)
;  x = x_coord(xDimen/2,1)
;  y = transpose(x)
  sz = xDimen
  x_coord = (findgen(sz, sz) mod sz)-sz/2
  y_coord = transpose(x_coord)
  r_coord = sqrt(x_coord*x_coord+y_coord*y_coord)
  a = r_coord lt limfreq
  mask = (1-a) and (abs(x_coord) gt xDimen/6) and (abs(y_coord) gt xDimen/6)

  for i = 0, Nim-1 do begin

    p = float(f[*, *, i]*conj(f[*, *, i]))

    noise_power = total(mask * shift(p, xDimen/2, xDimen/2)) / total(mask)
    
    NoiseLevels[i] = sqrt(noise_power) * xDimen
    
  endfor

  if Ndim le 2 then return, NoiseLevels[0] else return, NoiseLevels

end                             ; NoiseLevel
