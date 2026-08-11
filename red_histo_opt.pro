FUNCTION red_histo_opt, image, cutoff, ix, top_only=top, bot_only=bot $
                        , cmin = cmin $
                        , cmax = cmax
;+
; NAME:
;       HISTO_OPT
; PURPOSE:
;       Clip image values which are the CUTOFF brightest or darkest,
;       resp. 
; CATEGORY:
;       
; CALLING SEQUENCE:
;       CLIP_IMAGE = HISTO_OPT ( IMAGE [, CUTOFF [, IX]] [,<keywords>])
; INPUTS:
;       IMAGE : Array with data. may be 1 to 3dim
; OPTIONAL PARAMETERS:
;       IX    : (Output) Contains indices of the clipped values
; KEYWORDS:
;       CMAX     : (Output) Upper cutoff.
;       CMIN     : (Output) Lower cutoff.
;       TOP_ONLY : (Flag) Clip only the upper values
;       BOT_ONLY : (Flag)  "     "  "   lower   "
; OUTPUTS:
;       CLIP_IMAGE : Image with the CUTOFF fraction lowest and highest
;                    values set to the value of the next highest/lowest
;                    point. 
; RESTRICTIONS:
;       Maybe this should be a procedure, as it uses a lot of memory
;       for big arrays. OTOH it is used mainly for displaying, so you
;       wouldn't want to change the real data.
; PROCEDURE:
;       Compute histogram, evaluate the boundaries and return
;       IMAGE>LOW<HIGH 
; MODIFICATION HISTORY:
;       06-Jul-1993  P.Suetterlin, KIS
;       16-Feb-1995  P.Suetterlin, KIS: Take care for float
;                    arrays. Histogram doesn't like them.
;       15-Apr-2021  M.Löfdahl, ISP: Deal with NaNs.
;       10-Jun-2024  M.Löfdahl, ISP: New keywords cmin, cmax.
;-

  on_error, 2

  IF n_params() EQ 0 THEN BEGIN
    message, 'Usage: RESULT = HISTO_OPT ( IMAGE [,CUTOFF] )', /cont
    return, undefined
ENDIF

IF n_params() LT 2 THEN cutoff = 1e-3

s = size(image)
 ;;;
 ;;; If the image is in a float format, then histogram() doesn't know
 ;;; what to do. In that case, convert to fix. But then you have to be
 ;;; shure that the range is ok (especially for normalized images with
 ;;; a range from 0. to 1.). 
 ;;;
indx = where(finite(image), complement = indx_missing, Ncomplement = Nmissing)
IF s(s(0)+1) GT 3 THEN BEGIN
  fak = 10000./(max(image[indx], min = hmin)-hmin)
  h = histogram(fix((image[indx]-hmin)*fak), /nan)
ENDIF ELSE BEGIN
  h = histogram(image[indx], /nan)
  hmin = min(image[indx])
  fak = 1
ENDELSE

nh = n_elements(h)
 ;;;
 ;;; Integrate the histogram so that h(i) holds the number of points
 ;;; with equal or lower intensity.
 ;;;
FOR i = 1l, nh-1 DO h(i) = h(i)+h(i-1)
 ;;;
 ;;; and normalize it to unity
 ;;;
h = float(h)/h(nh-1)
 ;;;
 ;;; As CUTOFF is in percent and h is normalized to unity,
 ;;; cmin/cmax are the indices of the point where the number of pixels
 ;;; with lower/higher intensity reach the given limit. This has to be
 ;;; converted to a real image value by dividing by the scalefactor
 ;;; FAK and adding the min value of the image
 ;;;
cmin = max(where(h LE cutoff))/fak+hmin
cmax = min(where(h GE (1.-cutoff)))/fak+hmin
 ;;;
 ;;; Where is slow. Only compute if requested.
 ;;;
IF n_params() EQ 3 THEN ix = where((image LE cmin) OR (image GE cmax))

if Nmissing gt 0 then begin
  if keyword_set(top) then $
     return, image < cmax $
  else if keyword_set(bot) then $
     return, image > cmin $
  else $
     return, image > cmin < cmax
endif else begin
  outimage = image
  if keyword_set(top) then $
     outimage <= cmax $
  else if keyword_set(bot) then $
     outimage >= cmin $
  else $
     outimage = outimage > cmin < cmax
  outimage[indx_missing] = !values.f_nan
  return, outimage
endelse

end
