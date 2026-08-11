; docformat = 'rst'

;+
; Calculate the centroid/center-of-mass of a 2D array.
;
; Based on http://www.idlcoyote.com/tips/centroid.html
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
;    The centraid pixel coordinates, (xc,yc).
; 
; :Params:
; 
;    array : in, type="2D array"
; 
;      The input array.
; 
; 
; :History:
; 
;    2018-05-04 : MGL. First version.
; 
;-
function red_centroid, array
  
  if size(array, /n_dimensions) ne 2 then begin
    message, 'Array must be two-dimensional. Returning...', /informational
    return, -1
  endif

  s = size(array, /dimensions)
  totalMass = total(array)

  xcm = total( total(array, 2) * indgen(s[0]) ) / totalMass
  ycm = total( total(array, 1) * indgen(s[1]) ) / totalMass

  return, [xcm, ycm]

end
