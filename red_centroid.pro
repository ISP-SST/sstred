; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
;    David Fanning
; 
; 
; :returns:
;    The arithmetic mean coordinate of the points in the input (2D) array.
;    ("center of gravity")
; 
; :Params:
; 
;    array : 
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2014-01-20 :  Imported into the crispred namespace
; 
; 
;-
function red_centroid, array
  s = Size(array, /Dimensions)
  totalMass = Total(array)
  xcm = Total( Total(array, 2) * Indgen(s[0]) ) / totalMass
  ycm = Total( Total(array, 1) * Indgen(s[1]) ) / totalMass                                                                        
  return, [xcm, ycm]
end

