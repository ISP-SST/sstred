; docformat = 'rst'

;+
; Wrapper for IDL's rotate() function, adding an inverse flag. 
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
;    Returns the rotated and/or transposed array.
; 
; :Params:
; 
;    array : in, type="2D array"
; 
;       The array to be rotated.
; 
;    direction : in, type=integer, default=0
; 
;       Defines the operation to be performed, see documentation for
;       rotate().
; 
; 
; :Keywords:
; 
;    dir : out, optional, type=integer
;   
;       The value of direction that would give the same operation
;       without /inverse.
; 
;    inverse : in, optional, type=boolean
;
;       Do the inverse operation.
; 
; :History:
; 
;   2019-12-09 : MGL. First version. 
;
;-
function red_rotate, array, direction, inverse = inverse, dir = dir

  if n_elements(direction) eq 0 then direction = 0
  
  if keyword_set(inverse) then begin
    inv_directions = [0, 3, 2, 1, 4, 5, 6, 7]
    dir = inv_directions[direction mod 8]
  endif else begin
    dir = direction
  endelse

  return, rotate(array, dir)

end
