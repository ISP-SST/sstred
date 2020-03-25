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
;       rotate(). If negative, the inverse operation is performed.
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
;       Do the inverse operation. Note: if direction is negative and
;       the inverse keyword is set, this is equivalent to calling with
;       abs(direction) and not setting inverse.
; 
; :History:
; 
;   2019-12-09 : MGL. First version. 
; 
;   2020-03-16 : MGL. Accept negative direction.
;
;-
function red_rotate, array, direction, inverse = inverse, dir = dir

  if n_elements(direction) eq 0 then direction = 0

  if keyword_set(inverse) xor direction lt 0 then begin
    inv_directions = [0, 3, 2, 1, 4, 5, 6, 7]
    dir = inv_directions[abs( direction)]
  endif else begin
    ;; No inverse, just do it normally
    dir = direction
  endelse

  return, rotate(array, dir)

end

;; Testing

a = indgen(3, 5)+10

if n_elements(direction) eq 0 then direction = 5

print, a
print
print, rotate(a, direction)
print
print, red_rotate(a, direction)
print
print, red_rotate(red_rotate(a, direction), -direction)
print

print, red_rotate(red_rotate(a, direction), direction, /inverse)


end

