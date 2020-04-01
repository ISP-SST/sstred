; docformat = 'rst'

;+
; Wrapper for IDL's rotate() function, adding an inverse flag.
; Also handling image cubes.
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
;    Returns the rotated and/or transposed image(s).
; 
; :Params:
; 
;    array : in, type=array
; 
;       The image to be rotated. For arrays with more than two
;       dimensions, the first two dimensions will be rotated.
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
;   2020-04-01 : MGL. Handle image cubes.
;
;-
function red_rotate, array, direction, inverse = inverse, dir = dir

  if n_elements(direction) eq 0 then direction = 0

  if keyword_set(inverse) xor direction lt 0 then begin
    inv_directions = [0, 3, 2, 1, 4, 5, 6, 7]
    dir = inv_directions[abs( direction)]
  endif else begin
    ;; No inverse, just do it normally. Note that if both
    ;; keyword_set(inverse) and direction lt 0, they cancel. 
    dir = direction
  endelse

  dims = size(array, /dim)
  
  case n_elements(dims) of
    2 : return, rotate(array, dir)
    3 : begin
      if max(direction eq [1, 3, 4, 6]) eq 1 then begin
        ;; X and Y switched
        rotarray = translate(array, dims[[1, 0, 2]])
      endif else begin
        rotarray = array
      endelse
      for i = 0, dims[2]-1 do $
         rotarray[0, 0, i] = rotate(array[*, *, i], dir)
    end
    4 : begin
      if max(direction eq [1, 3, 4, 6]) eq 1 then begin
        ;; X and Y switched
        rotarray = translate(array, dims[[1, 0, 2, 3]])
      endif else begin
        rotarray = array
      endelse
      for i2 = 0, dims[2]-1 do $
         for i3 = 0, dims[3]-1 do $
            rotarray[0, 0, i2, i3] = rotate(array[*, *, i2, i3], dir)
    end
    5 : begin
      if max(direction eq [1, 3, 4, 6]) eq 1 then begin
        ;; X and Y switched
        rotarray = translate(array, dims[[1, 0, 2, 3, 4]])
      endif else begin
        rotarray = array
      endelse
      for i2 = 0, dims[2]-1 do $
         for i3 = 0, dims[3]-1 do $
            for i4 = 0, dims[4]-1 do $
               rotarray[0, 0, i2, i3, i4] = rotate(array[*, *, i2, i3, i4], dir)
    end
    else : stop
  endcase

  return, rotarray
  
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

