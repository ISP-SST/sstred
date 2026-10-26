; docformat = 'rst'

;+
; Add an element to an array, create array if needed.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Tomas Hillberg.
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    array : in, out, type=array
; 
;      The array we are appending.
; 
;    data : in, type=scalar 
; 
;      The new element to append to array.
; 
; :Keywords:
; 
;    count : out, optional, type=integer
;   
;      The number of elements of array on exit.
; 
;    ifnotthere : in, optional, type=boolean
; 
;      If this is set, don't do anything if array already exists.
; 
; :History:
; 
;    2016-09-13 : MGL. New keyword ifnotthere.
; 
;    2016-09-21 : THI. Use logical OR, otherwise the second condition
;                 will always fail for structures.
; 
;    2016-09-13 : MGL. New keyword count.
; 
;-
pro red_append, array, data, ifnotthere = ifnotthere, count = count

  if ~n_elements(data) then return

  if n_elements(array) eq 0 then array = [ data ] else begin
    
    if ~keyword_set(ifnotthere) || total(data eq array) eq 0 then array = [ array, data ]
    
  endelse

  count = n_elements(array)

end
