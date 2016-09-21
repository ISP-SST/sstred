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
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2016-09-13 : MGL. New keyword ifnotthere.
; 
;    2016-09-21 : THI. Use logical OR, otherwise the second condition will always fail for structures.
; 
;-
pro red_append, array, data, ifnotthere = ifnotthere

    if ~n_elements(data) then return

    if n_elements(array) eq 0 then array = [ data ] else begin
    
       if ~keyword_set(ifnotthere) || total(data eq array) eq 0 then array = [ array, data ]
      
    endelse

end
