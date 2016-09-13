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
;-
pro red_append, array, data, ifnotthere = ifnotthere

    if ~n_elements(data) then return

    if n_elements(array) eq 0 then array = [ data ] else begin

       if ~keyword_set(ifnotthere) or total(data eq array) eq 0 then array = [ array, data ]

    endelse

end
