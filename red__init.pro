; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
; 
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
;    develop : in, optional, type=boolean
; 
;       Run in developer mode.
; 
;   
;   
;   
; 
; 
; :History:
;
;   2017-03-09 : MGL. New keyword "develop".
;
; 
; 
;-
function red::init, filename, develop = develop

    ; This function is called implicitly when an instance is created. 
        
    if n_elements(filename) eq 0 then filename = 'config.txt'
    
    if file_test(filename) then self->initialize, filename, develop = develop
    
    return,1
    
end
