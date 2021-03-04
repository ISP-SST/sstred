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
;    no_db : in, optional, type=boolean
;
;       Do not use metadata database.
; 
; 
; :History:
;
;   2017-03-09 : MGL. New keyword "develop".
;
;   2021-03-03 : MGL. New keyword no_db.
; 
;-
function red::init, filename, develop = develop, no_db = no_db

  ;; This function is called implicitly when an instance is created. 
  
  if n_elements(filename) eq 0 then filename = 'config.txt'
  
  if file_test(filename) then self->initialize, filename, develop = develop, no_db = no_db
  
  return,1
  
end
