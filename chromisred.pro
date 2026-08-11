; docformat = 'rst'

;+
; Class CHROMIS
;
; :Author:
; 
;    Tomas Hillberg
;
; :Params:
; 
;   filename : in, optional, type=string, default="config.txt"
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
; :History:
;
;   2016-04-29 : THI. Split class RED into a base-class (instrument
;                independent parts) and derived classes (CRISP/CHROMIS). 
;
;   2017-03-09 : MGL. New keyword "develop".
;
;   2021-03-03 : MGL. New keyword no_db.
;
;-
function chromisred, filename, develop = develop, no_db = no_db
  
  return, obj_new('chromis', filename, develop = develop, no_db = no_db)
  
end
