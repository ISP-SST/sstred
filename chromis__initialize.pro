; docformat = 'rst'

;+
; Class CHROMIS
;
; :Author:
; 
;    Tomas Hillberg
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
;   2016-04-29 : THI. Split class RED into a base-class (instrument
;                independent parts) and derived classes (CRISP/CHROMIS). 
;
;   2017-03-09 : MGL. New keyword "develop".
;
;   2021-03-03 : MGL. New keyword no_db.
;
;-
pro chromis::initialize, filename, develop = develop, no_db = no_db
  
  ;; Call initialize of the base-class first to load common parameters
  self->RED::initialize, filename, develop = develop, no_db = no_db
  
  ;; Then load CHROMIS specific stuff
  self.chromis_dummy = 1
  
end
