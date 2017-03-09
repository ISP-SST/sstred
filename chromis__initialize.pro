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
;
; :History:
;
;   2016-04-29 : THI. Split class RED into a base-class (instrument
;                independent parts) and derived classes (CRISP/CHROMIS). 
;
;   2017-03-09 : MGL. New keyword "develop".
;
;-
pro chromis::initialize, filename, develop = develop

  ;; Call initialize of the base-class first to load common parameters
  self->RED::initialize, filename , develop = develop
  
  ;; Then load CHROMIS specific stuff
  self.chromis_dummy = 1
  
end
