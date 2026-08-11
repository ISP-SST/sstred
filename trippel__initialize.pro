; docformat = 'rst'

;+
; Class TRIPPEL
;
; :Author:
; 
;    Mats Löfdahl, ISP
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
;   2017-03-13 : MGL. First version.
;
;-
pro trippel::initialize, filename, develop = develop

  ;; Call initialize of the base-class first to load common parameters
  self->RED::initialize, filename , develop = develop
  
  ;; Then load TRIPPEL specific stuff
  self.trippel_dummy = 1
  
end
