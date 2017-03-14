; docformat = 'rst'

;+
; Class SLITJAW
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
pro slitjaw::initialize, filename, develop = develop

  ;; Call initialize of the base-class first to load common parameters
  self->RED::initialize, filename , develop = develop
  
  ;; Then load SLITJAW specific stuff
  self.slitjaw_dummy = 1
  
end
