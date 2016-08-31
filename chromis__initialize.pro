; docformat = 'rst'

;+
; Class CHROMIS
;
; :Author:
; 
;    Tomas Hillberg
;
;
; :History:
;
;   2016-04-29 : THI. Split class RED into a base-class (instrument
;                independent parts) and derived classes (CRISP/CHROMIS). 
;
;-
pro chromis::initialize, filename

    self->RED::initialize, filename     ; Call initialize of the base-class first to load common parameters

    ; Then load CHROMIS specific stuff
    self.chromis_dummy = 1
    
end