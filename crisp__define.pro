; docformat = 'rst'

;+
; Class CRISP
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
pro crisp__define
                                
    struct = { CRISP, inherits RED, $
               crisp_dummy:0B $     ; temporary dummy, move CRISP specific content here
             }

end
