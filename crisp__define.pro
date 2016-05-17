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
               camt:'', $              ;
               camr:'', $              ;
               camwb:'', $             ;
               docamt:0B, $            ;
               docamr:0B, $            ;
               docamwb:0B, $           ;
               camttag:'', $           ;
               camrtag:'', $           ;
               camwbtag:'' $           ;
             }

end
