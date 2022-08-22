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
;   2016-05-19 : THI. Define state structures
;
;   2017-06-28 : MGL. Added a few states from CHROMIS_STATE. Removed
;                cam{wb,r,t}-based states.
;
;   2022-07-29 : MGL. Define the crisp2red class.
;
;-
pro crisp2__define

  struct = { CRISP2, inherits RED, $
             crisp2_dummy:0B $  ; temporary dummy, move CRISP2 specific content here
           }

  nb = { CRISP2_STATE, $
         inherits RED_STATE, $
         scannumber:0, $
         framenumber:-1L, $
         tuning:'', $
         prefilter:'', $
         pf_wavelength:0.0, $
         tun_wavelength:0.0D, $
         exposure:0.0D, $
         gain:0.0, $
         cam_settings:'', $
         is_wb:0B, $
         lc:0B $
       }
  
  pc = { CRISP2_POLCAL_STATE, $
         inherits CRISP2_STATE, $
         lp:0.0, $
         qw:0.0 $
       }
  
end
