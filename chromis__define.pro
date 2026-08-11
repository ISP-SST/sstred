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
;                independent parts) and derived classes
;                (CRISP/CHROMIS). 
;
;   2016-05-24 : MGL. Removed LC from CHROMIS_STATE, added gain and
;                exposure.  
;
;   2016-05-25 : MGL. Make the gain in CHROMIS_STATE a float.
;
;-
pro chromis__define
                 
  struct = { CHROMIS, inherits RED, $
             chromis_dummy:0B $ ; temporary dummy, move CHROMIS specific content here
           }

  nb = { CHROMIS_STATE, $
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
  
  pc = { CHROMIS_POLCAL_STATE, $
         inherits CHROMIS_STATE, $
         lp:0.0, $
         qw:0.0 $
       }
  
end
