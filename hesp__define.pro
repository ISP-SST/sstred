; docformat = 'rst'

;+
; Class HESP
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
;   2023-11-21 : MGL. New version for HeSP based on the chromis
;                version.
;
;-
pro hesp__define
  
  struct = { HESP, inherits RED, $
             hesp_dummy:0B $    ; temporary dummy, move CHROMIS specific content here
           }

  nb = { HESP_STATE, $
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
         is_wb:0B $
       }
  
  pc = { HESP_POLCAL_STATE, $
         inherits HESP_STATE, $
         lp:'', $
         qw:'' $
       }
  
end
