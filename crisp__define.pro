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

    nb = { CRISP_STATE, $
           inherits RED_STATE, $
           scannumber:0, $
           framenumber:-1L, $
           tuning:'', $
           prefilter:'', $
           pf_wavelength:0.0, $
           tun_wavelength:0.0D, $
           lc:'' $
    }
                                
    pc = { CRISP_POLCAL_STATE, $
           inherits CRISP_STATE, $
           lp:'', $
           qw:'' $
    }
                                
end
