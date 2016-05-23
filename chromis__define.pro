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
pro chromis__define
                                
    struct = { CHROMIS, inherits RED, $
               chromis_dummy:0B $     ; temporary dummy, move CHROMIS specific content here
             }

    nb = { CHROMIS_STATE, $
           inherits RED_STATE, $
           scannumber:0, $
           framenumber:-1L, $
           tuning:'', $
           prefilter:'', $
           pf_wavelength:0.0, $
           tun_wavelength:0.0D, $
           lc:'' $
    }
                                
    pc = { CHROMIS_POLCAL_STATE, $
           inherits CHROMIS_STATE, $
           lp:'', $
           qw:'' $
    }
                                
end
