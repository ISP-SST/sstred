; docformat = 'rst'

;+
; Class TRIPPEL.
;
; :Author:
; 
;    Mats Löfdahl, ISP
;
;
; :History:
;
;   2017-03-13 : MGL. First version.
;
;-
pro trippel__define
                                
    struct = { TRIPPEL, inherits RED, $
               trippel_dummy:0B $     ; temporary dummy, move TRIPPEL specific content here
             }

    nb = { TRIPPEL_STATE, $
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
                                
    pc = { TRIPPEL_POLCAL_STATE, $
           inherits TRIPPEL_STATE, $
           lp:'', $
           qw:'' $
    }
                                
end
