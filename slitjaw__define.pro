; docformat = 'rst'

;+
; Class SLITJAW.
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
pro slitjaw__define
                                
    struct = { SLITJAW, inherits RED, $
               slitjaw_dummy:0B $     ; temporary dummy, move SLITJAW specific content here
             }

    nb = { SLITJAW_STATE, $
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
                                
    pc = { SLITJAW_POLCAL_STATE, $
           inherits SLITJAW_STATE, $
           lp:'', $
           qw:'' $
    }
                                
end
