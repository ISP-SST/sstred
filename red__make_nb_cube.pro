; docformat = 'rst'

;+
; Wrapper method for make_nb_cube with or without demodulation.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;    nopolarimetry : in, optional, type=boolean
;
;       For a polarimetric dataset, don't make a Stokes cube.
;       Instead combine all LC states for both cameras into a single
;       NB image per tuning, producing a cube similar to that for a
;       data set without polarimetry. (For a nonpolarlimetric dataset,
;       no effect.)
;
  
;   
;   
; 
; 
; :History:
; 
; 
; 
;-
pro red::make_nb_cube, wcfile $
                       , nopolarimetry = nopolarimetry $
                       , _ref_extra = _extra

  instrument = ((typename(self)).tolower())
  polarimetric_data = self -> polarimetric_data()

  case 1 of

    ;; CHROMIS data with WDN cameras
    ~polarimetric_data : self -> make_nb_cube_intensity, wcfile, _extra = _extra

    ;; CRISP or CHROMIS with W(D)TR cameras, but we don't want Stokes
    keyword_set(nopolarimetry) : self -> make_nb_cube_intensity, wcfile, _extra = _extra

    ;; CRISP or CHROMIS with W(D)TR cameras
    else : self -> make_nb_cube_stokes, wcfile, _extra = _extra

  endcase
  
end
