; docformat = 'rst'

;+
; Class crispred 
;
; :History:
;
;   2013-12-12 : Move class definition to own file
;
;   2014-01-08 : MGL. Added members: isodate, log_dir, telog, pinhole_spacing.
; 
;   2014-01-09 : MGL. Added some documentation.
; 
;   2014-01-10 : PS  New config variable filtype
;
;   2014-03-21 : MGL. Allow for multiple dark_dir.
;
;   2016-05-19 : THI. Define state structures
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;-
PRO red__define

    st = { RED_STATE, $
           filename:'', $
           detector:'', $
           camera:'', $
           fullstate:'', $
           skip:0B $
    }
                                
  done = {done, sumdark:0B, sumflat:0B, cleandata:0B, sumpolcal:0B, polcal:0B, makegains:0B}
                                
  struct = {red, $
            dark_dir:ptr_new(),$    ; The directories where raw dark frames are stored
            flat_dir:ptr_new(),$    ; The directories where raw flat fields are stored
            data_dirs:ptr_new(),$   ; The directories where raw science data is stored
            current_data_dir:0B,$   ; Index of data directory to process (default is 0)
            pinh_dirs:ptr_new(),$   ; The directories where raw pinhole array data is stored
            prefilter_dir:'',$      ;
            polcal_dir:'',$         ; The directory where raw polcal data is stored
            detectors:ptr_new(),$   ; List of detector identifiers (e.g. camXX)
            cameras:ptr_new(),$     ; List of camera locations/channels (e.g. Crisp-R)
            refcam:0B, $            ; Selected reference channel
            out_dir:'',$            ; The directory where all output is stored
            filename:'',$           ;
            filetype:'', $          ;
            dopolcal:0B, $          ;
            dodata:0B, $            ;
            doflat:0B, $            ;
            dodark:0B, $            ;
            dopinh:0B, $            ;
            descatter_dir:'',$      ;
            root_dir:'',$           ;
            dodescatter:0B,$        ;
            telescope_d:'',$        ; The diameter of the telescope pupil in meters
            image_scale:'',$        ; The image scale in arxsec per pixel
            pixel_size:'',$         ; The pixel size in meters on the detector
            camsz:'',$              ;
            done:done, $            ;
            isodate:'', $           ; The date of observations in ISO format YYYY-MM-DD
            log_dir:'', $           ; Path to directory where log files are downloaded
            telog:'' , $            ; Path to telescope pointing log
            pinhole_spacing:0.0 $   ; Pinhole array spacing in arcsec
           }
END
