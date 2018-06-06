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
;-
PRO red__define

  done = {done, sumdark:0B, sumflat:0B, cleandata:0B, sumpolcal:0B, polcal:0B, makegains:0B}
                                
  struct = {red, $
            dark_dir:ptr_new(),$    ; The directory where raw dark frames are stored
            flat_dir:ptr_new(),$    ; The directories where raw flat fields are stored
            data_dir:'',$           ; The directory where raw science data is stored
            data_list:ptr_new(),$   ; The directories where raw science data is stored
            ndir:0B, $              ;
            pinh_dir:'',$           ; The directory where raw pinhole array data is stored
            prefilter_dir:'',$      ;
            polcal_dir:'',$         ; The directory where raw polcal data is stored
            camt:'', $              ;
            camr:'', $              ;
            camwb:'', $             ;
            out_dir:'',$            ; The directory where all output is stored
            filename:'',$           ;
            filetype:'', $          ;
            dopolcal:0B, $          ;
            dodata:0B, $            ;
            doflat:0B, $            ;
            dodark:0B, $            ;
            dopinh:0B, $            ;
            docamt:0B, $            ;
            docamr:0B, $            ;
            docamwb:0B,$            ;
            camttag:'',$            ;
            camrtag:'',$            ;
            camwbtag:'',$           ;
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
