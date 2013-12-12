; docformat = 'rst'

;+
; Class crispred 
;
;   2013-12-12 : move class definition to own file
; 
;-

PRO red__define
  done = {done, sumdark:0B, sumflat:0B, cleandata:0B, sumpolcal:0B, polcal:0B, makegains:0B}
                                ;
  struct = {red, $
            dark_dir:'',$
            flat_dir:ptr_new(),$
            data_dir:'',$
            data_list:strarr(100),$
            ndir:0B, $
            pinh_dir:'',$
            prefilter_dir:'',$
            polcal_dir:'',$
            camt:'', $
            camr:'', $
            camwb:'', $
            out_dir:'',$
            filename:'',$
            dopolcal:0B, $
            dodata:0B, $
            doflat:0B, $
            dodark:0B, $
            dopinh:0B, $
            docamt:0B, $
            docamr:0B, $
            docamwb:0B,$
            camttag:'',$
            camrtag:'',$
            camwbtag:'',$
            descatter_dir:'',$
            root_dir:'',$
            dodescatter:0B,$
            telescope_d:'',$
            image_scale:'',$
            pixel_size:'',$
            camsz:'',$
            done:done $
           }
END
