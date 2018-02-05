; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    img : 
;   
;   
;   
; 
; :Keywords:
; 
;    clip  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;   2018-02-05 : THI: pass all unrecognized keywords to mozaic.
; 
;-
function red_mozaic, momfbd_struct, _REF_EXTRA = ex
  
  return, mozaic( momfbd_struct.patch.img, $
                  momfbd_struct.patch[*,0].xl, $
                  momfbd_struct.patch[*,0].xh, $
                  momfbd_struct.patch[0,*].yl, $
                  momfbd_struct.patch[0,*].yh, $
                  _EXTRA=ex)
  
end
