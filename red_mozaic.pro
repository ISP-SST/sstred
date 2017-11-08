; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
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
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_mozaic, img, clip = clip
                                ;
  return, mozaic(img.patch.img, $
                 img.patch[*,0].xl, $
                 img.patch[*,0].xh, $
                 img.patch[0,*].yl, $
                 img.patch[0,*].yh, clip= clip)
                                ;
end
