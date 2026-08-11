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
;    x : 
;   
;   
;   
;    y : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_getclips, img, x, y
                                ;
  dx = -img.patch[x,y].offx
  dy = -img.patch[x,y].offy
                                ;
  xoff = img.clip[0,0,0]
  yoff = img.clip[0,1,0]

  clip = [xoff + img.patch[x,y].xl + img.patch[x,y].dx - 2,$
          xoff + img.patch[x,y].xh + img.patch[x,y].dx - 2,$
          yoff + img.patch[x,y].yl + img.patch[x,y].dy - 2,$
          yoff + img.patch[x,y].yh + img.patch[x,y].dy - 2]
                                ;
  clip[0:1] += dx
  clip[2:3] += dy

  return, clip
end
