; docformat = 'rst'

;+
; Clip and optionally mirror an image the same way as the momfbd
; program would, using its align_clip info.
; 
; :categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    Pit Sütterlin (as momfbd_clip.pro)
; 
; 
; :returns:
; 
;    A clipped and optionally mirrored image.
; 
; :params:
; 
;    pic : in, type="2D array"
;   
;      The image to be clipped.
;   
;    cl : in, type="integer array(4)"
;   
;      The align_clip info. In its simplest form, just
;      [llx,urx,lly,ury] where llx<urx and lly<ury, specifying that
;      the returned image should be pic[llx-1:urx-1,lly-1:ury-1].
;      (Here, "ll" means "lower left" and "ur" means "upper right".)
;
;      But if the order of llx and urx is switched, it means the image
;      should be mirrored in the X direction. Similar for lly, ury and
;      the Y direction.
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-06-10 : Added some documentation. Mats Löfdahl (MGL). 
;
;-
function red_clipim, pic, cl

  if cl[0] GT cl[1] then begin
     ;; Should mirror in X.
     mi_x = cl[1]-1
     ma_x = cl[0]-1
     fx = 1
  endif else begin
     ;; Should not mirror in X.
     mi_x = cl[0]-1
     ma_x = cl[1]-1
     fx = 0
  endelse

  if cl[2] GT cl[3] then begin
     ;; Should mirror in Y.
     mi_y = cl[3]-1
     ma_y = cl[2]-1
     fy = 1
  endif else begin
     ;; Should not mirror in Y.
     mi_y = cl[2]-1
     ma_y = cl[3]-1
     fy = 0
  endelse

  ;; Translate the mirror settings fx and fy to the orientations
  ;; recogniced by IDL's rotate command.
  case fx+2*fy of 
     0:  r = 0
     1:  r = 5
     2:  r = 7
     3:  r = 2
  endcase

  ;; Return the clipped image with the specified orientation.
  return, rotate(pic[mi_x:ma_x, mi_y:ma_y], r)

end

