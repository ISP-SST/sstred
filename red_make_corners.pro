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
;    Tomas Hillberg, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
;
;   2016-06-13 : MGL. Renamed from make_corners.
; 
; 
;-
function red_make_corners, clip
    dim = size(clip,/dim)
    corners = intarr(4,3)
    if dim(0) gt 3 then begin
        corners(*,2) = 1
        corners([0,2], 0) = clip(0)
        corners([1,3], 0) = clip(1)
        corners([0,1], 1) = clip(2)
        corners([2,3], 1) = clip(3)
    endif
    return, corners
end

