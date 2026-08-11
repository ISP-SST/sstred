; docformat = 'rst'

;+
; Find the boundingbox for non-zero pixels in a mask.
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
; :Returns:
;
;   The boundingbox as [llx,lly,urx,ury]
; 
; :Params:
; 
;    mask : in, type="2d array" 
; 
;       A binary 2d mask.
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
;    2024-07-11 : MGL. First version.
; 
;-
function red_boundingbox, mask

  xproj = total(mask, 2)
  yproj = total(mask, 1)

  llx = min(where(xproj gt 0))
  urx = max(where(xproj gt 0))
  lly = min(where(yproj gt 0))
  ury = max(where(yproj gt 0))

  return, [llx, lly, urx, ury]

end

  
