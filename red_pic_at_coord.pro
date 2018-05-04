; docformat = 'rst'

;+
; Return an xSize by ySize version of Pic, that is centered on x,y
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
; 
; :Params:
; 
;   pic : array
; 
;      The array from which the subfield is to be extracted.
; 
;   x : in, type=integer
; 
;      The center x coordinate of the subfield.
; 
;   y : in, type=integer
; 
;      The center y coordinate of the subfield. 
; 
;   xSize :   in, type=integer
; 
;      The x size of the subfield.
; 
;   ySize : in, type=integer
; 
;      The y size of the subfield. 
; 
; 
; :Keywords:
; 
;    mode : in, optional, type=boolean
;   
;      Set to force correct, even size.
; 
; 
; :History:
; 
;   1998-12-03 : MGL. First ANA version as pic_at_coord.
; 
;   2014-02-27 : MGL. Allow subfields that are not completely within
;                the original image.
; 
;   2018-05-04 : MGL. Add to pipeline repository as red_pic_at_coord.
;
;-
function red_pic_at_coord, pic, x, y, xSize, ySize, mode = mode 
  
  llx = x - xSize/2 
  lly = y - ySize/2 
  urx = x + xSize/2 - 1 
  ury = y + ySize/2 - 1
  
  IF urx - llx + 1 LT xSize THEN urx++
  IF ury - lly + 1 LT ySize THEN ury++

  dims = size(pic, /dim)
  
  llx = llx > 0	   
  lly = lly > 0	   
  urx = urx < (dims[0]-1)
  ury = ury < (dims[1]-1)
  
  if keyword_set(mode) then begin
    if odd(urx-llx+1) then urx--
    if odd(ury-lly+1) then ury--
  end
  
  ppic = pic(llx:urx,lly:ury)
  
;  IF mode THEN BEGIN
;    box = FltArr(xSize,ySize)+mean(ppic)
;    ;    help,ppic,box
;    insert,box,ppic,llx-(x-xSize/2),lly-(y-ySize/2)
;    
;    Return,box
;  END ELSE BEGIN
  return,ppic 
;  END
  
END                             ; Pic_at_coord

