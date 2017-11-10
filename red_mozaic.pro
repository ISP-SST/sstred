; docformat = 'rst'

;+
; Make a mosaic image from the subfields in a momfbd struct.
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
;   The mosaicked image.
; 
; 
; :Params:
; 
;    momfbd_struct : in, type=struct
;   
;      A momfbd struct as returned by momfbd_read(filename) or
;      momfbd_read(filename,/img). 
;   
; 
; :Keywords:
; 
;    clip : in, optional, type=boolean
;   
;       Crop the mosaicked image to remove the default border as well
;       as any lines and columns that are all zero.
;   
;    crop : in, optional, type=boolean
; 
;       Crop the mosaicked image to remove the default border.
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2017-11-08 : THI: new keyword crop passed to mozaic, which will crop
;                away the default border around the mozaic.
; 
;-
function red_mozaic, momfbd_struct, clip = clip, crop = crop
                                ;
  return, mozaic(momfbd_struct.patch.img, $
                 momfbd_struct.patch[*,0].xl, $
                 momfbd_struct.patch[*,0].xh, $
                 momfbd_struct.patch[0,*].yl, $
                 momfbd_struct.patch[0,*].yh, $
                 clip=clip, crop=crop)
                                ;
end
