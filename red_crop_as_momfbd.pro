; docformat = 'rst'

;+
; Crop a raw-sized image as the momfbd-restored image returned by
; rdx_mozaic(momfbd_struct,/crop).
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
;    The cropped image.
;
; :Params:
; 
;    im : in, type="array(Nx,Ny)"
; 
;      The image to be cropped.
; 
;    momfbd_struct : in, type=struct
;
;      A struct as returned by momfbd_read(), in particular including
;      momfbd_struct.roi and momfbd_struct.margin.
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
;    2022-09-12 : MGL. First version.
; 
;-
function red_crop_as_momfbd, im, momfbd_struct

  x0 = momfbd_struct.roi[0] + momfbd_struct.margin
  x1 = momfbd_struct.roi[1] - momfbd_struct.margin
  y0 = momfbd_struct.roi[2] + momfbd_struct.margin
  y1 = momfbd_struct.roi[3] - momfbd_struct.margin
  
  ;; May have to take the version into account for older momfbd files.
  ;; In case the coordinates have to be transposed.
  
  return, im[x0:x1,y0:y1]

end
