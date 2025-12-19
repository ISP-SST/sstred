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
;    inverse : in, type=boolean
;   
;      Do the inverse operation: the input im is embedded in an
;      original size array, at the appropriate position.
; 
; 
; :History:
; 
;    2022-09-12 : MGL. First version.
; 
;    2025-12-19 : MGL. New keyword inverse.
; 
;-
function red_crop_as_momfbd, im, momfbd_struct, inverse = inverse

  ;; May have to take the version into account for older momfbd files.
  ;; In case the coordinates have to be transposed.
  
  x0 = momfbd_struct.roi[0] + momfbd_struct.margin
  x1 = momfbd_struct.roi[1] - momfbd_struct.margin
  y0 = momfbd_struct.roi[2] + momfbd_struct.margin
  y1 = momfbd_struct.roi[3] - momfbd_struct.margin

  if keyword_set(inverse) then begin

    ;; Make a larger array and put im in it at the appropriate
    ;; position.

    ;; Get the original dimensions
    hdr = red_readhead(momfbd_struct.name[0])
    naxis = fxpar(hdr, 'NAXIS*')
    dims = naxis[0:1]

    im_out = make_array(dims, type = size(im, /type))
    im_out[x0:x1,y0:y1] = im

    return, im_out
    
  endif else begin
    
    return, im[x0:x1,y0:y1]

  endelse
  
end
