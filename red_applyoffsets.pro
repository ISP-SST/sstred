; docformat = 'rst'

;+
; Apply xoffs and yoffs corrections to an image.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
;
;    Jaime de la Cruz Rodriguez (ISP-SU 2014)
; 
; 
; :Returns:
; 
;     A version of the input image where the offsets (from pinhole
;     calibrations) have been applied.
; 
; 
; 
; 
; :Params:
; 
;    img : in, type="2D array"
; 
;       The image.
; 
; 
;    xoff : in, , type="2D array"
; 
;       The X offsets.
; 
;    yoff : in, , type="2D array"
; 
;       The Y offsets.
; 
; 
; 
; 
; 
; :Keywords:
; 
; 
; :History:
;
;    2014-10-16 : MGL. Proper comment section. Added optional
;                 align-clipping. 
;
;-
function red_applyoffsets, img, xoff, yoff, clips = clips
  
  ;; Get function name
  inam = red_subprogram(/low, calling = inam1)                                    

  if n_elements(clips) ne 0 then imc = red_clipim(img, clips) else imc = img

  ;; compare dimensions
  dim = size(imc,/dim)
  di  = size(xoff, /dim)
  
  if max(dim-di) ne 0 then begin
     print, inam +'Error, the image and xoffs and yoffs grids have different dimensions:'
     print, '   imc -> ', string(dim[0], format='(I0)')+'x'+string(dim[1], format='(I0)')
     print, '   {x,y}offs -> ', string(di[0], format='(I0)')+'x'+string(di[1], format='(I0)')
     print, inam + 'Corrections NOT APPLIED!'
     return, img
  endif

  ;; create interpolation grids (offsets given in 1/100 of a pixel)
  xx = dindgen(dim[0]) # replicate(1.0d0, dim[1]) + 1.0d-2 * xoff
  yy = replicate(1.0d0, dim[0]) # dindgen(dim[1]) + 1.0d-2 * yoff
  
  ;; interpolate data
  return, float(interpolate(temporary(imc), xx, yy, cubic=-0.5))

end
