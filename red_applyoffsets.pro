;; function red_applyoffsets
;; Purpose: apply xoffs and yoffs corrections to an image.
;; Author: Jaime de la Cruz Rodriguez (ISP-SU 2014)
;;
function red_applyoffsets, img, xoff, yoff
  
  ;; Get function name
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])+' : '


  ;; compare dimensions
  dim = size(img,/dim)
  di    = size(xoff, /dim)
  if max(dim-di) ne 0 then begin
     print, inam +'Error, the image and xoffs and yoffs grids have different dimensions:'
     print, '   img -> ', string(dim[0], format='(I0)')+'x'+string(dim[1], format='(I0)')
     print, '   {x,y}offs -> ', string(di[0], format='(I0)')+'x'+string(di[1], format='(I0)')
     print, inam + 'Corrections NOT APPLIED!'
     return, img
  endif

  ;; create interpolation grids (offsets given in 1/100 of a pixel)
  xx = dindgen(dim[0]) # replicate(1.0d0, dim[1]) + 1.0d-2 * xoff
  yy = replicate(1.0d0, dim[0]) # dindgen(dim[1]) + 1.0d-2 * yoff
  
  ;; interpolate data
  return, float(interpolate(temporary(img), xx, yy, cubic=-0.5))
end
