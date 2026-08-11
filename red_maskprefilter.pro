; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; 
; 
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
;    2017-05-04 : MGL. Give more visual feedback when mouse-clicking
;                 in the plot. 
; 
; 
; 
; 
;-
function red_maskprefilter, wav, imean, img = img, mask = mask

  window, 2, xsize=1024, ysize=500
                                ;
  cgplot, wav, imean, psym = -1, xstyle = 3, ystyle = 3, $
          xtitle = 'Wavelength [A]', ytitle = 'Intensity', $
          title = 'Mask out data points from the fit (e.g., line core).'

  cgtext, /normal, 0.5, 0.85, align = 0.5 $
          , 'Left mouse click once to select beginning of mask, once more to select end of mask.'
  cgtext, /normal, 0.5, 0.80, align = 0.5 $
          , 'Right click exits.'

  w = dblarr(n_elements(wav)) + 1.0d0
  
  !MOUSE.button = 0
  fclick = 1B
  while(!MOUSE.button NE 4) do begin
    cursor, t1, t2, /down, /data
    if(fclick) then begin
      tw = t1
      range = [tw, max(wav)+1]
      range = range[sort(range)]
      pos = where((wav LE range[1]) AND (wav GE range[0]), count)
      if(count GE 1) then begin
        cgplot, /over, wav, imean, color = 'blue', psym = -1, thick = 2
        cgplot, /over, wav[pos], imean[pos], color = 'red', psym = -1, thick = 2
      endif
      fclick = 0B
    endif else begin
      range = [tw, t1]
      range = range[sort(range)]
      pos = where((wav LE range[1]) AND (wav GE range[0]), count)
      if(count GE 1) then begin
        w[pos] = 0.0d0
        cgplot, /over, wav, imean, color = 'blue', psym = -1, thick = 2
        cgplot, /over, wav[pos], imean[pos], color = 'red', psym = -1, thick = 2
      endif
      fclick = 1B
    endelse
  endwhile
  !MOUSE.button = 0
  wdelete, 2

  ;; Now the FOV
  if(n_elements(img) gt 0) then begin
    tmp = red_select_area(img)
    mask = where(tmp gt 0, count)
  endif
   
  return, w

end
