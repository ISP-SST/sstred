; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    wav : 
;   
;   
;   
;    imean : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_select_spoints, wav, imean
  !p.multi = [0,1,1]
                                ;
  cgplot, wav, imean, psym = -16, xstyle = 3, ystyle = 3, $
          xtitle = 'Wavelength [A]', ytitle = 'Stokes I', $
          title = 'Select points regions to calculate cross-talk from I'
                                ;
                                ; Loop
                                ;
  !MOUSE.button = 0
  fclick = 1B
  idx = [-1L]
  while(!MOUSE.button NE 4) do begin
     cursor, t1, t2, /down, /data
     if(fclick) then begin
        tw = t1
        fclick = 0B
     endif else begin
        range = [tw, t1]
        range = range[sort(range)]
        pos = where((wav LE range[1]) AND (wav GE range[0]), count)
        if(count GE 1) then begin
           idx = [temporary(idx), pos]
;           loadct,3, /silent
           cgplot, /over, wav[pos], imean[pos], color = 'red', psym = -16
;           loadct,0, /silent   
         endif
        fclick = 1B
     endelse
  endwhile
  !MOUSE.button = 0
                                ;
  return, idx[1:*]
end
