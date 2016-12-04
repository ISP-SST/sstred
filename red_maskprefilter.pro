function red_maskprefilter, wav, imean

  !p.multi = [0,1,1]
  window, 2, xsize=1024, ysize=500
                                ;
  plot, wav, imean, psym = -1, xstyle = 3, ystyle = 3, $
        xtitle = 'Wavelength [A]', ytitle = 'Intensity', $
        title = 'Mask out regions from the fit (e.g., line core). Right mouse-click exits.'
                                ;
                                ; Loop

  w = dblarr(n_elements(wav)) + 1.0d0
  
  !MOUSE.button = 0
  fclick = 1B
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
           w[pos] = 0.0d0
           loadct,3, /silent
           oplot, wav[pos], imean[pos], color = 160, psym = -1
           loadct,0, /silent   
        endif
        fclick = 1B
     endelse
  endwhile
  !MOUSE.button = 0
  wdelete, 2
  return, w
end
