function red_select_area, im, noedge=noedge

  
  im1 = fix(im)*0
  if(keyword_set(noedge)) then im1[*,*] = 1 else im1[50:-50,50:-50] = 1

  imo = red_histo_opt(im*im1,2.e-3)
  mma = max(imo)
  mmi = min(imo)
  
  dim = size(im, /dim)

  
  window, 2, xsize=1280*1.1, ysize=(1280*1.1*dim[1])/dim[0], title='Mask out regions with active patches, right button exits!'
  red_tvimg, bytscl(imo, mmi, mma), /sam, /nosc
  
  fclick = 1B
  
  !MOUSE.button = 0
  while(!MOUSE.button NE 4) do begin
     cursor, t1, t2, /down, /data
     t1 = (t1>0)<(dim[0]-1)
     t2 = (t2>0)<(dim[1]-1)
     
     if(fclick) then begin
        tt1 = t1
        tt2 = t2
        fclick = 0B
     endif else begin
        im1[tt1:t1, tt2:t2] = 0
        im2 = im
        idx = where(im1 eq 0, count)
        if(count gt 0) then im2[idx] *= 0.35
        red_tvimg, bytscl(im2, mmi, mma), /sam, /nosc
        fclick = 1B
     endelse
  endwhile
  
  wdelete, 2
  
  return, im1
end
