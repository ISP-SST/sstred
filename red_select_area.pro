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
;    Jaime de la Cruz Rodriguez, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
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
;    2017-04-07 : MGL. New keyword xroi, use XROI GUI to select area.
; 
; 
;-
function red_select_area, im, noedge=noedge, xroi = xroi

  margin = 50
  
  if keyword_set(xroi) then begin
    
    imshow = im
    dims = size(im, /dim)

    mask = bytarr(dims) 
    if(keyword_set(noedge)) then begin
      mask[margin:-margin,margin:-margin] = 1
      imshow[margin:-margin,margin:-margin] = 1.5*imshow[margin:-margin,margin:-margin]
    endif else begin
      mask[*, *] = 1
    endelse

    print
    print, 'Use the XROI GUI to select areas that are NOT to be used.'
    if(keyword_set(noedge)) then $
       print, 'The outermost '+strtrim(margin, 2)+' pixels are automatically deselected (as indicated).'

    
    mn = median(im)-2.5*robust_sigma(im)
    mx = median(im)+2.5*robust_sigma(im)
    XROI, bytscl(imshow, mn, 1.5*mx), Regions_Out=ROIs, /Block, title = 'Deselect areas.'

    if ROIs ne !null then $
       for iroi = 0, n_elements(ROIs)-1 do $
          mask = mask and (ROIs[iroi] -> ComputeMask(Dimensions=dims, Mask_Rule=2) eq 0)

    Obj_Destroy, ROIs

    scrollwindow, xs = dims[0], ys = dims[1]
    tv, bytscl(mask*im, mn, mx)

    return, mask

  endif else begin

    im1 = fix(im)*0
    if(keyword_set(noedge)) then im1[*,*] = 1 else im1[margin:-margin,margin:-margin] = 1

    imo = red_histo_opt(im*im1,2.e-3)
    mma = max(imo)
    mmi = min(imo)
    
    dim = size(im, /dim)

    msg = 'Mask out regions with active patches, right button exits!'
    window, 2, xsize=1280*1.1, ysize=(1280*1.1*dim[1])/dim[0], title=msg
    red_tvimg, bytscl(imo, mmi, mma), /sam, /nosc
    
    print
    print, msg

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

  endelse
  

end
