; docformat = 'rst'

;+
; Use XROI GUI to select areas that are NOT to be used.
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
;    2021-05-01 : MGL. Missing data, already deselected pixels in yellow.
; 
; 
;-
function red_select_area, im_in, noedge=noedge, xroi = xroi

  im = im_in

  if keyword_set(xroi) then begin
    
    imshow = im
    dims = size(im, /dim)

    mask = finite(im)
    indx_data = where(mask, Ndata, complement = indx_missing, ncomplement = Nmissing)

    print
    print, 'Use the XROI GUI to select areas that are NOT to be used.'
    if(keyword_set(noedge)) then $
       print, 'Yellow pixels are already deselected.'
    print
    
;    indx_missing = where(~finite(im), Nwhere, complement = indx_data)
    
    mn = median(im[indx_data])-2.5*robust_sigma(im[indx_data])
    mx = median(im[indx_data])+2.5*robust_sigma(im[indx_data])
;    imshow = bytscl(imshow, mn, 1.5*mx)
    imshow = bytscl(imshow, mn, mx)
    rg = imshow
    b = imshow

    if Nmissing gt 0 then begin
      ;; Yellow for missing data
      rg[indx_missing] = 255
      b[indx_missing] = 0
    endif
    XROI, [[[rg]],[[rg]],[[b]]], Regions_Out=ROIs, /Block, title = 'Deselect areas.'

    if ROIs ne !null then $
       for iroi = 0, n_elements(ROIs)-1 do $
          mask = mask and (ROIs[iroi] -> ComputeMask(Dimensions=dims, Mask_Rule=2) eq 0)

    Obj_Destroy, ROIs

    
    scrollwindow, xs = dims[0], ys = dims[1]
    tv, bytscl(mask*im, mn, mx)
    im2 = im
    indx_deselected = where(~mask, Ndeselect)
    if Ndeselect gt 0 then im2[indx_deselected] = !values.f_nan
    cgimage, im2 $
             , missing_color = 'yellow' $
             , missing_index = 0 $
             , stretch = 1  
    
    return, mask

  endif else begin

    margin = 15
    
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
