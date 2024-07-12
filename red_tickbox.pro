; docformat = 'rst'

;+
; Put a box with tickmarks around an image.
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
; :Returns:
; 
; 
; :Params:
; 
;   im : in, type=image
; 
;      Input image.
; 
;   sc : in, type=float
;
;      Scaling of tickmarks in arcsec (or Mm) per pixel 
;
; 
; :Keywords:
;
;   margin : in, optional, type=float
;
;     Margin around box in pixels, a small fraction of the image
;     dimensions. Negative margin gives inward pointing marks
;
;   arc5 : in, optional, type=boolean
;
;     Put longer ticks at intervals of 5.
;
;   arc10 : in, optional, type=boolean
;
;     Put even longer ticks at intervals of 10.
;
;   inv	: in, optional, type=boolean
;
;     Invert; white box with black ticks.
;
;   thickness	: in, optional, type=integet, default=1
;
;     Thickness of box and tickmarks, should be an odd number.
; 
; 
; :History:
; 
;   2024-07-12 : MGL, based on ANA code from 2000 by Luc Rouppe van
;                der Voort.
; 
;-
function red_tickbox, im, sc $
                      , margin = margin $
                      , arc5 = arc5 $
                      , arc10 = arc10 $
                      , inv = inv $
                      , rgb = rgb $
                      , thickness = thickness $
                      , tickcolor = tickcolor
  
  if not defined(thickness) then thickness = 1
  IF not defined(inv) THEN inv = 0
  if not defined(rgb) then rgb = 0
  
  if rgb then begin
     Ncolors = 3
     ;; Assume rgb is first dimension
     nx = dimen(im,1)	
     ny = dimen(im,2)
     imm = im
  endif else begin
     Ncolors = 1
     nx = dimen(im,0)	
     ny = dimen(im,1)
     if (size(imm))[(size(imm))[0]+1] eq 1 then begin ; Byte array
        imm = bytarr(1, nx, ny)
     endif else begin
        imm = fltarr(1, nx, ny)
     endelse
     imm[0, *, *] = im
  endelse
  
  if not defined(margin) then begin
     marg = round(max([nx,ny])/15.)
  end else begin
     marg = margin
  end
  
  if marg lt 0 then begin
     marg = -marg
     ticks_direction = -1
  end else begin
     ticks_direction =  1
  end
  
  marg = marg + odd(marg)
  
  if ticks_direction eq 1 then begin
     hmarg = marg/2
  end else begin
     hmarg = 0
  end
  
  if (size(imm))[(size(imm))[0]+1] eq 1 then begin ; Byte array
     if inv ne 0 then begin
        mi = byte(255)
        ma = byte(0)
     end else begin
        mi = byte(0)
        ma = byte(255)
     end
     box = bytarr(Ncolors, nx+2*hmarg+2*thickness, ny+2*hmarg+2*thickness) + mi
  end else begin
     if inv ne 0 then begin
        mi = max(imm)    
        ma = min([0,min(imm)])
     end else begin
        mi = min(imm)	
        ma = max(imm)
     endelse
     box = fltarr(Ncolors, nx+2*hmarg+2*thickness, ny+2*hmarg+2*thickness) + mi
  end
  
;  pix = fix(1/sc)
  pix = 1/sc
  lx = nx*sc	&	lx = fix(lx)
  ly = ny*sc	&	ly = fix(ly)

  if n_elements(tickcolor) eq 0 then tickcolor = ma
  
  for icolor = 0, Ncolors-1 do begin
    
;     print, 'Thickness:', thickness
    
    ;; Insert image in box
    box(icolor, hmarg+thickness, hmarg+thickness) = imm[icolor, *, *]
    
    ;; Draw box around image
;     help, box
;     print, hmarg, hmarg+nx+2*thickness-1, hmarg, hmarg+thickness-1
;     print, hmarg, hmarg+nx+2*thickness-1, hmarg+ny+1*thickness, hmarg+ny+2*thickness-1
;     print, hmarg, hmarg+thickness-1, hmarg, hmarg+ny+2*thickness-1
;     print, hmarg+nx+1*thickness, hmarg+nx+2*thickness-1, hmarg, hmarg+ny+2*thickness-1
    
    box(icolor, hmarg:hmarg+nx+2*thickness-1, hmarg:hmarg+thickness-1) = tickcolor ; Horizontal bottom
    box(icolor, hmarg:hmarg+nx+2*thickness-1, hmarg+ny+1*thickness:hmarg+ny+2*thickness-1) = tickcolor ; Horizontal top
    box(icolor, hmarg:hmarg+thickness-1, hmarg:hmarg+ny+2*thickness-1) = tickcolor                     ; Vertical left
    box(icolor, hmarg+nx+1*thickness:hmarg+nx+2*thickness-1, hmarg:hmarg+ny+2*thickness-1) = tickcolor ; Vertical right
    
    ;; Draw tick marks
    ;;for i = 0, lx do begin 
    for i = 1, lx do begin      ; Don't draw first tick mark!
;      fac = 0.4 * (2 + ((i mod 5) eq 0) + ((i mod 10) eq 0))/2.
      fac = 2
      if keyword_set(arc5)  and (i mod 5)  eq 0 then fac += 1
      if keyword_set(arc10) and (i mod 10) eq 0 then fac += 1
      fac *= 0.2
      
      dd = fix(marg/2 * fac * ticks_direction)
      if ticks_direction gt 0 then begin
        ;; Outward
        ;; Along horizontal bottom
          box(icolor, hmarg+i*pix, hmarg-dd:hmarg) = tickcolor
          ;; Along horizontal top
          box(icolor, hmarg+i*pix, hmarg+ny:hmarg+ny+dd) = tickcolor
         end else begin
           ;;print, 'Inward'
           ;;print, dd
           x = hmarg+thickness-1+i*pix
           ;; Along horizontal bottom
           box(icolor, x-thickness/2:x+thickness/2, hmarg+thickness-1:hmarg+thickness-1-dd) = tickcolor
           ;; Along horizontal top
           box(icolor, x-thickness/2:x+thickness/2, hmarg+ny+thickness-1+dd:hmarg+ny+thickness-1) = tickcolor
         end
     end
     for j = 1, ly do begin
;        fac = 0.4 * (2 + ((j mod 5) eq 0) + ((j mod 10) eq 0))/2.
       fac = 2
       if keyword_set(arc5)  and (j mod 5)  eq 0 then fac += 1
       if keyword_set(arc10) and (j mod 10) eq 0 then fac += 1
       fac *= 0.2

       dd = fix(marg/2 * fac * ticks_direction)
       if ticks_direction gt 0 then begin
         ;; Outward
          box(icolor, hmarg-dd:hmarg, hmarg+j*pix) = tickcolor
          box(icolor, hmarg+nx:hmarg+nx+dd, hmarg+j*pix) = tickcolor
         end else begin
           ;;print, 'Inward'
           y = hmarg+thickness-1+j*pix
           ;; Along vertical left
           box(icolor,  hmarg+thickness-1:hmarg+thickness-1-dd, y-thickness/2:y+thickness/2) = tickcolor
           ;; Along vertical right
           box(icolor, hmarg+nx+thickness-1+dd:hmarg+nx+thickness-1, y-thickness/2:y+thickness/2) = tickcolor
         end
     end
  endfor

  if rgb then begin
     return, box
  endif else begin
     ;; Remove extra dimension before returning
    return, reform(box)
  endelse

end
