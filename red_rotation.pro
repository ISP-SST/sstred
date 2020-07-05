; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
;    The shifted and rotated image.
; 
; 
; :Params:
; 
;    img : in, type=array 
;   
;       The image to be shifted and rotated.
;   
;    angle : in, type=scalar
;   
;      Rotation angle in radians.
;   
;    sdx :  in, type=scalar
;   
;      Shift in the X direction.
;   
;    sdy :  in, type=scalar
;   
;      Shift in the Y direction.
;   
; 
; :Keywords:
;
;    full : in, optional, type=array
;
;      Use this keyword to make full FOV cubes without clipping of
;      corners, otherwise the rotated and shifted frame will be of the
;      original size. Should be array with the following elements:
;      [maxangle, mdx0, mdx1, mdy0, mdy1 [, angles]], where elements
;      [1:4] are the max shifts in all directions.
; 
;    linear : in, optional, type=boolean
;
;      Use bilinear interpolation if this is set. Otherwise bicubic
;      interpolation with cubic = -0.5 is used.
;
;    background : in, optional, type=scalar, default="median(img)"
;
;      The full-FOV array is padded with this value.
;
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-06-10 : MGL. Made using bicubic interpolation with -0.5 the
;                default and introduced a new flag for linear
;                interpolation. 
;
;   2014-05-05 : JdlCR. Error in the definition of dx, dy found. The
;                shifts should be applied inside the parenthesis.
;
;   2014-11-29 : JdlCR. Added support for fullframe cubes (aka,
;                despite rotation and shifts, the entire FOV is inside
;                the image)
; 
;   2019-03-19 : MGL. New keyword background.
; 
;   2020-06-22 : MGL. Use entire range of angles to calculate
;                dimensions of rotated image iff keyword full is given
;                with more than 5 elements. Also fix bug so shifts are
;                added before rotation for the same calculation.
; 
;-
function red_rotation, img, angle, sdx, sdy $
                       , background = background $
                       , full = full $
                       , linear = linear 

  if n_elements(sdx) eq 0 then sdx = 0.0
  if n_elements(sdy) eq 0 then sdy = 0.0

  if n_elements(background) eq 0 then background = median(img)
  
  dim = size(img, /dim)

  if n_elements(full) gt 5 then begin

    ;; If "full" is only five elements, use old method below, kept for
    ;; backward compatibility. Otherwise use the angles part when
    ;; calculating the size of the rotated array. Also add shifts
    ;; before rotation, this is not done correctly in the old method.

    angles = full[5:*]

    ;; Corner coordinates, origin at center, with shifts
    x0 = -dim[0]/2. + full[1]
    x1 =  dim[0]/2. + full[2]  
    y0 = -dim[1]/2. + full[3]
    y1 =  dim[1]/2. + full[4] 

    ;; Rotated corner coordinates
    x00 = x0*cos(angles) - y0*sin(angles)
    x01 = x0*cos(angles) - y1*sin(angles)
    x10 = x1*cos(angles) - y0*sin(angles)
    x11 = x1*cos(angles) - y1*sin(angles)
    
    y00 = x0*sin(angles) + y0*cos(angles)
    y01 = x0*sin(angles) + y1*cos(angles)
    y10 = x1*sin(angles) + y0*cos(angles)
    y11 = x1*sin(angles) + y1*cos(angles)

    ;; Extreme rotated corner coordinates
    xmin = min([x00, x01, x10, x11])
    xmax = max([x00, x01, x10, x11])
    ymin = min([y00, y01, y10, y11])
    ymax = max([y00, y01, y10, y11])
    
    ;; Extreme rotated and shifted corner coordinates, origin at [0,0]
    ;; of old grid
    xmin = (xmin + dim[0]/2.) <0
    xmax = (xmax + dim[0]/2.) >dim[0]
    ymin = (ymin + dim[1]/2.) <0
    ymax = (ymax + dim[1]/2.) >dim[1]
    
    dim1 = round([xmax - xmin + 1, ymax - ymin + 1])
    
    ima = fltarr(dim1) + background
    ima[-round(xmin), -round(ymin)] = img

    ;; Apply rotation and shifts
    if keyword_set(linear) then begin
      return, red_rotshift(ima, angle, sdx, sdy, background = background, /linear)
    endif else begin
      return, red_rotshift(ima, angle, sdx, sdy, background = background, cubic = -0.5)
    endelse
    
  endif

  ;; Old mechanisms
  
  xsi = dim[0] * 0.5
  ysi = dim[1] * 0.5

  ;; Get the index of each matrix element and create an array for the
  ;; output image.
  xgrid = findgen(dim[0]) # (fltarr(dim[1]) + 1.0)
  ygrid = (fltarr(dim[0]) + 1.0) # findgen(dim[1])
  if n_elements(full) eq 0 then begin

    ;; Make the size of the rotated image the same as the original image
    
    dx = cos(angle) * (xgrid - xsi - sdx) - sin(angle) * (ygrid - ysi-sdy) + xsi
    dy = sin(angle) * (xgrid - xsi - sdx) + cos(angle) * (ygrid - ysi-sdy) + ysi

    ima = img

  endif else begin

    ;; Make the size of the rotated image larger so that it all fits.
    
    
    
    ;; For small angles, this should work.
    
    dx = cos(full[0]) * (xgrid - xsi) - sin(full[0]) * (ygrid - ysi) + xsi
    dy = sin(full[0]) * (xgrid - xsi) + cos(full[0]) * (ygrid - ysi) + ysi

    dim1 = dim
    xmin = round(abs(min(dx         )) + abs(full[1]))
    xmax = round(abs(max(dx-dim[0]-1)) + abs(full[2]))
    ymin = round(abs(min(dy         )) + abs(full[3]))
    ymax = round(abs(max(dy-dim[1]-1)) + abs(full[4]))

    dim1[0] += xmin + xmax
    dim1[1] += ymin + ymax

    ;; Force an odd number to make things easier for flipthecube to
    ;; find a correct dimension.

    if((dim1[0]/2)*2 ne dim1[0]) then begin
      dim1[0] += 1
      xmax += 1
    endif 

    if((dim1[1]/2)*2 ne dim1[1]) then begin
      dim1[1] += 1
      ymax += 1
    endif
    
    
    xgrid = (findgen(dim1[0]) # (fltarr(dim1[1]) + 1.0)) * float(dim[0])/float(dim[0])
    ygrid = ((fltarr(dim1[0]) + 1.0) # findgen(dim1[1])) * float(dim[1])/float(dim[1])
    xgrid *= float(dim1[0]) / max(xgrid)
    ygrid *= float(dim1[1]) / max(ygrid)
    xsi = dim1[0] * 0.5
    ysi = dim1[1] * 0.5
    
    dx = cos(angle[0]) * (xgrid - xsi - sdx) - sin(angle[0]) * (ygrid - ysi - sdy) + xsi
    dy = sin(angle[0]) * (xgrid - xsi - sdx) + cos(angle[0]) * (ygrid - ysi - sdy) + ysi

    ima = fltarr(dim1) + background
    ima[xmin, ymin] = img

  endelse

  ;; Interpolation onto new grid
  if keyword_set(linear) then begin
    return, interpolate(ima, dx, dy, missing = background)
  endif else begin
    return, interpolate(ima, dx, dy, missing = background, cubic = -0.5)
  endelse
  
end
