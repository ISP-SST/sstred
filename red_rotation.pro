; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author: J. de la Cruz Rodriguez & M. Löfdahl, Institute for Solar Physics
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
;    stretch_grid : in, optional, type=array
;
;       stretch grid computed with red_gridmatch (2, nTiles_x, nTiles_y)
;
;    nthreads : in, optional, type=int
;   
;       Number of threads to use during the interpolation
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
;   2020-10-01 : JdlCR. Modified this routine to a) only allow
;                bilinear o nearest neighbor interpolation and
;                b) to allow including the stretch correction
;                in this interpolation so all corrections are 
;                applied at once. Use our own interpolation routines.
;
; 
;-
function red_rotation, img, angle, sdx, sdy $
                       , background = background $
                       , full = full $
                       , stretch_grid = stretch_grid $
                       , nthreads = nthreads $
                       , nearest = nearest $
                       , stretch_matrix = stretch_matrix

  if n_elements(sdx) eq 0 then sdx = 0.0
  if n_elements(sdy) eq 0 then sdy = 0.0

  if n_elements(background) eq 0 then background = median(img)
  
  dim = size(img, /dim)

  if n_elements(full) ge 5 then begin

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
    
    ima = dblarr(dim1) + background
    ima[-round(xmin), -round(ymin)] = img
    

    return, red_rotshift(ima, angle, sdx, sdy, background = background, stretch_grid = stretch_grid $
                         , nthreads=nthreads, nearest = nearest, stretch_matrix = stretch_matrix, original_dimensions = dim)
    
  endif else begin

    ;; Old mechanisms
    
    ;; xsi = dim[0] * 0.5
    ;; ysi = dim[1] * 0.5
    
    ;; ;; Get the index of each matrix element and create an array for the
    ;; ;; output image.
    ;; xgrid = dindgen(dim[0]) # (dblarr(dim[1]) + 1.0d)
    ;; ygrid = (dblarr(dim[0]) + 1.0d) # dindgen(dim[1])
    
    
    ;; ;; Add destretch grid correction
    
    ;; if(n_elements(stretch) ne 0) then begin
    ;;   smat = red_get_full_stretch_matrix(dim[0], dim[1], stretch, original_size=dim, /only_shifts)
    ;;   smatx = smat[*,*,0]; - xgrid 
    ;;   smaty = smat[*,*,1]; - ygrid
    ;; endif else begin
    ;;   smatx = xgrid*0
    ;;   smaty = ygrid*0
    ;; endelse
    ;; ;; Make the size of the rotated image the same as the original image

    ;; dx = cos(angle) * (xgrid - xsi - sdx + smatx) - sin(angle) * (ygrid - ysi - sdy + smaty) + xsi 
    ;; dy = sin(angle) * (xgrid - xsi - sdx + smatx) + cos(angle) * (ygrid - ysi - sdy + smaty) + ysi
    ;; ;;ima = img

    ;; return, red_interpolate2D(xgrid[*,0], reform(ygrid[0,*]), img, dx, dy, nthreads=nthreads, nearest = nearest)
    
    return, red_rotshift(img, angle, sdx, sdy, background = background, stretch_grid = stretch_grid $
                         , nthreads=nthreads, nearest = nearest, stretch_matrix = stretch_matrix, original_dimensions = dim)
    
  endelse
end
