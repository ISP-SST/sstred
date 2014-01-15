FUNCTION Red_remap_image, ref, pic, cs, sw, ALIGN=al, MAP=d_map, PLOT=plot, QUINT=quint, $
                          RANGE = range, SAME_SCREEN = same, SMOOTH = sm, MEDIAN = med, $
                          DONT_FIX_BORDER = dfb, REFINE_SHIFTS = rs, rlim = rlim, RCS = rcs
;+
; NAME:
;       RED_REMAP_IMAGE
; PURPOSE:
;       Do a 1:1 alignment of two images using warping techniques
; CALLING SEQUENCE:
;       Res = Map_Image( Ref, Image, Cell [, Step ] [, /ALIGN]
; INPUTS:
;       Ref:    Reference image
;       Image:  Original image
;       Cell:   (integer) 
; OPTIONAL PARAMETERS:
;       Step:   (int) Stepwidth of the correlation grid.  If omitted,
;               same value as Cell is used
; KEYWORDS:
;       ALIGN   (Flag) If set, an allover crosscorrelation of the
;               images will be performed to remove a full-frame shift.
;
;       MAP:    name of a variable to pass back the distortion map.
;
;       PLOT:   (float) Display image and shift vectors, latter scaled
;               by a factor of PLOT 
;
;       RANGE:  (integer) Allowed range for displacements.
;               Correlations will only be looked for in +/-RANGE.
;               The parameter is just passed to shc.
;       REFINE_SHIFTS: (flag) If the computed shift is larger than RLIM then
;               recompute it with the properly centered subimage
;       RLIM:   (integer) Limit when to recompute shifts (DEFAULT: 1)
;       RCS:    (integer) Cell size for shift refinement. (DEFAULT: cell)
;       SAME_SCREEN: (flag) Use the currently active window for
;               plotting.  Default is to open a new one.
;       SMOOTH: (integer) Smooth images with a boxcar of width SMOOTH
;               before computing the displacement map.
; OUTPUTS:
;       Res:    Warped image of same type as input image, size is the
;               min of the sizes of the two images.
; RESTRICTIONS:
;       
; PROCEDURE:
;       Divide the two images into subfields of cellsize Cell, using
;       stepwidth Step.  Compute the shifts between the subimages and
;       set-up a polynomial warping matrix.  Call WARP_TRI to do the
;       warping. 
; MODIFICATION HISTORY:
;       10-Jun-1997  P.Suetterlin, USG
;       11-Sep-2001  call empty after plot to assure all arrows are there
;       05-Feb-2002  Rename to remap_image - there's already a
;                    map_image in the standard IDL.
;       04-Mar-2002  Smoothing of images before calculation of shifts
;       13-Jan-2014  Option to refine (large) shifts
;       15-Jan-2014  include in crispred as red_remap_image
;-

 ;;; raw centering of images
IF keyword_set(al) THEN BEGIN
    sh = shc(ref, pic)
    tmp = shift(pic, sh(0), sh(1))
ENDIF ELSE $
  tmp = pic

 ;;; find grid

s1 = size(ref)
s2 = size(pic)
 
sx = s1(1) < s2(1)
sy = s1(2) < s2(2)

 ;;; cell size must be integer
cs = fix(cs)
cs2 = cs/2

IF n_params() LT 4 THEN sw = cs2

 ;;; number of cells that fit in image
nx1 = (sx-cs)/sw+1
ny1 = (sy-cs)/sw+1

 ;;; startpoint of grid
x0 = (sx-(nx1-1)*sw-cs)/2+cs2
y0 = (sy-(ny1-1)*sw-cs)/2+cs2
 ;;; generate index matrices

xref = intarr(nx1+2, ny1+2)
yref = intarr(nx1+2, ny1+2)

xref(1:nx1, 1:ny1) = (indgen(nx1) # replicate(sw, ny1)) + x0
yref(1:nx1, 1:ny1) = (replicate(sw, nx1) # indgen(ny1)) + y0
xref(*, ny1+1) = xref(*, 1)
xref(*, 0) = xref(*, 1)
xref(nx1+1, *) = sx-1
yref(0, *) = yref(1, *)
yref(nx1+1, *) = yref(1, *)
yref(*, ny1+1) = sy-1

nn = (nx1+2)*(ny1+2)
xout = float(xref)
yout = float(yref)

IF keyword_set(sm) THEN BEGIN
    sh_ref = smooth(ref, sm, /edge)
    sh_tmp = smooth(tmp, sm, /edge)
ENDIF ELSE BEGIN
    sh_ref = ref
    sh_tmp = tmp
ENDELSE

FOR i=1, nx1 DO BEGIN
    FOR j=1, ny1 DO BEGIN
        sh = red_shc(sh_ref(xref(i, j)-cs2:xref(i, j)+cs2-1, $
                            yref(i, j)-cs2:yref(i, j)+cs2-1), $
                     sh_tmp(xref(i, j)-cs2:xref(i, j)+cs2-1, $
                            yref(i, j)-cs2:yref(i, j)+cs2-1), $
                     /interpolate, /filt, range = range)
        xout(i, j) = xref(i, j)+sh(0)
        yout(i, j) = yref(i, j)+sh(1)
    ENDFOR
ENDFOR
IF keyword_set(rs) THEN BEGIN
      ;;; recompute shifts if they are larger than rlim pixel
    IF NOT keyword_set(rlim) THEN rlim = 1
    rcs2 = (keyword_set(rcs) ? rcs/2 : cs2)
    shx = round(xout-xref)
    shy = round(yout-yref)
    xr1 = xref-shx
    yr1 = yref-shy
    nosh = (abs(shx) LT rlim) AND (abs(shy) LT rlim)
    FOR i=1, nx1 DO BEGIN
        FOR j=1, ny1 DO BEGIN
            IF nosh(i, j) THEN CONTINUE
            sh = red_shc(sh_ref(xref(i, j)-rcs2:xref(i, j)+rcs2-1, $
                                yref(i, j)-rcs2:yref(i, j)+rcs2-1), $
                         sh_tmp(xr1(i, j)-rcs2:xr1(i, j)+rcs2-1, $
                                yr1(i, j)-rcs2:yr1(i, j)+rcs2-1), $
                         /interpolate, /filt, range = range)
            xout(i, j) = xref(i, j)+shx(i, j)+sh(0)
            yout(i, j) = yref(i, j)+shy(i, j)+sh(1)
        ENDFOR
    ENDFOR
ENDIF

IF keyword_set(med) THEN BEGIN
    xout = xref+median(xout-xref, med)
    yout = yref+median(yout-yref, med)
ENDIF

IF keyword_set(dfb) THEN BEGIN
    xref = xref(1:nx1, 1:ny1)
    yref = yref(1:nx1, 1:ny1)
    xout = xout(1:nx1, 1:ny1)
    yout = yout(1:nx1, 1:ny1)
    nn = nx1*ny1
ENDIF

xref = xref(*)
yref = yref(*)
xout = xout(*)
yout = yout(*)

IF keyword_set(plot) THEN BEGIN
    IF keyword_set(same) THEN $
      tvscl, tmp $
    ELSE $
      show, tmp, tit='Distortion Map'
    
    FOR i=0l, nn-1 DO BEGIN
        plots, [xref(i), xref(i)+plot*(xout(i)-xref(i))], $ 
          [yref(i), yref(i)+plot*(yout(i)-yref(i))], /dev 
    ENDFOR
    empty
ENDIF

d_map = [[xref], [yref], [xout], [yout]]

return, warp_tri(xout, yout, xref, yref, tmp, QUINT=quint)

END
