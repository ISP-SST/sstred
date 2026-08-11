PRO red_tvimg, image, xax, yax, position=pos, box=box, noerase=noerase, $
           nolabels=nolabels, noscale=noscale, ASPECT=aspect, $
           SAMPLE=sample, TRUE=true, USE_COLORTABLE=usect, mn=mn, mx=mx, _EXTRA=extra

;+
; NAME:
;       TVIMG
; PURPOSE:
;       Display an image in an window with border and axis
; CALLING SEQUENCE:
;       TVIMG, IMAGE [, XAX [, YAX] [, KEYWORDS]]
; INPUTS:
;       IMAGE : Image to display
; OPTIONAL INPUT PARAMETER:
;       XAX   : (input) Xaxis values
;       YAX   : (input) Yaxis values
; KEYWORDS:
;       POSITION: Position of Plot area corners in norm coordinates. 
;                 Standard graphics keywords
;       BOX     : (input) Draw the koordinate axes with linethick
;                 BOX. Makes the plot look a bit Displa-Style.
;       ASPECT  : (Flag) Keep aspect ratio of the image
;       NOERASE : (Flag) Don't clear plot window before drawing
;       NOLABELS: (Flag) Don't draw axis annotations
;       NOSCALE : (Flag) By default, TVSCL is used to display the
;                 image. Set this keyword to use TV instead. Remember
;                 to do bytescaling yourself!
;       SAMPLE  : (Flag) use nearest neighbour sampling for image
;                 scaling in X-mode.
;       MIN     : Min value for bytescaling
;       MAX     : Max value for bytescaling
;       USE_COLORTABLE: (Flag) use previous loaded colortable to display image
;       Additionally, all valid keywords for the plot routine are
;       allowed.
; RESTRICTIONS:
;       
; PROCEDURE:
;       Eventualy rescale the image to fit the display, then oplot an
;       empty co-ordinate system. Draw thick bounding box if required.
; MODIFICATION HISTORY:
;       24-Sep-1994  P.Suetterlin, KIS
;       23-Feb-1995  Correct handling of multiple-plot styles (!P.style<>0)
;       25-Jul-1996  Use Keyword _EXTRA to pass keywords to plot
;                    routine. 
;       29-Oct-1997  Add keyword ASPECT to preserve the aspect ratio of
;                    the image
;       27-Nov-2001  Remove separate call to axis, so that keywords
;                    affecting the axis layout really work.
;       23-Sep-2002  Keyword SAMPLE
;       13-May-2003  Implement true color images; replace rescale with congrid
;-

on_error, 2

IF n_params() LT 1 THEN $
  message, 'Syntax: TVIMG, Image [, xax, yax]'

s = size(image)
IF keyword_set(true) THEN BEGIN
    CASE true OF
        1: BEGIN & sx = s(2) & sy=s(3) & END
        2: BEGIN & sx = s(1) & sy=s(3) & END
        3: BEGIN & sx = s(1) & sy=s(2) & END
    ENDCASE
ENDIF ELSE BEGIN
    sx = s(1) & sy=s(2)
ENDELSE

IF n_params() LT 3 THEN $
  yax = indgen(sy)

IF n_params() LT 2 THEN $
  xax = indgen(sx)

;;;
;;; We do an empty plot to
;;;  1) clear the screen (except NOERASE is set)
;;;  2) Set !x.window and !y.window if it is a multiple plot
;;;

plot, [0, 1], /nodata, xsty = 4, ysty = 4, xtit = '', ytit = '', $
 subtit = '', tit = '', noerase=noerase

IF NOT keyword_set(pos) THEN BEGIN 
  IF (!X.window(1)-!X.window(0)) EQ 0 THEN $
    pos = [0.1, 0.1, 0.95, 0.95] $
  ELSE $
    pos = [!X.window(0), !Y.window(0), !X.window(1), !Y.window(1)]
ENDIF

IF keyword_set(aspect) THEN BEGIN
     ;;; current aspect ratio
    asp = float(sx)/sy
     ;;; aspect ratio of the plot area
    asp1 = (pos(2)-pos(0))*!d.x_size/((pos(3)-pos(1))*!d.y_size)
    IF asp LT asp1 THEN BEGIN
         ;;; area is broader than pic -> shrink area horizontaly
        nw = asp/asp1*(pos(2)-pos(0))
        pos(0) = pos(0)+((pos(2)-pos(0))-nw)/2
        pos(2) = pos(0)+nw
    ENDIF ELSE BEGIN
         ;;; area is higher than pic -> shrink area vertically
        nw = asp1/asp*(pos(3)-pos(1))
        pos(1) = pos(1)+((pos(3)-pos(1))-nw)/2
        pos(3) = pos(1)+nw
    ENDELSE
ENDIF
    
IF !D.name NE 'PS' THEN GOTO, x

IF keyword_set(mn) and keyword_set(mx) THEN BEGIN
  noscale = 0
  image = BYTSCL(image, min = mn, max = mx)
ENDIF

IF keyword_set(noscale) THEN BEGIN
    tv, image, pos(0), pos(1), xsize=pos(2)-pos(0), $
      ysize=pos(3)-pos(1), /norm, TRUE=true
ENDIF ELSE BEGIN
    tvscl, image, pos(0), pos(1), xsize=pos(2)-pos(0), $
      ysize=pos(3)-pos(1), /norm, TRUE=true
ENDELSE

IF keyword_set(nolabels) THEN BEGIN
    plot, xax([0, sx-1]), yax([0, sy-1]), /nodata, /noerase, pos=pos, $
          /xsty, /ysty, xtickname=replicate(' ', 15), $
          ytickname=replicate(' ', 15), _EXTRA=extra
ENDIF ELSE BEGIN
    plot, xax([0, sx-1]), yax([0, sy-1]), /nodata, /noerase, pos=pos, $
          /xsty, /ysty, _EXTRA=extra
ENDELSE

IF keyword_set(box) THEN BEGIN
    lu = [!X.crange(0), !Y.crange(0)]
    ro = [!X.crange(1), !Y.crange(1)]
    plots, [lu(0), ro(0), ro(0), lu(0), lu(0)], $
      [lu(1), lu(1), ro(1), ro(1), lu(1)], thick = box
ENDIF 


return

X:
nx = fix((pos(2)-pos(0))*!D.x_size+0.5)+1
ny = fix((pos(3)-pos(1))*!D.y_size+0.5)+1

IF keyword_set(true) THEN BEGIN
    CASE true OF
        1: im = congrid(image, 3, nx, ny, INTERP=keyword_set(sample) EQ 0)
        2: im = congrid(image, nx, 3, ny, INTERP=keyword_set(sample) EQ 0)
        3: im = congrid(image, nx, ny, 3, INTERP=keyword_set(sample) EQ 0)
    ENDCASE
ENDIF ELSE $
  im = congrid(image, nx, ny, INTERP=keyword_set(sample) EQ 0)

IF NOT keyword_set(noscale) THEN im = bytscl(temporary(im))

IF keyword_set(usect) THEN BEGIN
    tvlct, r, g, b, /get
    tv, r(im), pos(0), pos(1), 1, /norm
    tv, g(im), pos(0), pos(1), 2, /norm
    tv, b(im), pos(0), pos(1), 3, /norm
ENDIF ELSE $
  tv, im, pos(0), pos(1), /norm, TRUE=true

IF keyword_set(nolabels) THEN BEGIN
    plot, xax([0, sx-1]), yax([0, sy-1]), /nodata, /noerase, pos=pos, $
          /xsty, /ysty, xtickname=replicate(' ', 15), $
          ytickname=replicate(' ', 15), _EXTRA=extra
ENDIF ELSE BEGIN
    plot, xax([0, sx-1]), yax([0, sy-1]), /nodata, /noerase, pos=pos, $
          /xsty, /ysty, _EXTRA=extra
ENDELSE

;IF NOT keyword_set(nolabels) THEN BEGIN
;    axis, xax=0, /xsty
;    axis, yax=0, /ysty
;ENDIF

IF keyword_set(box) THEN BEGIN
    lu = [!X.crange(0), !Y.crange(0)]
    ro = [!X.crange(1), !Y.crange(1)]
    plots, [lu(0), ro(0), ro(0), lu(0), lu(0)], $
      [lu(1), lu(1), ro(1), ro(1), lu(1)], thick = box
ENDIF 

END

