; docformat = 'rst'

;+
; Define distinct colors for plotting.
; 
; Default, without keyword, is black (0), white (255) and 9 colours
; that are as distinct as possible in both normal and colour-blind
; vision, but also match well together. They stay distinct when
; printed on paper. Optimized preselections are available with colnum
; set to a value in the range 1 to 12.
; 
; :Categories:
;
;    SST pipeline
; 
; :Author:
; 
;     P.J.J. Tol, SRON
; 
; :Returns:
; 
;    Color codes in various formats, by default as an array of RGB triples.
; 
; :Params:
; 
;    colnum : in, type=integer
; 
;      The number of distinct colors wanted.
; 
; 
; :Keywords:
; 
;   gray_background : in, optional, type=boolean
; 
;     Sets colours 0, 255 and up to 4 main colours as used in SRON
;     presentations.
; 
;   hex : in, optional, type=boolean
; 
;     Return the RGB values as an array of hex strings.
; 
;   num : in, optional, type=boolean
; 
;     Return the RGB values as an array of long integers.
; 
; :History:
; 
;   2009-08 : P.J.J. Tol. Written.
; 
;   2009-11-23 : P.J.J. Tol. Added GRAY_BACKGROUND and SCREEN.
; 
;   2010-01-17 : P.J.J. Tol. Improved palette.
; 
;   2010-08-14 : P.J.J. Tol. Changed 8-colour set.
; 
;   2015-01-18 : MGL. Do not set colour table, return colors instead. 
; 
;   2024-09-24 : MGL. Make the number of colors a parameter instead of
;                a keyword.
; 
;   2024-10-07 : MGL. Make 12 the default number of colors if you
;                specify it out of range.
; 
;   2024-10-09 : MGL. Remove unused keyword screen.
; 
;-
function red_distinct_colors, colnum $
                              , gray_background = gray $
                              , hex = hex $
                              , num = num

  IF n_elements(colnum) gt 0 THEN BEGIN
    colnum = ROUND(colnum)
    IF KEYWORD_SET(gray) THEN BEGIN
      IF (colnum LT 1 OR colnum GT 4) THEN BEGIN
        MESSAGE, 'Color numbers for a gray background is 1..4, set to 4.', /INFORMATIONAL
        colnum = 4
      ENDIF
    ENDIF ELSE BEGIN
      print, colnum
      IF (colnum LT 1 OR colnum GT 12) THEN BEGIN
        MESSAGE, 'Color numbers are 1..12, colnum set to 12.', /INFORMATIONAL
        colnum = 12
      ENDIF
      print, colnum
    ENDELSE
  ENDIF ELSE BEGIN
    IF KEYWORD_SET(gray) THEN colnum = 4 ELSE colnum = 9
  ENDELSE
  
  ;; colour coordinates
  xarr = $
   [[12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], $
    [12,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], $
    [12,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0], $
    [12,  6,  5,  3,  0,  0,  0,  0,  0,  0,  0,  0], $
    [ 0,  1,  3,  5,  6,  0,  0,  0,  0,  0,  0,  0], $
    [ 0,  1,  3,  5,  6,  8,  0,  0,  0,  0,  0,  0], $
    [ 0,  1,  2,  3,  5,  6,  8,  0,  0,  0,  0,  0], $
    [ 0,  1,  2,  3,  4,  5,  6,  8,  0,  0,  0,  0], $
    [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  0,  0,  0], $
    [ 0,  1,  2,  3,  4,  5,  9,  6,  7,  8,  0,  0], $
    [ 0, 10,  1,  2,  3,  4,  5,  9,  6,  7,  8,  0], $
    [ 0, 10,  1,  2,  3,  4,  5,  9,  6, 11,  7,  8]]
  x = xarr[0:colnum-1,colnum-1]
  
  IF KEYWORD_SET(gray) THEN BEGIN
    red   = [128,255,255,100]
    green = [155,102,204,194]
    blue  = [200,102,102,4]
    red   = red[0:colnum-1]
    green = green[0:colnum-1]
    blue  = blue[0:colnum-1]
  ENDIF ELSE BEGIN
    red   = [ 51, 136,  68,  17, 153, 221, 204, 136, 170, 102, 102, 170,  68]
    green = [ 34, 204, 170, 119, 153, 204, 102,  34,  68,  17, 153,  68, 119]
    blue  = [136, 238, 153,  51,  51, 119, 119,  85, 153,   0, 204, 102, 170]
    red   = red[[x]]
    green = green[[x]]
    blue  = blue[[x]]
  ENDELSE

  N = n_elements(red)
  
  rgb = lonarr(3, N)
  rgb[0, *] = Red
  rgb[1, *] = Green
  rgb[2, *] = Blue
  
  
  ;; Make RGB hexadecimal strings.
  if keyword_set(hex) then begin  
    rgbhex = strarr(N)
    for i = 0, N-1 do rgbhex[i] = string(rgb[0, i]*65536L + rgb[1, i]*256L + rgb[2, i] $
                                         , format='(Z06)')
    return, rgbhex
  endif
  
  ;; Make long integer array, suitable for plotting.
  if keyword_set(num) then begin
    rgbnum = lonarr(N)
    for i = 0, N-1 do rgbnum[i] = rgb[0, i] + rgb[1, i]*256L + rgb[2, i]*65536L
    return, rgbnum
  endif
  
  ;; Return the triples
  return, rgb
  
END


n_colors = 12

colors = red_distinct_colors(n_colors, /num)


cgplot, [0], [0], /nodata, xrange = [0, 1], yrange = [0, n_elements(colors)+1]
for i = n_elements(colors)-1, 0, -1 do begin
  cgplot, !x.crange, [i, i]+1, color = colors[i], /over, thick = 15
endfor

end
