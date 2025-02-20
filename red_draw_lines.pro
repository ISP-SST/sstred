; docformat = 'rst'

;+
; Draw lines in a graphics window.
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
; 
; :Params:
; 
;    x : in, type=array
; 
;      The horizontal coordinates of points, between which the lines
;      are to be drawn. (But see keywords for how this can be
;      modified.)
; 
;    y : in, type=array
; 
;      The vertical coordinates of points, between which the lines
;      are to be drawn. (But see keywords for how this can be
;      modified.)
; 
; 
; :Keywords:
; 
;    axes : in, optional, type=boolean
;
;       Interpret x and y as the coordinates at which to draw axes.
;
;    close : in, optional, type=boolean
;
;       Add a line to close the curve.
;
;    hull : in, optional, type=boolean
;
;       Draw the convex hull of the nodes defined by x and y.
;
;    minmax : in, optional, type=boolean
;
;       Draw a rectangle with corners determined by the minimum and
;       maximum values of x and y.
;
;    _extra : in, optional
;   
;      Other keywords to be sent to cgPolyLines.
; 
; 
; :History:
; 
;   2020-04-24 : MGL. First version.
; 
;-
pro red_draw_lines, x, y $
                    , axes = axes $
                    , close = close $
                    , hull = hull $
                    , minmax = minmax $
                    , _extra = extra

  case 1 of

    keyword_set(axes) : begin
      if keyword_set(extra.device) then begin
        minx = 0
        maxx = !d.x_size
        miny = 0
        maxy = !d.y_size
      endif else begin
        ;; Not device coordinates, assume data coordinates.
        minx = !x.crange[0]
        maxx = !x.crange[1]
        miny = !y.crange[0]
        maxy = !y.crange[1]
      endelse
      x_draw = [minx, maxx, x, x, x]
      y_draw = [y, y, y, miny, maxy]
    end

    keyword_set(close) : begin
      x_draw = [x, x[0]]
      y_draw = [y, y[0]]
    end

    keyword_set(hull) : begin
      ;; The documentation of qhull suggests that one should use the
      ;; bounds keyword to get the indices of the convex hull nodes.
      ;; That seems to just get the value -1 but the tr parameter
      ;; seems to have the info needed.
      qhull, x, y, tr           ; Calculate the convex hull
      indx = tr[0, *]           ; Indices of the convex hull nodes
      ;; The convex hull nodes are not sorted, so we have to do that
      ;; here. Calculate a point indside the hull and then sort in
      ;; angular order around that point.
      xc = mean(x[indx])
      yc = mean(y[indx])
      theta = atan(x[indx]-xc, y[indx]-yc)
      indx = indx[sort(theta)]
      x_draw = [x[indx], x[indx[0]]]
      y_draw = [y[indx], y[indx[0]]]
    end

    keyword_set(minmax) : begin
      x_draw = [min(x), max(x), max(x), min(x), min(x)]
      y_draw = [min(y), min(y), max(y), max(y), min(y)]
    end

    else : begin
      x_draw = x
      y_draw = y
    end
    
  endcase

  if n_elements(extra) eq 0 || total(tag_names(extra) eq 'DEVICE') eq 0 then begin
    cgPolygon, x_draw, y_draw, _strict_extra = extra, color = 'opposite'
  endif else begin 
    cgPolygon, x_draw, y_draw, _strict_extra = extra
  endelse
  
end
