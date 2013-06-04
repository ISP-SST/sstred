function pol_getborder, var, x0, x1, y0, y1, square = square
  inam = 'pol_getborder : '
  dim = size(var, /dimension)
  ver = total(var, 1) / dim[0]
  hor = total(var, 2) / dim[1]
                                ; Vertical axis
  pos = where(ver ne 0.)
  y0 = min(pos)
  y1 = max(pos)
                                ; Horzontal axis
  pos = where(hor ne 0.)
  x0 = min(pos)
  x1 = max(pos)
                                ;
  nx = x1 - x0 + 1
  ny = y1 - y0 + 1
                                ;
                                ; Preferred square image
  if(keyword_set(square)) then begin
     if((nx gt ny)) then begin
        dif = nx - ny
        x0+= dif / 2L
        x1-= dif / 2L
        if((dif/2)*2 ne dif) then x1-=1L
        nx = x1 - x0 + 1L
     endif else begin
        if(ny gt nx) then begin
           dif = ny - nx
           y0+= dif / 2L
           y1-= dif / 2L
           if((dif/2)*2 ne dif) then y1-=1L
           ny = y1 - y0 + 1L
        endif  
     endelse
  endif
                                ;
  print, inam + 'Clipping info:'
  print, '   x0 = '+red_stri(x0)+',  x1 = ' + red_stri(x1)
  print, '   y0 = '+red_stri(y0)+',  y1 = ' + red_stri(y1)
                                ;
  return, [nx, ny]
end
