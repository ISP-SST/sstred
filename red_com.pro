function red_com, data

  on_error, 2
  
  s = size(data)
  
  n = s[0]
  sx = s[1]
  sy = s[2]
  sz = s[3]
  rx = findgen(sx)
  
  td = total(data)
  
  IF n EQ 1 THEN return, total(data*rx)/td
  ry = transpose(findgen(sy)) 
  IF n EQ 2 THEN BEGIN
     x1 = total(data*rebin(rx, sx, sy))/td
     y1 = total(data*rebin(ry, sx, sy))/td
     cm= [x1, y1]
  ENDIF
  IF n EQ 3 THEN BEGIN
     rz = fltarr(1, 1, sz) & rz[*] = findgen(sz)
     x1 = total(data*rebin(rx, sx, sy, sz))/td
     y1 = total(data*rebin(ry, sx, sy, sz))/td
     z1 = total(data*rebin(rz, sx, sy, sz))/td
     cm =  [x1, y1, z1]
  ENDIF
  return, cm
END
;
pro red_readclips, file, tc = tc, rc = rc, wc = wc
  inam = 'red_readclips : '
  if(~file_test(file)) then begin
     print, inam + 'ERROR -> file not found '+file
     STOP
  endif
  openr, lun, file, /get_lun
  tc = ' '
  rc = ' '
  wc = ' '
  readf, lun, wc
  readf, lun, tc
  readf, lun, rc
  free_lun, lun 
  return
end
