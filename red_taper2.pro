function red_taper2, nx, ny, sm
                                ;
                                ; Coded by Jaime, inspired on Mats Loefdahl's makewindow.pro
                                ;
  smx = long(sm * nx)
  smy = long(sm * ny)
  x = fltarr(nx) + 1.0
  y = fltarr(ny) + 1.0
                                ;
  xa = findgen(smx) 
  ya = findgen(smy)
                                ;
  spx = 0.5*(1-cos(2.0*!pi*xa/(2.0*smx-1.0)))
  spy = 0.5*(1-cos(2.0*!pi*ya/(2.0*smy-1.0)))
  x[0:smx-1] = spx
  x[nx-smx:*] = reverse(spx)
  y[0:smy-1] = spy
  y[ny-smy:*] = reverse(spy)
                                ;
  return, x#y
end
