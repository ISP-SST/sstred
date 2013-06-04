function red_initcmap, wav, dat, x0 = x0, x1 = x1
  dim = size(dat, /dim)
  nt3 = dim[0]-2
  cmap = fltarr(dim[1], dim[2])
                                ;
  if(n_elements(x0) eq 0) then x0 = 1
  if(n_elements(x1) eq 0) then x1 = nt3
                                ;
  for jj = 0L, dim[2] -1 do for ii = 0L, dim[1] - 1 do begin
     if(total(dat[*,ii,jj])  le 1.e-3) then continue
     dum = min(dat[x0:x1,ii,jj], p)
     
                                ;p = (p>x0)<x1
     idx = [p-1, p, p+1] + x0
     cc = parab_fit(wav[idx], dat[idx,ii,jj])
                                ;cc = poly_fit(wav[p-2:p+2], dat[p-2:p+2,ii,jj], 2)
     cmap[ii,jj] = -0.5 * cc[1] / cc[2]
  endfor
                                ;
  return, cmap
end
