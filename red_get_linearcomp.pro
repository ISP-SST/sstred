function red_get_linearcomp, iwav, pp, npar
                                ;
  dim = size(pp, /dimension)
  res = fltarr(dim[1], dim[2]) + 1.0
  if(npar le 2) then return, res
                                ;
  wavf = 1.0
                                ;
  for ii = 2L, npar - 1 do begin
     wavf *= iwav
     res += pp[ii,*,*] * wavf
  endfor
                                ;
  return, res
end
