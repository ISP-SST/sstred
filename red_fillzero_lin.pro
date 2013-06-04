function red_fillzero_lin, var
  idx = where(var LE 1.d-5, count)
  if count eq 0 then begin
     print, 'not fixing data!'
     return, var
  endif
;
  for ii = 0L, count-1 do begin
     jj = idx[ii]
                                ;print,max(var[[jj-1,jj+1]])

     if(max(var[[jj-1,jj+1]]) GT 1.e-5) then var[jj] = 0.5 * (var[jj-1] + var[jj+1])
  endfor

  return, var
end
