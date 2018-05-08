function red_get_ctalk, im, idx=idx

  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  res = fltarr(4)
  dim = size(im,/dim)

  if(n_elements(idx) eq 0) then idx = indgen(dim[3])
  for ss = 1, 3 do begin
    res[ss] = total(im[*,*,ss,idx] * im[*,*,0, idx], /double) / total( im[*,*,0, idx]^2, /double)
  endfor
  
  print, inam+': crosstalk from I -> Q,U,V=',res[1],', ',res[2],', ',res[3], format='(A,F8.5,A,F8.5,A,F8.5)'

  return, res

end
