function red_get_ctalk, im, idx=idx, mask=mask
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if(n_elements(mask) eq 0) then mask=1.0
  
  dim = size(im[*,*,0], /dim)
  tmask = fltarr([dim, n_elements(idx)])
  for ii=0, n_elements(idx)-1 do tmask[*,*,ii] = mask
  
  
  res = fltarr(4)
  dim = size(im,/dim)

  if(n_elements(idx) eq 0) then idx = indgen(dim[3])
  for ss = 1, 3 do begin
     res[ss] = total(reform(im[*,*,ss,idx] * im[*,*,0, idx]) * tmask, /double) / total(tmask*reform(im[*,*,0, idx])^2, /double)
  endfor
  
  print, inam+': crosstalk from I -> Q,U,V=',res[1],', ',res[2],', ',res[3], format='(A,F8.5,A,F8.5,A,F8.5)'

  return, res
end
