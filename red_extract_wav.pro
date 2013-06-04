function red_extract_wav, state, lc = lc
  nt = n_elements(state)
  res = dblarr(nt)
  lc = strarr(nt)
  for ii = 0L, nt - 1 do begin
     tmp = strsplit(state[ii],'._',/extract)
     res[ii] = double(tmp[1]) + double(tmp[2]) * 0.001d0
     lc[ii] = tmp[3]
  endfor
  return, res
end
