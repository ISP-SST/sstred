function red_getstates_pinh, files, lam = lam
  nt = n_elements(files)
  state = strarr(nt)
  lam = strarr(nt)
                                ;
  for ii = 0L, nt - 1 do begin
     tmp = strsplit(file_basename(files[ii]),'.',/extract)
     state[ii] = tmp[1] + '.' + tmp[2] + '.' +tmp[3]
     lam[ii] = strmid(string(float(tmp[1]) * 1.e-10), 2)
  endfor
                                ;
  return, state
end
