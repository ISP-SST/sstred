function red_sortfiles, files
  nt = n_elements(files)
  num = lonarr(nt)
  for ii = 0L, nt -1 do begin
     tmp = strsplit(files[ii],'.',/extract)
     num[ii] = long(tmp[n_elements(tmp) - 1])
  endfor
  pos = sort(num)
  res = files[pos]
  return, res
end
