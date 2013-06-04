function red_inum, file_list
  res = strarr(n_elements(file_list))
                                ;
  for i = 0L, n_elements(res) - 1 do begin
     tmp = strsplit(file_list[i],'.',/extract)
     res[i] = tmp[n_elements(tmp)-1]
  endfor
                                ;
  return, res
end
