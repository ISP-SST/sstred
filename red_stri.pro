function red_stri, var, ni = ni
  if(n_elements(ni) eq 0) then begin
     res = strcompress(string(var), /remove_all)
  endif else begin
     res = string(var, FORMAT = ni)
  endelse
  return, res
end
