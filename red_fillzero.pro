function red_fillzero, var
                                ;
  idx = where(var lt 1.e-5, count, complement = idx1)
  if(count gt 1) then var[idx] = median(var[idx1])
                                ;
  return, var
end
