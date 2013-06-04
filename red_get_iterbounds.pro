function red_get_iterbounds,tot,np
  maxi=double(tot) / double(np)
  maxit=fix(maxi)
  jot=fltarr(maxit+1)+np
  jot[maxit]=(maxi-maxit)*np
  if jot[maxit] eq 0 then jot=jot[0:maxit-1]
  return,jot
end
