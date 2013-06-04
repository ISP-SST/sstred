function red_flagchange, st
  nt = n_elements(st)
  star = bytarr(nt)
                                ;
  os = st[0]
  for ii = 0L, nt - 1 do begin
     if(st[ii] ne os) then begin
        star[ii] = 1B
        os = st[ii]
     endif
  endfor
                                ;
  return, star
end
