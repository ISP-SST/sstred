pro red_flagtunning, stat, nremove
  if(n_elements(nremove) eq 0) then nremove=1
  if nremove eq 0 then return
                                ;
  os = stat.wav[0]
                                ;
  for ii = 1L, n_elements(stat.wav) - 1 do begin
     if(stat.wav[ii] ne os) then begin
                                ;print, os, stat.state[ii]
        os = stat.wav[ii]
        stat.star[ii:ii:ii+nremove-1] = 1B
     endif
  endfor
                                ;
  return
end
