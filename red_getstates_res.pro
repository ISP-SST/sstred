function red_getstates_res, files
  nt = n_elements(files)
  stat = {pref:strarr(nt), wav:strarr(nt), scan:strarr(nt), lcs:strarr(nt), lc:bytarr(nt),$
          cam:strarr(nt), wb:bytarr(nt), state:strarr(nt), camwb:' '}
                                ;
  for ii = 0L, nt - 1 do begin
     tmp = strsplit(file_basename(files[ii]), '.',/extract)
     n = n_elements(tmp)
     stat.cam[ii] = tmp[0]
     stat.scan[ii] = tmp[1]
     stat.pref[ii] = tmp[2]
     if(n gt 4) then begin
        stat.wav[ii] = tmp[3]
        stat.lcs[ii] = tmp[4]
        stat.lc[ii] = strmid(tmp[4], 2, 1)
        stat.state[ii] = tmp[1]+'.'+tmp[2]+'.'+tmp[3]
     endif else begin
        stat.wb[ii] = 1B
        stat.camwb = tmp[0]
     endelse
  endfor
  return, stat
end
