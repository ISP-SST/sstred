function red_getstates_polcal_out, files
                                ;
  nt = n_elements(files)
                                ;
  stat = {state:strarr(nt), lp:fltarr(nt), qw:fltarr(nt), pref:strarr(nt), cam:strarr(nt),$
          wav:strarr(nt), lc:bytarr(nt), lps:strarr(nt), qws:strarr(nt), lcs:strarr(nt)}
                                ;
  for ii = 0L, nt  -1 do begin
     tmp = strsplit(file_basename(files[ii]), '.', /extract)
     stat.state[ii] = tmp[1] + '.' + tmp[2] + '.' + tmp[3] + '.' + tmp[5]
     stat.lp[ii] = float(strmid(tmp[1], 2, 3))
     stat.qw[ii] = float(strmid(tmp[2], 2, 3))
     stat.pref[ii] = tmp[3]
     stat.wav[ii] = tmp[4]
     stat.lc[ii] = strmid(tmp[5],2,1)
     stat.lcs[ii] = tmp[5]
     stat.lps[ii] = tmp[1]
     stat.qws[ii] = tmp[2]
     stat.cam[ii] = tmp[0]
  endfor
                                ;
  return, stat
end
