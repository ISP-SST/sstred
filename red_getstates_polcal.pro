function red_getstates_polcal, files
                                ;
  nt = n_elements(files)
                                ;
  pstat = {state:strarr(nt), lp:strarr(nt), qw:strarr(nt), pref:strarr(nt),$
           wav:strarr(nt), lc:strarr(nt), nums:strarr(nt), files:files, star:bytarr(nt)}
                                ;
  for ii = 0L, nt - 1 do begin
     tmp = strsplit(file_basename(files[ii]), '.', /extract)
     pstat.lp[ii] = tmp[4]
     pstat.qw[ii] = tmp[5]
     pstat.pref[ii] = tmp[6]
     pstat.wav[ii] = tmp[7]
     pstat.lc[ii] = tmp[8]
     pstat.nums[ii] = tmp[10]
     pstat.state[ii] = tmp[4]+'.'+tmp[5]+'.'+tmp[6]+'.'+tmp[7]+'.'+tmp[8]
  endfor
                                ;
  return, temporary(pstat)
end
