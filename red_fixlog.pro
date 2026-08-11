function readlog,logn,nrow
 ;OPENS LOG FILE AND CHECKS THE NUMBER OF ROWS 
  openr, lu, logn, /get_lun
  i=0L
  a=''
  while not EOF(lu) do begin 
     i=i+1L
     readf,lu,a
  endwhile
  nrow=i
  point_lun,lu,0
  print,logn+': '+strcompress(string(i),/remove_all)+' LINES'

  log=strarr(i)
  readf,lu,log
  free_lun,lu
  return,log
end

pro red_fixlog, inlog, outlog
  inam = 'red_fixlog: '

  ;; Read logfile
  if(~file_test(inlog)) then begin
     print, inam+'ERROR, file not found -> '+inlog
     return
  endif
  log1 = readlog(inlog, ct)
  log = log1
  
  ;; split into columns
  for ii = 0L, ct-1 do begin
     tmp = strsplit(log[ii], ' ', /extract)
     tmp[9] = string(float(tmp[8]) - float(tmp[9]), format='(F8.3)')
     tmp[10] = string((float(tmp[10]) + 180 ) mod 360. , format='(F10.5)')
     log[ii] = strjoin(tmp, ' ')
  endfor

  print, inam+'Saving '+outlog
  openw, lun, outlog, /get_lun
  printf, lun, log
  free_lun, lun
  
end
