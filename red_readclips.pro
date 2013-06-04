pro red_readclips, file, tc = tc, rc = rc, wc = wc
  inam = 'red_readclips : '
  if(~file_test(file)) then begin
     print, inam + 'ERROR -> file not found '+file
     STOP
  endif
  openr, lun, file, /get_lun
  tc = ' '
  rc = ' '
  wc = ' '
  readf, lun, wc
  readf, lun, tc
  readf, lun, rc
  free_lun, lun 
  return
end
