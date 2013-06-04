function red_unpol_lpheader, nx, ny, nz, float = float


  ;;
  ;; File header
  ;;
  if(keyword_set(float)) then begin
     typestring = '(float)'
     datatype = 4
  endif else begin
     typestring = '(integer)'
     datatype = 2
  endelse
  header = 'datatype='+strtrim(datatype,2)+' '+typestring
  header = header + ', dims='+strtrim(3,2)
  header = header + ', nx='+strtrim(nx,2)
  header = header + ', ny='+strtrim(ny,2)
  header = header + ', nt='+strtrim(nz,2)
  if ((byte(1L, 0, 1))[0] eq 1) then endianstr = 'endian=l'  $ ; little endian
  else endianstr = 'endian=b'                                  ; big endian
  header = header + ', '+endianstr
  
  hh = bytarr(512)
  hh1 = byte(header)
  hh[0:n_elements(hh1)-1] = hh1


  return, hh
end
