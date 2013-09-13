function findthedim, ny, maxs
  res = 0L
  for ii = 1L, maxs do begin
     if((ny / ii) * ii EQ ny) then res = ii
  endfor
  print, 'findthedim: using ny1 = '+string(res, format='(I5)')
  return, res
end
pro lp_header, filename, header=header, datatype=datatype, $
               dims=dims, nx=nx, ny=ny, nt=nt ,endian=endian
;
; extracts data parameters from header
; header should contain following entries :
;   datatype=4 (float), dims=2, nx=2027, ny=2042
; see lp_write.pro
;
 IF n_params() EQ 0 THEN BEGIN
    message, /info, 'lp_header, filename, header=header, datatype=datatype $'
    print, '                      , dims=dims, nx=nx, ny=ny, nt=nt'
    retall
 ENDIF

 if ARG_PRESENT(header) then printheader=0 else printheader=1

 openr, lur, filename, /get_lun
 rec = assoc(lur, bytarr(512))   ; first 512 bytes is header info
 header = string(rec[0])
 free_lun, lur
 if printheader eq 1 then print, header
 ; extract info from header
 len = strlen(header)
 ; datatype
 searchstring = 'datatype='
 pos = strpos(header, searchstring)
 if pos eq -1 then begin
     message, /info, 'unknown datatype'
     print, '  header: '+header
     retall
 endif
 datatype = fix(strmid(header, pos+strlen(searchstring), 1))
 ; number of dimensions :
 searchstring = 'dims='
 pos = strpos(header, searchstring)
 if pos eq -1 then begin
     message, /info, 'unknown number of dimensions'
     print, '  header: '+header
     retall
 endif
 dims = fix(strmid(header, pos+strlen(searchstring), 1))
 if (dims lt 2) or (dims gt 3) then begin
     message, /info, 'number of dimensions not supported'
     print, '  dimensions: ', dims
     retall
 endif
 ; number of pixels in x-direction
 searchstring = 'nx='
 pos = strpos(header, searchstring)
 if pos eq -1 then begin
     message, /info, 'unknown number of pixels in x-direction'
     print, '  header: '+header
     retall
 endif
 nx = fix(strmid(header, pos+strlen(searchstring), 4))
 ; number of pixels in y-direction
 searchstring = 'ny='
 pos = strpos(header, searchstring)
 if pos eq -1 then begin
     message, /info, 'unknown number of pixels in y-direction'
     print, '  header: '+header
     retall
 endif
 ny = fix(strmid(header, pos+strlen(searchstring), 4))
 ; number of pixels in t-dimension
; if dims eq 3 then begin
     searchstring = 'nt='
     pos = strpos(header, searchstring)
     if pos eq -1 then begin
         message, /info, 'unknown number of pixels in t-direction'
         print, '  header: '+header
         return
     endif else nt = long(strmid(header, pos+strlen(searchstring), 5))
; endif
     searchstring = 'endian='
     pos = strpos(header, searchstring)
     if pos eq -1 then begin
         message, /info, 'unknown endianness'
         print, '  header: '+header
         return
     endif else endian = strmid(header, pos+strlen(searchstring), 1)

end


pro red_flipthecube_unpol, file, nw = nw, nt = nt, maxsize=maxsize, icube = icube

  inam = 'flipthecube_unpol: '

  nt = long64(nt)
  nw = long64(nw)

  if(~file_test(file)) then begin
     print, inam + 'file not found -> ', file
     return
  endif
  lp_header, file, nx = nx, ny = ny, nt = tmp
  
  ntot = nt * nw 
  if(ntot ne tmp) then begin
     print, inam + 'wrong number of elements: nw * nt  != tmp'
     print, '  ', ntot, tmp
     stop
  endif
  
  nx1 = nx
  if(~keyword_set(maxsize)) then maxsize = 256L
  ny1 = findthedim(ny, maxsize)
  

  
  openr,lun,file, /get_lun
  if(~keyword_set(icube)) then a = assoc(lun, fltarr(nx1,ny1,/noze), 512) $
  else a = assoc(lun, intarr(nx1,ny1,/noze), 512)

  ext = strsplit(file,'.',/extract)
  ext = ext[n_elements(ext)-1]

  odir = file_dirname(file)+'/'
  ofile = odir+file_basename(file,'.'+ext)+'_sp.' + ext
  print, inam + 'saving result -> '+ofile
  openw, lun1, ofile, /get_lun
  if(~keyword_set(icube)) then b = assoc(lun1, fltarr(nw,nt,/noze), 512) $
  else b = assoc(lun1, intarr(nw,nt,/noze), 512)

  ntimes = ny / ny1
  k = 0LL
  for ii = 0L, ntimes - 1 do begin
     print, 'processing part '+string(ii,format='(A,I3)')+' of '+string(ntimes-1,format='(A,I3)')
     if(~keyword_set(icube)) then cub = fltarr(nx1,ny1,nw,nt) $
     else cub = intarr(nx1,ny1,nw,nt)
     
     ;
     ; read
     ;
     for tt = 0LL, nt-1 do begin
        for ww = 0LL, nw-1 do begin
           ele = tt * nw * ntimes + ww * ntimes + ii
           cub[*,*,ww,tt] = a[ele]             
        endfor
     endfor
     
     cub = transpose(temporary(cub), [2,3,0,1])
     for yy = 0L, ny1 - 1 do for xx=0L, nx1-1 do begin
        b[k] = cub[*,*,xx,yy]
        k += 1L
     endfor
     cub = 0B
  endfor
  free_lun, lun
  free_lun, lun1

  ;
  ; Write header
  ;
  openu, lun1, ofile, /get_lun
  extra = '';'stokes=[I,Q,U,V], ns=4'
  if(~keyword_set(icube)) then typestring = '(float)' else typestring = '(integer)'
  if(~keyword_set(icube)) then datatype = 4 else datatype = 2
  header = 'datatype='+strtrim(datatype,2)+' '+typestring
  header = header + ', dims='+strtrim(3,2)
  header = header + ', nx='+strtrim(nw,2)
  header = header + ', ny='+strtrim(nt,2)
  header = header + ', nt='+strtrim(long64(nx)*long64(ny),2)
  if ((byte(1L, 0, 1))[0] eq 1) then endianstr = 'endian=l'  $ ; little endian
  else endianstr = 'endian=b'                                  ; big endian
  header = header + ', '+endianstr
  hh = bytarr(512)
  hh1 = byte(header)
  hh[0:n_elements(hh1)-1] = hh1
  writeu, lun1, hh
  free_lun, lun1
  ;
end
