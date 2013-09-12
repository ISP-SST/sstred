; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    file : 
;
;
; :Keywords:
; 
;   nw : 
;   
;   
;   nt : 
;   
;   
;   maxsize : 
;   
;   
;   icube : 
; 
; 
; 
; 
; :History:
; 
;    2013-09-12 : MGL. Brought into red_ namespace.
; 
;-
pro red_flipthecube, file, nw = nw, nt = nt, maxsize=maxsize, icube = icube

  inam = 'flipthecube: '

  nt = long64(nt)
  nw = long64(nw)

  if(~file_test(file)) then begin
     print, inam + 'file not found -> ', file
     return
  endif
  lp_header, file, nx = nx, ny = ny, nt = tmp
  
  ntot = nt * nw * 4LL
  if(ntot ne tmp) then begin
     print, inam + 'wrong number of elements: nw * nt * 4L != tmp'
     return
  endif
  
  nx1 = nx
  if(~keyword_set(maxsize)) then maxsize = 256L
  ny1 = red_findthedim(ny, maxsize)
  
  openr,lun,file, /get_lun
  if(~keyword_set(icube)) then a = assoc(lun, fltarr(nx1,ny1,/noze), 512) $
  else a = assoc(lun, intarr(nx1,ny1,/noze), 512)

  ext = strsplit(file,'.',/extract)
  ext = ext[n_elements(ext)-1]

  odir = file_dirname(file)+'/'
  ofile = odir+file_basename(file,'.'+ext)+'_sp.' + ext
  print, inam + 'saving result -> '+ofile
  openw, lun1, ofile, /get_lun
  if(~keyword_set(icube)) then b = assoc(lun1, fltarr(nw,nt,4,/noze), 512) $
  else b = assoc(lun1, intarr(nw,nt,4,/noze), 512)

  ntimes = ny / ny1
  k = 0LL
  for ii = 0L, ntimes - 1 do begin
     print, 'processing part '+string(ii,format='(A,I3)')+' of '+string(ntimes-1,format='(A,I3)')
     if(~keyword_set(icube)) then cub = fltarr(nx1,ny1,nw,4,nt) $
     else cub = intarr(nx1,ny1,nw,4,nt)
     
     ;; Read
     for tt = 0LL, nt-1 do begin
        for ss = 0LL, 3 do begin
           for ww = 0LL, nw-1 do begin
              ele = tt * 4L * nw * ntimes + ss * nw * ntimes + ww * ntimes + ii
              cub[*,*,ww,ss,tt] = a[ele]             
           endfor
        endfor
     endfor
     
     cub = transpose(temporary(cub), [2,4,3,0,1])
     for yy = 0L, ny1 - 1 do for xx=0L, nx1-1 do begin
        b[k] = cub[*,*,*,xx,yy]
        k += 1L
     endfor
     cub = 0B
  endfor
  free_lun, lun
  free_lun, lun1

  ;; Write header
  openu, lun1, ofile, /get_lun
  extra = 'stokes=[I,Q,U,V], ns=4'
  if(~keyword_set(icube)) then typestring = '(float)' else typestring = '(integer)'
  if(~keyword_set(icube)) then datatype = 4 else datatype = 2
  header = 'stokes=[I,Q,U,V], ns=4 : ' +' datatype='+strtrim(datatype,2)+' '+typestring
  header = header + ', dims='+strtrim(3,2)
  header = header + ', nx='+strtrim(nw,2)
  header = header + ', ny='+strtrim(nt,2)
  header = header + ', nt='+strtrim(long64(nx)*long64(ny)*4LL,2)
  if ((byte(1L, 0, 1))[0] eq 1) then endianstr = 'endian=l'  $ ; little endian
  else endianstr = 'endian=b'                                  ; big endian
  header = header + ', '+endianstr
  hh = bytarr(512)
  hh1 = byte(header)
  hh[0:n_elements(hh1)-1] = hh1
  writeu, lun1, hh
  free_lun, lun1

end
