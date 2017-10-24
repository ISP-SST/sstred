; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
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
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2016-10-31 : MGL. Use red_findthedim and red_lp_header.
; 
;    2016-11-01 : MGL. Transpose entire cube without asking if enough
;                 memory. 
; 
; 
;-
pro red_flipthecube_unpol, file, Nw = Nw, Nt = Nt, maxsize=maxsize, icube = icube
  
  inam = 'flipthecube_unpol: '

  Nt = long64(Nt)
  Nw = long64(Nw)

  if(~file_test(file)) then begin
    print, inam + 'file not found -> ', file
    return
  endif
  red_lp_header, file, Nx = Nx, Ny = Ny, Nt = tmp

  Ntot = Nt * Nw 
  if(ntot ne tmp) then begin
    print, inam + 'wrong number of elements: Nw * Nt  != tmp'
    print, '  ', Ntot, tmp
    stop
  endif

;  Nbytes = Ntot * long(Nx) * long(Ny) * 4.
;  print, 'Total number of bytes in cube:', Nbytes

  ;; Compare cube size with available memory (free + cashed) in GB
  spawn,'free -g ',free
  free = total(long((strsplit(free[1],/extract))[[3,6]]))
  cubesize = (float(Ntot) * float(Nx) * float(Ny)) * 1e-9
  if keyword_set(icube) then cubesize *= 2. else cubesize *= 4.
  print, 'Free memory and cube size in GB:', free, cubesize

  Nx1 = Nx
  if(~keyword_set(maxsize)) then maxsize = 256L
  Ny1 = red_findthedim(Ny, maxsize)
  if Ny1 eq 1 then begin

    if free*0.75 gt cubesize then begin

      print, inam + ' : Enough memory to transpose entire cube.'
      Ny1=Ny

    endif else begin
      ;; If Ny1 is 1, making the sp cube takes a lot of time. If you
      ;; have sufficient memory, you might want to transpose the whole
      ;; cube in one go provide that option here
      print, inam + 'No optimal factor found to transpose smaller cubes (ny1=1)'
      print, inam + 'Try to transpose the whole cube in one go? [y/n] '
      dum=''
      read,dum
      if dum eq 'y' then ny1=ny
    endelse
  endif 
  
  
  openr,lun,file, /get_lun
  if(~keyword_set(icube)) then  $
     a = assoc(lun, fltarr(Nx1, Ny1, /noze), 512) $
  else $
     a = assoc(lun, intarr(Nx1, Ny1, /noze), 512)

  ext = strsplit(file, '.', /extract)
  ext = ext[n_elements(ext)-1]

  odir = file_dirname(file)+'/'
  ofile = odir+file_basename(file,'.'+ext)+'_sp.' + ext
  print, inam + 'saving result -> '+ofile
  openw, lun1, ofile, /get_lun
  if(~keyword_set(icube)) then $
     b = assoc(lun1, fltarr(Nw,Nt,/noze), 512) $
  else $
     b = assoc(lun1, intarr(Nw,Nt,/noze), 512)

  Ntimes = Ny / Ny1
  k = 0LL
  for ii = 0L, Ntimes - 1 do begin
    ;;print, 'processing part '+string(ii,format='(A,I3)')+' of '+string(ntimes-1,format='(A,I3)')
    if Ntimes ne 1 then $
       red_progressbar, ii, Ntimes $
                        , 'Flipping the cube '+strtrim(ny1, 2)+' rows at a time'

    if(~keyword_set(icube)) then $
       cub = fltarr(Nx1,Ny1,Nw,Nt) $
    else $
       cub = intarr(Nx1,Ny1,Nw,Nt)
    
    ;; Read
    for it = 0LL, Nt - 1 do begin
      if Ntimes eq 1 then $
         red_progressbar, it, Nt, 'Reading the cube'
      for iw = 0LL, Nw - 1 do begin
        ele = it * Nw * Ntimes + iw * Ntimes + ii
        cub[*,*,iw,it] = a[ele]             
      endfor                    ; iw
    endfor                      ; it
    
    if Ntimes eq 1 then print, inam + ' : Transposing the cube (this may take some time)'
    cub = transpose(temporary(cub), [2,3,0,1])

    ;; Write
    for iy = 0L, Ny1 - 1 do begin
      if Ntimes eq 1 then $
         red_progressbar, iy, Ny1, 'Writing the transposed cube'
      for ix = 0L, Nx1 - 1 do begin
        b[k] = cub[*,*,ix,iy]
        k += 1L
      endfor                    ; ix
    endfor                      ; iy
    cub = 0B
  endfor                        ; ii

  free_lun, lun
  free_lun, lun1

                                
  ;; Write header
  openu, lun1, ofile, /get_lun
  extra = ''                    ;'stokes=[I,Q,U,V], ns=4'
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
