function red_lp_read, filename, header=header
;
; read La Palma data stored as IDL assoc file
;
; First 512 bytes are assumed to contain header info containing:
;    datatype=2 (integer), dims=2, nx=2027, ny=2042
;
; see lp_write for writing right header style
;
  if n_params() eq 0 then begin
    message, /info, 'image = lp_read(filename [, header=header])'
    retall
  endif

  if arg_present(header) then printheader=0 else printheader=1
  ;; print header if header is not asked for as keyword
  
  red_lp_header, filename, header=header, datatype=datatype, $
                 dims=dims, nx=nx, ny=ny, nt=nt, endian=endian_file

  if ((byte(1L, 0, 1))[0] eq 1) then endian = 'l' else endian='b'
  if(datatype gt 1) and (endian ne endian_file) then swap_endian=1 else swap_endian=0
  openr, lur, filename, /get_lun, swap_endian=swap_endian
; rec = assoc(lur, bytarr(512))   ; first 512 bytes is header info
; header = string(rec[0])
; ; extract info from header
; len = strlen(header)
; ; datatype
; searchstring = 'datatype='
; pos = strpos(header, searchstring)
; if pos eq -1 then begin
;     message, /info, 'unknown datatype'
;     print, '  header: '+header
;     free_lun, lur
;     return, 0
; endif
; datatype = fix(strmid(header, pos+strlen(searchstring), 1))
; ; number of dimensions :
; searchstring = 'dims='
; pos = strpos(header, searchstring)
; if pos eq -1 then begin
;     message, /info, 'unknown number of dimensions'
;     print, '  header: '+header
;     free_lun, lur
;     return, 0
; endif
; dims = fix(strmid(header, pos+strlen(searchstring), 1))
; if (dims lt 2) or (dims gt 3) then begin
;     message, /info, 'number of dimensions not supported'
;     print, '  dimensions: ', dims
;     free_lun, lur
;     return, 0
; endif
; ; number of pixels in x-direction
; searchstring = 'nx='
; pos = strpos(header, searchstring)
; if pos eq -1 then begin
;     message, /info, 'unknown number of pixels in x-direction'
;     print, '  header: '+header
;     free_lun, lur
;     return, 0
; endif
; nx = fix(strmid(header, pos+strlen(searchstring), 4))
; ; number of pixels in y-direction
; searchstring = 'ny='
; pos = strpos(header, searchstring)
; if pos eq -1 then begin
;     message, /info, 'unknown number of pixels in y-direction'
;     print, '  header: '+header
;     free_lun, lur
;     return, 0
; endif
; ny = fix(strmid(header, pos+strlen(searchstring), 4))
; ; number of pixels in t-dimension
; if dims eq 3 then begin
;     searchstring = 'nt='
;     pos = strpos(header, searchstring)
;     if pos eq -1 then begin
;         message, /info, 'unknown number of pixels in t-direction'
;         print, '  header: '+header
;         free_lun, lur
;         return, 0
;     endif
;     nt = fix(strmid(header, pos+strlen(searchstring), 4))
; endif
  if printheader then message, /info, header

                                ; read actual data
  if dims eq 2 then begin       ; 2D case
    case datatype of
      1: begin 
        rec = assoc(lur, bytarr(nx,ny), 512)
        image = rec[0]
      end
      2: begin
        rec = assoc(lur, intarr(nx,ny), 512)
        image = rec[0]
      end
      3: begin
        rec = assoc(lur, lonarr(nx,ny), 512)
        image = rec[0]
      end
      4: begin
        rec = assoc(lur, fltarr(nx,ny), 512)
        image = rec[0]
      end
      else: begin
        message, /info, 'datatype not supported '
        print, ' datatype = ', datatype
        free_lun, lur
        return, 0
      end
    end
    free_lun, lur
    return, image
  endif else begin              ; 3D case
    case datatype of
      1: begin 
        rec = assoc(lur, bytarr(nx,ny,nt), 512)
        cube = rec[0]
      end
      2: begin
        rec = assoc(lur, intarr(nx,ny,nt), 512)
        cube = rec[0]
      end
      3: begin
        rec = assoc(lur, lonarr(nx,ny,nt), 512)
        cube = rec[0]
      end
      4: begin
        rec = assoc(lur, fltarr(nx,ny,nt), 512)
        cube = rec[0]
      end
      else: begin
        message, /info, 'datatype not supported '
        print, ' datatype = ', datatype
        free_lun, lur
        return, 0
      end
    end
    free_lun, lur
    return, cube
  endelse

end
