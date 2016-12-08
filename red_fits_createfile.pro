; docformat = 'rst'

;+
; Create a fits file to write data cube slices into when you don't
; want to have the entire cube in memory at the same time.
; 
; Usage:
;
; mkhdr,head,type,naxisx   ; Make a FITS header appropriate for your data cube.
; sxaddpar,hdr,...         ; Add whatever information you want
; 
; red_fits_createfile, 'file.fits', head, lun, fileassoc [, tabhdu = tabhdu ]
; for i=0,N-1 do begin
;    im = ... ; Calculate image number i.
;    fileassoc[i]=im
; endfor
; free_lun,lun
;
; Binary extensions can be included by use of the tabhdu keyword.
;
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; 
; :Params:
; 
;    filename : in, type=string
;
;       The name of the file to create.
;
;
;    head : in, type=strarr
;
;       The complete fits header of the file. The NAXIS* and BITPIX
;       keywords are used to set the file up.
;
;
;    lun : out, type=integer
;
;       The logical unit of the opened file.
;
;
;    fileassoc : out, type="file assoc variable"
;
;       An assoc variable that can be used to write (and read!) data
;       slices. 
;
;    tabhdu : in, type=struct
; 
;       Binary extensions to be added to the fits file. See
;       documentation in red_writedata.pro.
; 
; 
; :History:
; 
;    2016-11-07 : MGL. First version.
; 
;    2016-12-07 : MGL. New keyword tabhdu.
; 
; 
; 
;-
pro red_fits_createfile, filename, head, lun, fileassoc, tabhdu = tabhdu

  ;; Open the new fits file
  openw, lun, filename, /get_lun, /swap_if_little_endian

  hdr = head                    ; Protect input
  if n_elements(tabhdu) gt 0 then begin
    ;; Binary extension header 
    extname = tag_names(tabhdu)
    for i = 0, n_elements(extname)-1 do red_create_extheaders, tabhdu, extname[i], hdr
  endif

  ;; Make byte-version of header and write it to the file.
  Nlines = n_elements(hdr)     ; Number of header lines
  bhdr = replicate(32B, 80L*Nlines)
  for n = 0L, Nlines-1 do bhdr[80*n] = byte(hdr[n])
  writeu, lun, bhdr

  ;; Pad to mod 2880 bytes
  Npad = 2880 - (80L*Nlines mod 2880) ; Number of bytes of padding
  if npad GT 0 then writeu, lun, replicate(32B,Npad)
  Nblock = (Nlines-1)*80/2880+1 ; Number of 2880-byte blocks
  offset = Nblock*2880          ; Offset to start of data

  ;; Now prepare for the data part

  naxis = sxpar(head,'NAXIS')
  dims = lonarr(2)
  Nelements = 1L
  Nslices = 1L
  for i = 1, 2 do begin
     dims[i-1] = sxpar(head,'NAXIS'+strtrim(i, 2))
     Nelements *= dims[i-1]
     if i gt 1 then Nslices *= dims[i-1]
  endfor
  
  bitpix = sxpar(head,'BITPIX' )
  case bitpix of
      16 : array_structure = intarr(dims)
     -32 : array_structure = fltarr(dims)
     else : stop
  endcase

  ;; Pad the file, including the data part, to mod 2880 bytes
  endpos = offset + Nelements*abs(bitpix)/8 ; End position of data
  Npad = 2880 - (endpos mod 2880)
  if npad GT 0 and npad LT 2880 then begin
     bpad = bytarr(Npad)
     rec = assoc(lun,bpad,endpos)
     rec[0] = bpad
  endif

  ;; Write binary extension if any.
  if n_elements(tabhdu) gt 0 then begin
    free_lun, lun               ; The file is expected to be closed for this.
    red_write_tabhdu, tabhdu, filename
    ;; Open the fits file again, so the image data can be written
    ;; outside of this subroutine.
    openu, lun, filename, /get_lun, /swap_if_little_endian
  endif

  ;; Set up an assoc variable for writing the slices
  fileassoc = assoc(lun, array_structure, offset)

end
