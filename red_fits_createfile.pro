; docformat = 'rst'

;+
; Create a fits file to write data cube slices into.
; 
; Once the file is created, write slices to the file with
; fileassoc[islice]=data, where data has to be of the correct type and
; size.
;
; Close the file with free_lun,lun when done.
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
;
; 
; 
; :History:
; 
;    2016-11-07 : MGL. First version.
; 
; 
; 
;-
pro red_fits_createfile, filename, head, lun, fileassoc

  ;; Open the new fits file
  openw, lun, filename, /get_lun, /swap_if_little_endian

  ;; Make byte-version of header and write it to the file.
  Nlines = n_elements(head)     ; Number of header lines
  bhdr = replicate(32B, 80L*Nlines)
  for n = 0L, Nlines-1 do bhdr[80*n] = byte(head[n])
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

  ;; Set up an assoc variable for writing the slices
  fileassoc = assoc(lun, array_structure, offset)

  ;; Pad the file, including the data part, to mod 2880 bytes
  endpos = offset + Nelements*abs(bitpix)/8 ; End position of data
  Npad = 2880 - (endpos mod 2880)
  if npad GT 0 and npad LT 2880 then begin
     bpad = bytarr(Npad)
     rec = assoc(lun,bpad,endpos)
     rec[0] = bpad
  endif

end
