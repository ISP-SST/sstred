; docformat = 'rst'

;+
; Return a fits-type header from an ANA fz format file.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    fname : in, type=string
;
;       The name of the data file.
;
; :Keywords:
;
;    date_beg : out, type=strarr
;
;       The timestamps for exposure start.
; 
;    framenumbers : out, type=intarr
;
;       The framenumbers extracted from the file metadata.
;    
; 
; :History:
; 
;    2017-03-10 : MGL. Moved reading of headers from ANA fz format
;                 files from red_readhead.pro.
; 
;    2018-06-04 : MGL. New keywords date_beg and framenumbers.
; 
;-
function red_readhead_ana, fname, $
                           date_beg = date_beg, $
                           framenumbers = framenumbers

  compile_opt idl2

  anaheader = fzhead(fname)

  if n_elements(anaheader) ne 0 then begin

    if strmatch( anaheader, "SIMPLE*" ) gt 0 then begin

      ;; It's already a fits-type header, split into strarr with length 80

      len = strlen(anaheader)
      for i=0,len-1, 80 do begin
        card = strmid(anaheader,i,80)
	if (nc=strlen(card)) lt 80 then card+=string(replicate(32b,80-nc))
        red_append, header, card
        if strmid( card, 0, 2 ) eq 'END' then break
      endfor                    ; i

      framenumbers = fxpar(header, 'FRAMENUM')

    endif else begin

      ;; Do what we can with a "raw" fz file header.

      ;; Get size info from the file
      openr, llun, fname, /get_lun
      ;; Read the first few bytes to get data type and
      ;; number of dimensions:
      ah = bytarr(192) 
      readu, llun, ah
      case ah[7] of
        0: dtyp = 1               ; ANA_BYTE
        1: dtyp = 2               ; ANA_WORD
        2: dtyp = 3               ; ANA_LONG
        3: dtyp = 4               ; ANA_FLOAT
        4: dtyp = 5               ; ANA_DOUBLE
        else: dtyp = 2            ; default?
      endcase

      ;; Read bytes 192-255 as longs to get the dimensions
      bh = lonarr(16)
      readu, llun, bh
      naxisx = bh[0:ah[8]-1]

      ;; Close the file
      free_lun, llun

      ;; Get what other info is possible from the anaheader:
      header = red_anahdr2fits(anaheader, datatype = dtyp, naxisx = naxisx)

      if arg_present(framenumbers) then begin
        lastpart = (strsplit(fname, '[._]', /extract))[-1]
        if strlen(lastpart) eq 7 then framenumbers = long(lastpart)
      endif
      
    endelse

    date_beg = fxpar(header, 'DATE-BEG')

  endif

  return, header

end

