; docformat = 'rst'

;+
;       Read a simple FITS-Formated File into an IDL-Array. Size and
;       scaling of the image are automatically recognized and set by
;       the information of the header.
; 
;       The header is read line-by-line until the END-Statement is
;       found. Size and scaling factor of the image are extracted.
;
;       Only works with Basic FITS files, i.e. the *first line* of the 
;       header must contain the statement
;          SIMPLE =  T
;
;
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;       Pit Sütterlin, 1991-05-04 (as rdfits)
;
; 
; 
; 
; :Returns:
; 
;       The read-in data
;
;
; :Params:
; 
;       file : in, type=string
;
;           The filename of the FITS-File to be read.
;
;
; :Keywords:
; 
;       image : out, optional, type=string
;
;               Variable to return image in.
; 
;       log : in, optional, type=boolean
;
;               If set, the header information is copied to stdout. If
;               LOG is a string, the header gets written to a file
;               with that name. Existing files get appended.
;
;       header : out, optional, type=strarr
;
;               FITS header is passed back.
;
;       swap : in, optional, type=boolean
;
;               Byte-swap the low/high byte of each integer. This
;               overrides autodetection.
;
;       nosign : in, optional, type=boolean
;
;               Data is unsigned integer
;
;       framenumber : in, optional, type=integer
;
;               Read just this frame from a data cube.
; 
; :History:
; 
;       1992-06-03 : PS. Allow for different blocking sizes (BS)
;
;       1992-12-10 : PS. The '92 version of the AT1-Software writes
;                    some zeroes into the FITS-header, which caused
;                    some difficulties with the byte-to-string
;                    conversion. The 0's are replaced by 32 (ASCII for
;                    ' ')
;
;       1992-15-10 : PS. Implement the possibility of reading also
;                    bytearrays. For this, the BITPIX header card must
;                    be set to `8'.
;
;       1994-27-04 : PS. Uncompress zipped files on the fly.
;
;       1995-17-07 : PS. Add some decent error checking for empty
;                    files when reading from exabyte (e.g. double
;                    filemarks -> empty file), new keyword n_retry
;
;       1995-22-08 : PS. New keyword NOFILEMARKS for reading tapefiles
;                    written without separating filemarks. Also added
;                    keyword hd_blocks, which is only needed for
;                    reading tapes w/o filemarks.
;
;       1998-06-08 : PS. If a file is not found, look also for file.gz
;
;       1999-07-01 : PS. Byteswap was wrong for 32bit data
;
;       2000-03-30 : PS. Add bzip2 compressed files
;
;       2001-12-03 : PS. add support for float/double type and UINT
;                    keyword for bad DOT files
;
;       2001-01-11 : PS Auto-swap based on machine endianess (assuming
;                    FITS standard big endian data). Can be overridden
;                    by swap=0/1.
;
;       2004-08-24 : PS. Error in header length computation if header
;                    exactly filled n*2880 bytes
;
;       2016-05-16 : MGL. Rename: red_rdfits. Changed format of
;                    documentation header. Removed code that has to do
;                    with reading tapes and streams. Remove keyword
;                    noread, make image an optional keyword instead.
;                    Look for .Z files if needed. Additional cosmetic
;                    changes.
;
;       2016-05-30 : MGL. New keyword framenumber.
; 
;-
PRO red_rdfits, file $
                , image = image $
                , log = log $
                , swap = swap $
                , header = header $
                , nosign = nosign $
                , uint = ui $
                , framenumber = framenumber

  on_error, 2

  hdr = bytarr(80)

  count = 0
  zipped = 0
  cd, '.', current = work_dir
  endian = (byte(1, 0, 1))[0]   ;;; system endianess: 1-little  2-big
                                ;;; FITS requires data in big endianess

  IF n_params() EQ 1 THEN BEGIN 
     files = findfile(file)
     sf = size(files)
     IF sf[0] EQ 0 THEN BEGIN
         ;;; see if it is compressed: try with extension
        files = findfile(file+'.gz '+file+'.bz2 '+file+'.Z')
        sf = size(files)
        IF sf[0] EQ 0 THEN message, "File doesn't exist !!"
     ENDIF
     IF sf[1] GT 1 THEN BEGIN
        message, 'Ambiguous filename. Possible Expansions are :', /cont
        print, files
        return
     ENDIF
     IF ((strpos(files, '.gz'))[0] GE 0) OR ((strpos(files, '.Z'))[0] GE 0) THEN $
        BEGIN
        origfile = file
        file = nnumber(abs(fix(1e5*randomn(undef))), 5)
        file = '/tmp/tmp'+file+'.fits'
        spawn, 'cat '+files +' | gunzip > '+file
        zipped = 1
     ENDIF ELSE IF ((strpos(files, '.bz2'))[0] GE 0) THEN BEGIN
        origfile = file
        file = nnumber(abs(fix(1e5*randomn(undef))), 5)
        file = '/tmp/tmp'+file+'.fits'
        spawn, 'cat '+files +' | bunzip2 > '+file
        zipped = 1
     ENDIF
  ENDIF


  openr, unit, file, /get_lun


  IF keyword_set(log) THEN BEGIN
     
     s_log = size(log)
     
     IF s_log(s_log[0]+1) EQ 7 THEN $
        openw, logunit, log, /append, /get $
     ELSE $
        logunit=-1
  ENDIF

  sx = 1 & sy = 1 & sz = 1

  header = bytarr(80*200)

  REPEAT BEGIN
     readu, unit, hdr
     header(80*count:80*count+79) = hdr
     count = count+1
     hd1 = string(hdr(0:77))
     IF keyword_set(log) THEN printf, logunit, hd1
     CASE strmid(hd1, 0, 6) OF
        'BITPIX'   :  BEGIN 
           bitpix = fix(strmid(hd1, 10, 21))
           CASE bitpix OF
                8  : refpix = 8b
               16  : IF keyword_set(ui) THEN refpix = uint(16) $
                                        ELSE refpix = fix(16)
               32  : refpix = 32L
              -32  : refpix = float(32)
              -64  : refpix = 64d
              ELSE : message, 'red_rdfits : No support for ' + strtrim(bitpix, 2) + $
                              ' bytes per pixel'
           ENDCASE
        END
        'NAXIS '   :  npic = fix(strmid(hd1, 10, 21))
        'NAXIS1'   :  sx = fix(strmid(hd1, 10, 21))
        'NAXIS2'   :  sy = fix(strmid(hd1, 10, 21))
        'NAXIS3'   :  sz = fix(strmid(hd1, 10, 21))
        'BSCALE'   :  bsc = float(strmid(hd1, 10, 21))
        'BZERO '   :  bze = float(strmid(hd1, 10, 21))
        ELSE       :
     ENDCASE
  ENDREP UNTIL strpos(hd1, 'END   ') EQ 0
  
  IF keyword_set(log) THEN IF logunit GT 0 THEN free_lun, logunit

;-----
;  Read until end of block
;-----
  nhead = count*80
  header = header[0:nhead-1]
  ix = where(header EQ 0, nix)
  IF nix GT 0 THEN header[ix] = 32
  tmp = string(header)
  header = strarr(count)
  FOR i = 0, count-1 DO header[i] = strmid(tmp, 80*i, 80)
  iblocks = ceil(nhead/2880.)
  dif = fix(2880.*iblocks-nhead)
  IF dif GT 0 THEN BEGIN
     tmp = bytarr(dif)
     readu, unit, tmp
  ENDIF

;-----
;  read Image
;-----
  if arg_present(image) then begin

     if n_elements(framenumber) ne 0 then begin
        ;; We want only a single frame from the data cube. 
        if n_elements(sz) eq 0 then sz = 1 ; Cube with only one frame
        if framenumber ge sz then stop     ; The framenumber is outside the cube
        image = replicate(refpix, sx, sy)  ; Variable to read image into
        if framenumber ne 0 then begin     ; Skip over unwanted data
           byte_elem = abs(bitpix)/8 
           Nskip = long64(framenumber)*sx*sy*byte_elem
           mrd_skip, unit, Nskip
        endif
        sxaddpar, header, 'NAXIS', 2, /savecomment ; 2D output
        sxdelpar, header, 'NAXIS3'                 ; 2D output
     endif else begin
        CASE npic OF
           1 : image = replicate(refpix, sx)
           2 : image = replicate(refpix, sx, sy)
           3 : image = replicate(refpix, sx, sy, sz)
           ELSE: stop
        ENDCASE
     endelse

     readu, unit, image

;-----
; Swap High- and Lowbyte?
;            If not explicitely forced, assume data is big endian and
;            swap if the system is not
;-----

     IF n_elements(swap) EQ 0 THEN swap = endian

     IF swap NE 0 THEN BEGIN
        CASE refpix OF
           8: 
           16: byteorder, image, /sswap
           32: byteorder, image, /lswap
           64: byteorder, image, /l64swap
        ENDCASE
     ENDIF

;-----
; Daten als vorzeichenlose Integer?
;-----
     IF keyword_set(nosign) THEN BEGIN
        image = long(image)
        image = image+32768
     ENDIF

;-----
;  scaling of the pixels?
;-----
     IF keyword_set(bsc) THEN IF bsc NE 1 THEN $
        image=image*bsc
     IF keyword_set(bze) THEN IF bze NE 0 THEN $
        image = image+bze

  endif                         ; Read the image?

  ;; Close the file
  free_lun, unit

;-----
;  delete file, if temporarily uncompressed
;-----
  IF zipped THEN BEGIN
     spawn, 'rm '+file
     file = origfile
  ENDIF


END

