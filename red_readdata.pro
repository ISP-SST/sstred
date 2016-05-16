; docformat = 'rst'

;+
; Camera-independent read routine.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Mats Löfdahl, ISP, 2016-05-16
; 
; 
; :Returns:
; 
;    1 for success, 0 for failure.
; 
; :Params:
; 
;    fname : in, type=string
;
;       The file name.
;
; :Keywords:
; 
;    data : out, optional
;   
;       The data.
;
;    header : out, optional, type="strarr or struct"
; 
;       The header.
;
; :History:
; 
;   2016-05-16 : MGL. First version
;
;-
function red_readdata, fname, data = data, header = header, filetype = filetype, structheader = structheader

  if n_elements(filetype) eq 0 then begin

     ;; Try to guess based on file name

  endif

  if file_test(fname) eq 0 then begin

     print, 'red_readdata : File does not exist:'
     print, fname
     return, 0B

  endif

  case filetype of

     'fz' : begin
        
        ;; Data stored in ANA fz format files.
        
        if arg_present(data) then begin
           fzread, data, fname, anaheader
        endif else begin
           anaheader = fzhead(fname)
        endelse

        if n_elements(anaheader) ne 0 then begin
           
           ;; Convert ana header to fits header

           header = strarr(80)
           ih = 0

           ;; Yes, this is simple.
           header[ih] = 'SIMPLE  = T' & ih += 1
           
           ;; Data type
           header[ih] = 'BITPIX = 16' & ih += 1

           ;; Dimensions 
           header[ih] = 'NAXIS   = 2' & ih += 1
           NAXIS1 = strtrim(long(strmid(anaheader, strpos(anaheader, ' W=')+3)), 2)
           header[ih] = 'NAXIS1  = ' + NAXIS1
           NAXIS2 = strtrim(long(strmid(anaheader, strpos(anaheader, ' H=')+3)), 2)
           header[ih] = 'NAXIS2  = ' + NAXIS2

           ;; Origin etc.
           ;;header[ih] = 'ORIGIN  = Institute for Solar Physics' & ih += 1
           ;;header[ih] = 'TELESCOP= Swedish Solar Telescope' & ih += 1

           ;; Time info
           Ts = strmid(anaheader, strpos(anaheader, 'Ts=')+3, 26)
           Te = strmid(anaheader, strpos(anaheader, 'Te=')+3, 26)
           exptime = red_time2double(strmid(Te, strpos(Te, ' '))) $
                     - red_time2double(strmid(Ts, strpos(Ts, ' ')))
           
           ;;header[ih] = 'DATE_OBS= ' + ;; We could get this from the
           ;;file name, if it includes the original path.
           header[ih] = 'DATE    = ' + red_strreplace(Ts, ' ', 'T') & ih += 1
           header[ih] = 'EXPTIME = ' + strtrim(exptime, 2)+ ' / [s]' & ih += 1
           ;;header[ih] = 'INTERVAL= ' & ih += 1
           
           ;; A CRISP image file header continues with 

           ;; End!
           header[ih] = 'END' & ih += 1

        endif

     end

     'ptgrey-fits' : begin

        ;; Data stored in fits files, but in the strange format
        ;; returned by the PointGrey cameras.

        if arg_present(data) then begin
           red_rdfits, fname, image = data, header = header, /uint, swap=0
           data = ishft(data, -4) ; 12 bits in 16-bit variables
        endif else begin
           red_rdfits, fname, header = header
        endelse

     end

  endcase

  if n_elements(header) ne 0 and keyword_set(structheader) then begin
     header = red_paramstostruct(header)
  endif

  return, 1B

end
