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
;    The image contained in the file fname.
; 
; :Params:
; 
;    fname : in, type=string
;
;       The file name.
;
; :Keywords:
;
;    header : out, type=string array/struct
;
;  	Output the header. 
; 
;    filetype : in, type=string
;
;	The type of file to read. Allowed values are: 
;		ptgrey-fits
;		fz
;	If not set auto-detection will be attempted.
;
;    structheader : in, type=flag
;
;	If set output the header as a structure.
;
;    status : out, type=signed int
;
;	Output the return status, 0 for success, -1 for failure.
;
; :History:
; 
;   2016-05-16 : MGL. First version
; 
;   2016-05-17 : MGL. Default for filetype based on file name.
;
;   2016-05-18 : JLF. Returns the image (like readfits). Added Status
;                flag keyword.
;
;   2016-05-18 : JLF. Use red_anahdr2fits to create FITS headers from
;                ANA headers.
;
;   2016-05-23 : JLF. Use red_filterchromisheaders to clean pt-grey
;                headers. Use rdx_filetype.
;
;   2016-05-26 : THI. Allow reading of ana headers which are actually
;                fits-headers. 
;
;-
function red_readdata, fname $
                       , header = header $
                       , filetype = filetype $
                       , structheader = structheader $
                       , status = status

  if file_test(fname) eq 0 then begin

    message, 'File does not exist: '+fname,/info
    status = -1
    return, 0B

  endif

  if n_elements(filetype) eq 0 then begin

    filetype = rdx_filetype(fname)
  
    if filetype eq '' then begin
      message,'Cannot detect filetype. Pass it manually as',/info
      message,"img = red_readdata('"+fname+"',filetype='fits')",/info
      status = -1
      return, 0B
    endif
        
  endif                         ; filetype


  case strupcase(filetype) of

     'ANA' : begin
        
        ;; Data stored in ANA fz format files.
        
        fzread, data, fname, anaheader

        if n_elements(anaheader) ne 0 then begin
           if strmatch( anaheader, "SIMPLE*" ) gt 0 then begin
               ;; it's actually a fits-header, split into strarr with length 80
               len = strlen(anaheader)
               for i=0,len-1, 80 do begin
                   card = strmid(anaheader,i,80)
                   red_append, tmpheader, card
                   if strmid( card, 0, 2 ) eq 'END' then break
               endfor
               header = tmpheader
           endif else begin
               ;; Convert ana header to fits header
               header = red_anahdr2fits( anaheader, img=data )
           endelse
        endif

     end

     'FITS' : begin
        red_rdfits, fname, header = header
        ;; Data stored in fits files, 
        bit_shift = 0
        if fxpar(header, 'SOLARNET') eq 0 then begin
            caminfo = red_camerainfo( red_camtag(fname) )
            if strmatch(caminfo.model,'PointGrey*') then begin     ; hack to load weird PointGrey data
                uint = 1
                swap = 0
                bit_shift = -4
            endif
        endif
        red_rdfits, fname, image = data, uint=uint, swap=swap
        if( bit_shift ne 0 ) then data = ishft(data, bit_shift)
     end

  endcase

  
  if n_elements(header) ne 0 then $
    header = red_filterchromisheaders(header,meta={filename:fname})
  
  if n_elements(header) ne 0 and keyword_set(structheader) then begin
     header = red_paramstostruct(header)
  endif
  
  status = 0
  
  return, data

end
