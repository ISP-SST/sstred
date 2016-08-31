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
;		fits
;		ana
;               momfbd
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
;    framenumber : in, optional, type=integer
;
;       Read just this frame from a data cube. (Implemented only for
;       FITS files for now.)
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
;   2016-05-31 : JLF. Added silent keyword to suppress informational messages.
;
;   2016-06-08 : JLF. Bugfix, don't process headers if the user isn't asking
;		 for them.
;
;   2016-05-30 : MGL. New keyword framenumber.
;
;   2016-08-12 : MGL. Read also momfbd-format files, return the
;                mosaicked image.
;
;   2016-08-15 : MGL. Don't make header to be returned in header
;                keyword, get it from red_readhead instead. Read only
;                the parts of .momfbd files that are needed to make
;                the mosaic.
;
;   2016-08-19 : MGL. Only test for .momfbd extension if filetype is
;                not specified.     
;
;   2016-08-31 : MGL. Use red_detectorname instead of red_getcamtag.
;                Swap endian if needed.
;
;-
function red_readdata, fname $
                       , header = header $
                       , filetype = filetype $
                       , structheader = structheader $
                       , status = status $
                       , silent = silent $
                       , framenumber = framenumber

  if file_test(fname) eq 0 then begin

    message, 'File does not exist: '+fname, /info
    status = -1
    return, 0B

  endif

  if n_elements(filetype) eq 0 then begin

    ;; Remove this when rdx_filetype can recognize .momfbd files:
    if file_basename(fname,'.momfbd') ne file_basename(fname) then begin
       filetype = 'momfbd'
    endif else begin
       filetype = rdx_filetype(fname)
    endelse
     
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

        if arg_present(header) then header = red_readhead(fname)

     end

     'FITS' : begin

        ;; Data stored in fits files, but what kind?
        red_rdfits, fname, header = header
        bit_shift = 0
        if fxpar(header, 'SOLARNET') eq 0 then begin
            caminfo = red_camerainfo( red_detectorname(fname) )
            if strmatch(caminfo.model,'PointGrey*') then begin 
                ;; This is the first version PointGrey data from
                ;; spring 2016. Hack to load it:
                uint = 1
                swap = 0
                bit_shift = -4
            endif

        endif
        
        ;; Now read the data
        red_rdfits, fname, image = data $
                    , uint=uint, swap=swap, framenumber = framenumber
        
        ;; Compensate for initial weird format
        if( bit_shift ne 0 ) then data = ishft(data, bit_shift)

        ;; Does it need to be byteswapped?
        doswap = 0
        endian = sxpar(header,'ENDIAN', count = Nendian)
        if Nendian eq 1 then if edian eq 'little' then doswap = 1
        byteordr = sxpar(header,'BYTEORDR', count = Nbyteordr)
        if Nbyteordr eq 1 then if byteordr eq 'LITTLE_ENDIAN' then doswap = 1

        if doswap then begin
           swap_endian_inplace, data
        endif 

     end

     'MOMFBD' : begin
        mr = momfbd_read(fname, /img)
        data = red_mozaic(mr, /clip)
     end

  endcase
  
  if arg_present(header) then $
     header = red_readhead(fname, structheader = structheader $
                           , framenumber = framenumber)

  status = 0
  
  return, data

end
