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
;    extension : in, optional, type=string
;
;	Read from a FITS extension with extname taken from the value of 
;	extension keyword. If it is an image extension the image is returned.
;	If it is a binary or ASCII table the column data are returned which 
; 	correspond to the value of tabkey (or by default the 1st column).
;	
;	The extension header is returned in the header keyword (if specified)
;
;    tabkey : in, optional, type=string
;
;	The name of the tabulated keyword to extract from the table.
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
;                Swap endian if needed. Remove byteorder keywords.
;
;   2016-11-04 : JLF. Added extension and tabkey keywords to support 
;		 reading from tab-HDUs and other extensions.      
;
;   2016-12-07 : MGL. Use readfits so we can read 4-dimensional cubes.
;
;
;-
function red_readdata, fname $
                       , header = header $
                       , filetype = filetype $
                       , structheader = structheader $
                       , status = status $
                       , silent = silent $
                       , framenumber = framenumber $
                       , extension = extension $
                       , tabkey = tabkey

  compile_opt idl2
  
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
            caminfo = red_camerainfo( red_detectorname(fname,head=header) )
            if strmatch(caminfo.model,'PointGrey*') then begin 
                ;; This is the first version PointGrey data from
                ;; spring 2016. Hack to load it:
                uint = 1
                swap = 0
                bit_shift = -4
            endif

        endif
        
        ;; Now read the data
        
        ;; primary HDU
        if n_elements(extension) eq 0 then begin
          if keyword_set(uint) then begin
            ;; readfits does not support uint data so use Pit's
            ;; red_rdfits instead. 
            red_rdfits, fname, image = data $
                        , uint=uint, swap=swap, framenumber = framenumber
            ;; Compensate for initial weird format
            if bit_shift ne 0 then data = ishft(data, bit_shift)
          endif else begin
            ;; By default use readfits so we can read 4-dimensional
            ;; cubes.
            data = readfits(fname, Nslice = framenumber)
          endelse

	  ;; Does it need to be byteswapped?
	  doswap = 0
	  endian = sxpar(header,'ENDIAN', count = Nendian)
	  if Nendian eq 1 then if endian eq 'little' then doswap = 1
	  byteordr = sxpar(header,'BYTEORDR', count = Nbyteordr)
	  if Nbyteordr eq 1 then if byteordr eq 'LITTLE_ENDIAN' then doswap = 1

	  if doswap then begin
	    swap_endian_inplace, data
	    sxdelpar, header, 'ENDIAN'
	    sxdelpar, header, 'BYTEORDR'
	  endif 
	endif else begin
	  header = red_readhead(fname,extension=extension,silent=silent)
	  exttype = fxpar(header,'XTENSION')
	  case exttype of 
	    'IMAGE   ': begin
	      fxread,fname,data,extension=extension
	    end
	    'TABLE   ': begin
	      if size(extension,/type) ne 7 then extno = extension
	      fits_read,fname,tab,exten=extno,extname=extension,/no_pdu
	      if n_elements(tabkey) ne 0 then $
		data = ftget(header,tab,tabkey) $
	      else data = ftget(header,tab,1) ; default to first field.
	    end
	    'BINTABLE': begin
	      fxbopen,tlun,fname,extension
	      if n_elements(tabkey) ne 0 then $
		fxbread,tlun,data,tabkey $
	      else fxbread,tlun,data,1 ; default to the first column
	      fxbclose,tlun
	    end
	  endcase
	endelse
	
     end

     'MOMFBD' : begin
        mr = momfbd_read(fname, /img)
        data = red_mozaic(mr, /clip)
     end

  endcase
  
  if arg_present(header) then $
     header = red_readhead(fname, structheader = structheader, $
                           framenumber = framenumber, silent=silent,$
                           extension=extension)

  status = 0
  
  return, data

end


cd,'/storage/sand02/Incoming/2016.09.19/CHROMIS-flats/11:21:22/Chromis-N',current=curdir
fname = 'sst_camXXX_00000_0000000_wheel00005_hrz32061.fits'

img = red_readdata(fname)
tab = red_readdata(fname,head=header,extension='tabulations',tabkey='date-beg')
tabdefault = red_readdata(fname,extension='tabulations')
tabextno = red_readdata(fname,/extension)

stop

cd,'/polar-scratch/mats/2016.09.19/CHROMIS/calib_tseries'
fname = 'wb_3950_2016-09-19T09:28:36_scans=68-72_corrected_im.fits'

img = red_readdata(fname)
tabdefault = red_readdata(fname,extension='tabulations')
tab = red_readdata(fname,extension='tabulations',tabkey='time',header=head)
tabextno = red_readdata(fname,/extension)

stop

cd,curdir

end
