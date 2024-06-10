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
;    dark : in, optional, type=array
;
;       The dark frame associated with the momfbd file to be read.
;       Needed if the keyword rawstatistics is used. 
;
;    direction : in, optional, type=integer, default=0
;
;       Apply red_rotate(...,direction) to the data.
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
;    rawstatistics : in, optional, type=boolean
;
;       Calculate statistics for the raw data that went into the
;       restored image in the momfbd file to be read.
;
;    structheader : in, type=flag
;
;	If set output the header as a structure.
;
;    status : out, type=signed int
;
;	Output the return status, 0 for success, -1 for failure.
;
;    nslice : in, optional, type=integer
;
;       Read just this frame from a data cube. (Implemented only for
;       FITS files for now, and that needs some work, too.)
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
;   2017-11-10 : MGL. Use the new /crop option with red_mozaic to get
;                predictable dimensions.
;
;   2018-06-04 : MGL. Renamed keyword framenumber to nslice. Support
;                compressed fits files.
;
;   2020-01-31 : MGL. New keywords rawstatistics and dark.
;
;   2020-03-31 : MGL. New keyword direction.
;
;-
function red_readdata, fname $
                       , dark = dark $
                       , extension = extension $
                       , direction = direction $
                       , filetype = filetype $
                       , header = header $
                       , nslice = nslice $
                       , rawstatistics = rawstatistics $
                       , silent = silent $
                       , status = status $
                       , structheader = structheader $
                       , tabkey = tabkey

  compile_opt idl2
  
  if file_test(fname) eq 0 then begin

    message, 'File does not exist: '+fname, /info
    status = -1
    return, 0B

  endif

  if n_elements(filetype) eq 0 then begin
    filetype = rdx_filetype(fname)
    if filetype eq '' then begin
      message,'Cannot detect filetype. Pass it manually as, e.g.,',/info
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

      ;; Read the original header
      header = headfits(fname, errmsg = errmsg, silent = silent)
      if errmsg ne '' then begin

        message, errmsg + ' : '+fname, /info
        status = -1
        return, 0B

      endif

      
      if fxpar(header, 'SOLARNET') ne 0 then begin
        ;; This is a SOLARNET file 

        ;; Read the data and get the header that is consistent with
        ;; the actual data (e.g., if compressed).
        data = rdx_readdata(fname, status = status)
        header = rdx_readhead(fname, status = status $
                              , date_beg = date_beg, framenumber = framenumber)

        if n_elements(nslice) ne 0 then begin
          
          dims = size(data, /dim)
          Ndims = n_elements(dims)
          if Ndims gt 3 then data = reform(data, dims[0], dims[1] $
                                           , round(product(dims[2:*])), /overwrite)
          data = data[*, *, Nslice]

          ;; Change the header so it is consistent with the selected
          ;; slice
          for iaxis = 1, fxpar(header, 'NAXIS') do sxdelpar, header, 'NAXIS'+strtrim(iaxis, 2)
          fxaddpar, header, 'NAXIS', n_elements(dims),'NAXIS of slice'
          for iaxis = 1, fxpar(header, 'NAXIS') do fxaddpar, header, 'NAXIS'+strtrim(iaxis, 2), dims[iaxis-1]
          ;; And set the FRAMENUM keyword.
          fxaddpar, header, 'FRAMENUM', fxpar(header, 'FRAMENUM', comment = comment)+Nslice $
                    , 'Frame number of slice '+strtrim(Nslice, 2)
          fxaddpar, header, 'DATE-BEG', date_beg[nslice] $
                    , 'DATE-BEG of slice '+strtrim(Nslice, 2)

          date = (strsplit(date_beg[nslice], 'T', /extract))[0]
          time_beg = red_time2double((strsplit(date_beg[nslice], 'T', /extract))[1])
          fxaddpar, header, 'DATE-END' $
                    , date + 'T' + red_timestring(time_beg + fxpar(header, 'XPOSURE')) $
                    , 'DATE-END of slice '+strtrim(Nslice, 2)
          fxaddpar, header, 'DATE-AVG' $
                    , date + 'T' + red_timestring(time_beg + fxpar(header, 'XPOSURE')/2d) $
                    , 'DATE-AVG of slice '+strtrim(Nslice, 2)
          
          status = 0
        endif

;        return, data

      endif else begin          

        ;; Now take care of non-SOLARNET data files.
        
        ;; Data stored in fits files, but what kind?
        header = headfits(fname)
;      red_rdfits, fname, header = header
        bit_shift = 0

        caminfo = red_camerainfo( red_detectorname(fname,head=header) )
        if strmatch(caminfo.model,'PointGrey*') then begin 
          ;; This is the first version PointGrey data from
          ;; spring 2016. Hack to load it:
          uint = 1
          swap = 0
          bit_shift = -4
        endif

        ;; Now read the data
        
        ;; primary HDU
        if n_elements(extension) eq 0 then begin
          if keyword_set(uint) then begin
            ;; readfits does not support uint data so use Pit's
            ;; red_rdfits instead. 
            red_rdfits, fname, image = data $
                        , uint=uint, swap=swap, framenumber = nslice
            ;; Compensate for initial weird format
            if bit_shift ne 0 then data = ishft(data, bit_shift)
          endif else begin
            ;; By default use readfits so we can read 4-dimensional
            ;; cubes.
            data = readfits(fname, Nslice = nslice, silent = silent)
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
      endelse       
    end

    'MOMFBD' : begin
      if arg_present(rawstatistics) then begin
        mr = momfbd_read(fname)
        if max(strmatch(tag_names(mr),'NAME')) eq 1 then begin
          ;; There is raw file info in the momfbd file!
          Nfiles = n_elements(mr.name)
          if min(file_test(mr.name)) eq 1 then begin

            ;; The raw files all exist!
            if n_elements(dark) eq 0 then dark = 0.
            undefine, rawmedians
            for ifile = 0, Nfiles-1 do begin
              rawims = red_readdata(mr.name[ifile], h = rawhdr)
              rawdims = size(rawims, /dim)
              rawframenos = fxpar(rawhdr, 'FRAMENUM')+indgen(rawdims[2])
              for iframe = 0, rawdims[2]-1 do red_append, rawmedians, median(rawims[*, *, iframe]-dark)
            endfor              ; ifile
            ;; Return only medians for now!
            rawstatistics = {medians:rawmedians, framenumbers:rawframenos}
          endif
        endif
      endif else begin
        mr = momfbd_read(fname, /img)
      endelse
      if n_elements(mr.patch) eq 1 then begin
        ;; If just one patch, return it as is
        data = mr.patch.img
      endif else begin
        ;; If more than one patch, make mosaic
        data = red_mozaic(mr, /crop)
      endelse
    end

  endcase
  
  if arg_present(header) then $
     header = red_readhead(fname, structheader = structheader, $
                           silent=silent,$
                           extension=extension)

  status = 0

  if n_elements(direction) gt 0 then begin
    return, red_rotate(data, direction)
  endif else begin
    return, data
  endelse
  
end 


;; Test compressed file
fname = '/nadir-scratch/tomas/wfwfs/ffov_12_compressed/sst_camXXXVI_00000_0000000.fits'
d = red_readdata(fname, h = h, status = status)
d = red_readdata(fname, h = h, status = status, nslice = 11)

;; Test fitsheader file
fname = '/scratch/mats/2016.09.19/CHROMIS-jaime_recipe/momfbd/09:28:36/3950/cfg/results/camXXX_2016-09-19T09:28:36_00120_12.00ms_G10.00_3934_3934_-430.fitsheader'
dd = red_readdata(fname, h = hd, status = status)

date_beg = red_fitsgetkeyword(fname, 'DATE-BEG', variable_values = date_begs, count = count)

stop

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
