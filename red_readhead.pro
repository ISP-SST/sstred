; docformat = 'rst'

;+
; Read headers from data files in a camera independent way.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     J. Lewis Fox, ISP, 2016-05-18
; 
; 
; :Returns:
; 
;     A FITS compatible header (from a data file of any instrument)
; 
; :Params:
; 
;    fname : in, type=string
;
;       The data file.
;
; :Keywords:
;    
;    filetype : in, type=string
;
;	The type of file to read. Allowed values are: 
;		fits
;		fz
;               momfbd
;	If not set auto-detection will be attempted.
;
;    select_frame : in, optional, type=integer
;
;       Specify to get header corresponding to a single frame in a
;       multi-frame image file.
;
;    date_beg : out, type=strarr
;
;       The timestamps for exposure start.
; 
;    framenumbers : out, type=intarr
;
;       The framenumbers extracted from the file metadata.
; 
;    noexternal : in, optional, type=boolean
;
;       Set this to not read any external FITS header.
;
;    structheader : in, type=flag
;
;	If set output the header as a structure.
;
;    status : out, type=signed int
;
;	Output the return status, 0 for success, -1 for failure.
;
;    extension : in, type=string
;
;	The name of the extension (EXTNAME) that we are interested in.
;
; :History:
; 
;   2016-05-18 : JLF. Created.
;
;   2016-05-31 : JLF. Added keyword silent to suppress informational
;                messages. 
;
;   2016-08-12 : MGL. Make header for momfbd-format file. New keyword
;                framenumber. Call red_filterchromisheaders only for
;                FITS files without SOLARNET header keyword. Call new
;                function red_meta2head.
;
;   2016-08-19 : MGL. Only test for .momfbd extension if filetype is
;                not specified.      
;
;   2016-08-23 : MGL. Deal with ANA headers lacking data size
;                information.      
;
;   2016-08-31 : MGL. Use red_detectorname instead of red_getcamtag.
;                Add rudimentary handling of tabulated DATE-BEG.       
;
;   2016-09-01 : MGL. Change FRAME1 to FRAMENUM if needed.       
;
;   2016-09-08 : MGL. Allow for headers without prefilter info (like
;                darks).        
;
;   2016-09-15 : MGL. Fix bug in construction of DATE-END from table.
;                Fix bug in removal of empty lines. Also construct
;                DATE-AVE. Construct DATE-BEG, DATE-END, and DATE-AVE
;                only if they do not already exist. ANA files: get
;                dimensions and data type from the first 256 bytes of
;                the file.
;
;   2016-09-21 : MGL. Change filter tags to four-digit tags
;                representing approximate filter wavelength.
;
;   2016-10-13 : MGL. Removed unnecessary check for W and H in ANA
;                headers.
;
;   2016-10-19 : MGL. Look for data under /storage/ rather than under
;                /mnt/ in Stockholm.
;
;   2016-11-05 : JLF. Read in extension headers using the extension
;                keyword. 
;
;   2017-03-10 : MGL. Move file-format specific code to their own
;                subprograms. Trust rdx_filetype() to recognize momfbd
;                files.  
;
;   2017-03-16 : MGL. Look for (and optionally read) an external FITS
;                header.
;
;   2017-09-01 : THI. Get date_beg and framenumbers from file. Rename
;                keyword framenumber to select_frame.
; 
;   2017-12-01 : MGL. Use status from red_readhead_fits.
;
;
;-
function red_readhead, fname, $
                       filetype = filetype, $
                       date_beg = date_beg, $
                       framenumbers = framenumbers, $
                       select_frame = select_frame, $
                       noexternal = noexternal, $
                       structheader = structheader, $
                       status = status, $
                       silent = silent, $
                       extension = extension

  compile_opt idl2
  
  if( file_test(fname) eq 0 ) then begin
    message, 'File does not exist: '+fname, /info
    status = -1
    return, 0B
  endif
  
  if( n_elements(filetype) eq 0 ) then begin
    filetype = rdx_filetype(fname)
    if( filetype eq '' ) then begin
      message, 'Cannot detect filetype. Pass it manually as, e.g.,', /info
      message, "head = red_readhead('"+fname+"',filetype='fits')", /info
      status = -1
      return, 0B
    endif
  endif
  
  ;if( n_elements(filetype) eq 0 ) then filetype = rdx_filetype(fname)
  hdir = file_dirname(fname)
  
  ;; Read the header based on filetype and construct the name of an
  ;; external header file.
  case strupcase(filetype) of
    'ANA' : begin
      header = red_readhead_ana(fname)
      hname = file_basename(fname, '.fz')
      hname = hdir+'/'+file_basename(hname, '.f0') + '.fitsheader'
    end
    'FITS' : begin
      header = red_readhead_fits(fname, $
                                 date_beg = date_beg, $
                                 framenumbers = framenumbers, $
                                 select_frame = select_frame, $
                                 silent = silent, $
                                 status = status, $
                                 extension = extension)
      if status ne 0 then return, 0B
      hname = hdir+'/'+file_basename(fname, '.fits')   + '.fitsheader'
    end
    'MOMFBD' : begin
      header = red_readhead_momfbd(fname, version = momfbd_version)
      hname = hdir+'/'+file_basename(fname, '.momfbd') + '.fitsheader'
    end
    else     : begin
      message, 'Cannot detect filetype. Pass it manually as e.g.', /info
      message, "head = red_readhead('"+fname+"',filetype='fits')", /info
      status = -1
      return, 0B
    end
  endcase

  ;; Keyword FRAME1 changed to FRAMNUM. Rewrite old headers
  ;; to match the new standard.
  frnm = sxpar(header, 'FRAMENUM', count = Npar)
  if Npar eq 0 then begin
    fr1 = sxpar(header, 'FRAME1', count = Npar)
    if Npar eq 1 then sxaddpar, header, 'FRAMENUM', fr1
  endif

  if n_elements(extension) eq 0 then $
     header = red_meta2head(header, meta={filename:fname})
  
  if keyword_set(structheader) then begin
    header = red_paramstostruct(header)
  endif

  ;; Read the external fitsheader if it exists.
  if ~keyword_set(noexternal) && file_test(hname) then begin
    xhead = headfits(hname)
    ;; Add or replace all keywords from the external header except
    ;; DATE, BITPIX, NAXIS*, COMMENT, etc.
    Nlines = n_elements(xhead)
    keys = strtrim(strmid(xhead, 0, 8), 2)
    skipkeys = ['SIMPLE', 'BITPIX', 'NAXIS', 'DATE', 'LONGSTRN', 'CONTINUE', 'COMMENT', 'END']
    naxis = fxpar(xhead, 'NAXIS')
    for i = 1, naxis do red_append, skipkeys, 'NAXIS'+strtrim(i, 2)
    for iline = 0, Nlines-1 do begin
      if ~max(strmatch(skipkeys, keys[iline])) then begin
        value = fxpar(xhead, keys[iline], comment = comment)
        fxaddpar, header, keys[iline], value, comment
        if strmatch(keys[iline], 'PRLIB*') $
           && (strmatch(value, 'momfbd/redux*')) then begin
          if n_elements(momfbd_version) ne 0 then begin
            libindx = red_strreplace(keys[iline], 'PRLIB', '')
            fxaddpar, header, 'PRVER'+libindx, momfbd_version, after = keys[iline] $
                      , 'Library version/MJD of last update (From .momfbd file)'
          endif
        endif
      endif
    endfor                      ; iline
  endif
  
  status = 0

  ;; Remove blank lines
  header = header[where(header ne string(replicate(32B,80)))] 

  return, header

end


;; Test broken file

dir='/storage/sand05n/Incoming/2017.04.20/CHROMIS-flats/18:34:10/Chromis-N/'
fname='sst_camXXX_00004_0036400_wheel00006_hrz34410.fits'

header = red_readhead(dir+fname,date_beg=date_beg,framenumbers=framenumbers, status = status)
print, status

stop

cd,'/storage/sand02/Incoming/2016.09.19/CHROMIS-flats/11:21:22/Chromis-N',current=curdir
fname = 'sst_camXXX_00000_0000000_wheel00005_hrz32061.fits'

head = red_readhead(fname)
thead = red_readhead(fname,extension='tabulations')
ehead = red_readhead(fname,/extension)

stop

cd,'/polar-scratch/mats/2016.09.19/CHROMIS/calib_tseries'
fname = 'wb_3950_2016-09-19T09:28:36_scans=68-72_corrected_im.fits'

head = red_readhead(fname)
thead = red_readhead(fname,extension='tabulations')
ehead = red_readhead(fname,/extension)

stop

cd,curdir

end
