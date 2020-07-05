; docformat = 'rst'

;+
; Make an adjusted version of a fitscube that is suitable for
; exporting to an SVO, and copy it to its archival directory.
;
; Copy also the spectral cube, if there is one.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Params:
; 
;    filename : in, type=string
; 
;      The path of the fitscube file. 
; 
; 
; :Keywords:
;
;    smooth_width: in, optional, type=integer
;  
;      If present, blur all images in the exported cube with
;      smooth(frame, blur). This is to facilitate making sample cubes
;      without actually releasing data. If unity, interpreted as a
;      boolean and set to a default width of 11 pixels.
;
;    keywords: in, optional, type=struct
;
;      Keywords and values to be added to the header. Keywords are
;      upcased and checked against a list of SOLARNET-approved
;      keywords.
;
;    no_spectral_file: in, optional, type=boolean
;
;      Do not copy the sp version of the input file, if there is one.
;
;    outdir: in, out, optional, type=string, default='/storage_new/science_data/YYYY-MM-DD/'
;   
;      The directory in which to store the new file/to which the file
;      was exported.
;
;    outfile: out, optional, type=string
;   
;      The name of the file to which the data were exported.
;
;    overwrite: in, optional, type=boolean
;
;      Set this to overwrite an existing file.
;
;    release_date : in, optional, type=string
;
;      The value of the RELEASE keyword, a date after which the data
;      are not proprietary.
;
;    release_comment : in, optional, type=string
;
;      The value of the RELEASEC keyword, a comment with details about
;      the proprietary state of the data, and whom to contact.
;
;    svo_api_key : in, optional, type=string
;
;      Ingest the metadata of the cube into the prototype SVO using
;      this API_KEY.
;
;    svo_username : in, optional, type=string
;
;      Ingest the metadata of the cube into the prototype SVO using
;      this USERNAME.
;
; :History:
; 
;   2019-06-12 : MGL. First version.
; 
;   2019-06-17 : MGL. New keyword help.
; 
;   2019-06-18 : MGL. New keywords smooth_width and outfile. Add
;                prstep header info.
; 
;   2019-07-04 : MGL. Optionally (opt-out) add checksum and datasum.
; 
;   2019-10-18 : MGL. New keywords svo_api_key and svo_username.
; 
;-
pro red::fitscube_export, filename $
                          , help = help $
                          , keywords = keywords $
                          , no_spectral_file = no_spectral_file $
                          , no_checksum = no_checksum $
                          , outdir = outdir $
                          , outfile = outfile $
                          , overwrite = overwrite $
                          , release_date = RELEASE $
                          , release_comment = RELEASEC $
                          , smooth_width = smooth_width $
                          , svo_api_key = svo_api_key $
                          , svo_username = svo_username
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  ;; Make prpara
  red_make_prpara, prpara, keywords     
  red_make_prpara, prpara, no_spectral_file      
  red_make_prpara, prpara, release_date      
  red_make_prpara, prpara, release_comment      
  red_make_prpara, prpara, smooth_width
  
  ;; The keywords keyword will be checked against this list.
  approved_keywords = strarr(10)
  help_keywords     = strarr(10)
  approved_keywords[0] = ['AUTHOR']   & help_keywords[0] = 'Who designed the observation'
  approved_keywords[1] = ['DATATAGS'] & help_keywords[1] = 'Any additional search terms that do not fit in other keywords'
  approved_keywords[2] = ['CAMPAIGN'] & help_keywords[2] = 'Coordinated campaign name/number, including instance number'
  approved_keywords[3] = ['OBSERVER'] & help_keywords[3] = 'Who acquired the data'
  approved_keywords[4] = ['OBS_MODE'] & help_keywords[4] = 'A string uniquely identifying the mode of operation'
  approved_keywords[5] = ['PLANNER']  & help_keywords[5] = 'Observation planner(s)'
  approved_keywords[6] = ['PROJECT']  & help_keywords[6] = 'Name(s) of the project(s) affiliated with the data'
  approved_keywords[7] = ['REQUESTR'] & help_keywords[7] = 'Who requested this particular observation'
  approved_keywords[8] = ['SETTINGS'] & help_keywords[8] = 'Other settings - numerical values "parameter1=n, parameter2=m"'
  approved_keywords[9] = ['TELCONFG'] & help_keywords[9] = 'Telescope configuration'

  if n_elements(smooth_width) gt 0 then begin
    ;; We need to implement reading/writing frames from/to spectral
    ;; fitscubes for this to work easily for spectral cubes.
    no_spectral_file = 1
    ;; If smooth_width is set to the inconsequential value of unity,
    ;; interpret it as a boolean and set to a default width of 11.
    if smooth_width eq 1 then smooth_width = 11
  endif
  

  if keyword_set(help) then begin
    print
    print, inam + ' : These are the approved keywords:'
    print
    for i = 0, n_elements(approved_keywords)-1 do begin
      st = '         : ' + help_keywords[i]
      strput, st, approved_keywords[i], 0
      print, st
    endfor
    print
    return
  endif

  ;; This will be used for a new DATE header and also to indicate the
  ;; version of the cube as part of the file name.
  new_DATE = red_timestamp(/utc,/iso)

  indir  = file_dirname(filename)+'/'
  infile = file_basename(filename)

  hdr = red_readhead(filename)
  
  date_beg = fxpar(hdr,'DATE-BEG')
  if n_elements(outdir) eq 0 then outdir = '/storage_new/science_data/' $
                                           + (strsplit(date_beg,'T',/extract))[0] + '/' $
                                           + strtrim(fxpar(hdr,'INSTRUME'),2) + '/'

  ;; Any old versions of this file?
  oldfiles = file_search(outdir+'/'+red_strreplace(infile, '_im.fits', '_*.fits'), count = Noldfiles)
  if Noldfiles gt 0 then begin
    s = ''
    print, inam + ' : There are '+strtrim(Noldfiles, 2)+' older versions of this file in '+outdir+':'
    print, file_basename(oldfiles), format = '(a0)'
    read, '    Delete them [N]? ', s
    if s eq '' then s = 'N'
    if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
      file_delete, oldfiles
      ;; Should perhaps also add code that removes the deleted files'
      ;; metadata from the SVO?
    endif
  endif

  if n_elements(smooth_width) gt 0 then begin
    outfile = red_strreplace(infile, '_im.fits', '_blurred_export'+new_DATE+'_im.fits')
  endif else begin
    outfile = red_strreplace(infile, '_im.fits', '_export'+new_DATE+'_im.fits')
  endelse

  if file_same(filename, outdir+outfile) then begin
    print, inam + ' : This operation would overwrite the input file.'
    return
  endif
  
  if file_test(outdir+'/'+outfile) and ~keyword_set(overwrite) then begin
    print, inam + ' : The output file already exists: '
    spawn, 'ls -l '+outdir+'/'+outfile
    return
  endif

  ;; Any spectral file to copy?
  if keyword_set(no_spectral_file) then begin
    copy_spectral = 0
  endif else begin
    spfile = red_strreplace(infile, '_im.fits', '_sp.fits')
    copy_spectral = file_test(indir + spfile)
  endelse

  ;; Delete header keywords that should not be there
  red_fitsdelkeyword, hdr, 'STATE'

  ;; Add keywords after this one:
  anchor = 'TIMESYS'

  ;; Update the DATE keyword
  red_fitsaddkeyword, anchor = anchor, hdr, 'DATE', new_DATE
  red_fitsaddkeyword, anchor = anchor, hdr, 'FILENAME', outfile


  for itag = 0, n_tags(keywords)-1 do begin
    ;; Add FITS header keywords as requested
    if max(strmatch(approved_keywords, strupcase((tag_names(keywords))[itag]))) then begin
      red_fitsaddkeyword, anchor = anchor, hdr $
                          , strupcase((tag_names(keywords))[itag]), keywords.(itag)
    endif else begin
      print, inam + ' : The tag '+(tag_names(keywords))[itag]+' is not in the list of approved keywords: ' 
      print, approved_keywords
      return
    endelse
  endfor                        ; itag
  
  ;; Proprietary data?
  if n_elements(RELEASE) eq 0 then begin
    s = ''
    print, inam + ' : No release date given. Continue without? [Y]'
    read, s
    if s eq '' then s = 'Y'
    if strupcase(strmid(s, 0, 1)) ne 'Y' then return
    RELEASE = ''
    RELEASEC = 'Data is not proprietary.'
  endif else begin
    ;; Check that RELEASE is a proper date
    if ~(strmatch(RELEASE,'[12][90][0-9][0-9]-[012][0-9]-[0123][0-9]') $
         || RELEASE eq '') then begin
      print, inam + ' : Release date is not a proper ISO date: '+RELEASE
      return
    endif

  endelse
  ;; Proprietary data
  if n_elements(RELEASEC) eq 0 then begin
    releasec = 'Proprietary data, se RELEASE keyword for release date.'
  endif
  red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASE',  RELEASE
  red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASEC', RELEASEC



  ;; Add FITS keywords with info available in the file but not as
  ;; headers. 

  ;; WCS info: DATE-BEG and DATE-END take care of the temporal
  ;; coordinates, WAVEMIN and WAVEMAX of the spectral coordinates. But
  ;; we need to specify the pointing as an additional WCS coordinate.
  red_fitscube_getwcs, filename, coordinates = coordinates   
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'WCSNAMEA', 'AVERAGED APPROXIMATE HPLN-TAN/HPLT-TAN CENTER POINT'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRPIX1A', fxpar(hdr,'NAXIS1')/2. + 0.5, 'Center pixel of image array'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRPIX2A', fxpar(hdr,'NAXIS2')/2. + 0.5, 'Center pixel of image array' 
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRVAL1A', mean(coordinates.hpln), '[arcsec] Coordinates of center of image array'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CRVAL2A', mean(coordinates.hplt), '[arcsec] Coordinates of center of image array' 
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CDELT1A', 0.0, 'Zero FOV extent'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CDELT2A', 0.0, 'Zero FOV extent'

  ;; Add header info about this step
  prstep = 'HEADER-CORRECTION'
  if n_elements(smooth_width) then prstep += ', SMOOTHING'
  self -> headerinfo_addstep, hdr $
                              , prstep = prstep $
                              , prpara = prpara $
                              , prproc = inam

  
  ;; Copy the file and write the new header
  print, inam + ' : Copying the fitscube...'
  file_mkdir, outdir
  file_copy, filename, outdir+outfile
;  fxhmodify, outdir+outfile, new_header = hdr
  red_fitscube_newheader, outdir+outfile, hdr

  print, inam + ' : Wrote '+outdir+outfile

  if copy_spectral then begin
    ;; Copy the spectral file as well
    print, inam + ' : Copying the spectral cube...'
    spoutfile = red_strreplace(outfile, '_im.fits', '_sp.fits')
    file_copy, indir+spfile, outdir+spoutfile
    ;; Change some keywords in the spectral file
    sphdr = headfits(outdir+spoutfile)
    ;; Delete header keywords that should not be there
    red_fitsdelkeyword, sphdr, 'STATE'
    ;; Add keywords after this one:
    spanchor = 'TIMESYS'
    red_fitsaddkeyword, anchor = spanchor, sphdr, 'DATE', new_DATE
    red_fitsaddkeyword, anchor = spanchor, sphdr, 'FILENAME', spoutfile
    ;;fxhmodify, outdir+spoutfile, new_header = sphdr
    red_fitscube_newheader, outdir+spoutfile, sphdr
    print, inam + ' : Wrote '+outdir+spoutfile
  endif


  if n_elements(smooth_width) gt 0 then begin
    ;; Blur the frames in the fitscube(s)
    naxis = fxpar(hdr, 'NAXIS*')
    Nframes = round(product(naxis[2:*]))

    ;; The statistics will not be right after this!
    
    for iframe = 0, Nframes-1 do begin

      red_progressbar, iframe, Nframes, /predict, 'Blurring frames'

      red_fitscube_getframe, outdir+outfile, frame, iframe = iframe
      red_fitscube_addframe, outdir+outfile, smooth(frame, smooth_width), iframe = iframe

    endfor                      ; iframe
    
  endif

  if ~keyword_set(no_checksum) then begin

    ;; Don't modify the file(s) after adding datasum and checksum!
    
    ;; Checksums
    datasum = red_fitscube_datasum(outdir+outfile)

    red_fitsaddkeyword, anchor = anchor, hdr, 'DATASUM', datasum
    fits_add_checksum, hdr
    red_fitscube_newheader, outdir+outfile, hdr

    if copy_spectral then begin
      red_fitsaddkeyword, anchor = spanchor, sphdr, 'DATASUM', datasum
      fits_add_checksum, sphdr
      red_fitscube_newheader, outdir+spoutfile, sphdr
    endif
  endif
  
  if n_elements(svo_username) gt 0 and n_elements(svo_api_key) gt 0 then begin

    ;; Ingest metadata into the SVO database.
    
    ;; Look for the python script
    cmd = file_search(strsplit(!path,':',/extract)+'/submit_record.py', count=cnt)
    if cnt eq 0 then begin

      print, inam + " : Didn't find submit_record.py in !path"
      stop
      
    endif else begin

      ;; Keywords to the SVO ingestion script
      file_url = 'https://dubshen.isf.astro.su.se/'
      file_path = (strsplit(date_beg,'T',/extract))[0] + '/' $
                + strtrim(fxpar(hdr,'INSTRUME'),2) + '/'

      ;; Get thumbnail from keyword or generate one from the cube?
      ;; Sharpest from brightest tuning?
      
      ;; Generate OID = unique observation ID of the metadata?
      ;; oid=
      ;; Is date + timestamp enough or should this be a hash-like
      ;; string that changes with versioning etc?

      ;; Build the command string
      cmd = cmd[0]
      cmd += ' ' + strlowcase(strtrim(fxpar(hdr,'INSTRUME'),2)) ; The dataset ID in the SOLARNET Data Archive/SVO
      cmd += ' ' + outdir+outfile                               ; The FITS file to submit to the SVO
      cmd += ' --file-url ' + file_url                          ; The URL of the file
      cmd += ' --file-path ' + file_path                        ; The relative path of the file
      if n_elements(thumbnail_url) gt 0 then $                  ;
         cmd += ' --thumbnail-url ' + thumbnail_url             ; The URL of the thumbnail
      cmd += ' --username ' + svo_username ; The SVO username of the user owning the data
      cmd += ' --api-key ' + svo_api_key   ; The SVO API key of the user owning the data
;      if n_elements(oid) gt 0 then $                ;
;      cmd += ' --oid ' +        ; The unique observation ID of the metadata
      
      ;; Spawn running the script
      spawn, cmd, status
      print, cmd
      
    endelse
  endif else begin

    print, inam + " : Keywords svo_username or svo_api_key not provided."
    print, inam + " : Will not ingest metadata into the SVO database."

  endelse

  
end

;; Code for testing. Test only with smaller cubes because we have to
;; read the entire cube before calling fits_test_checksum below.

if 0 then begin
  cd, '/scratch/mats/2016.09.19/CRISP-aftersummer/'
  filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=12-16_stokes_corrected_im.fits'
  a = crispred(/dev)
endif else begin
  cd, '/scratch/mats/2016.09.19/CHROMIS-jan19/'
  filename = 'cubes_nb/nb_3950_2016-09-19T09:28:36_scans=69-75_corrected_im.fits'
  a = chromisred(/dev)
endelse



a -> fitscube_export, filename, outdir = outdir, outfile = outfile $
                      , /smooth $
                      , keywords = { observer:'Some person' $
                                     , requestr:'Boss person' $
                                   } $
                      , release_comment = 'These are test data, never to be released' $
                      , release_date = '2999-09-19'


search_keywords = [ 'STARTOBS' $ 
                    , 'DATE' $
                    , 'DATE-BEG' $
                    , 'DATE-AVG' $
                    , 'DATE-END' $
                    , 'RELEASE' $
                    , 'RELEASEC' $
                    , 'INSTRUME' $
                    , 'TELESCOP' $
                    , 'WCSNAMEA' $
                    , 'CRVAL1A' $
                    , 'CRVAL2A' $
                    , 'NAXIS1' $
                    , 'NAXIS2' $
                    , 'NAXIS3' $
                    , 'NAXIS4' $
                    , 'NAXIS5' $
                    , 'CTYPE1' $
                    , 'CTYPE2' $
                    , 'CTYPE3' $
                    , 'CTYPE4' $
                    , 'CTYPE5' $
                    , 'WAVEMIN' $
                    , 'WAVEMAX' $
                    , 'WAVEUNIT' $
                  ]

h = headfits(outdir+'/'+outfile)
for i = 0, n_elements(search_keywords)-1 do begin
  keyval = fxpar(h, search_keywords[i], count = n)
  if n gt 0 then begin
    print,search_keywords[i], ' : ', keyval
  endif else begin
    print,search_keywords[i], ' : -- '
  endelse
endfor

print, 'DATASUM:  ',fxpar(h, 'DATASUM')
print, 'CHECKSUM: ',fxpar(h, 'CHECKSUM')

d = red_readdata(outdir+'/'+outfile)
test = fits_test_checksum(h, d, ERRMSG = errmsg)
print, 'Test checksums: ', test
if test eq -1 then print, errmsg
end
