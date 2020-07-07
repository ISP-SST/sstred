; docformat = 'rst'

;+
; Correct various mistakes in the fitscube header, that have been made
; during the development of the pipeline. Set the SOLARNET compliance
; status to 1.
;
; Do it also to the spectral cube, if there is one.
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
;    keywords : in, optional, type=struct
;
;      Keywords and values to be added to the header. Keywords are
;      upcased and checked against a list of SOLARNET-approved
;      keywords.
;
;    no_checksum : in, optional, type=boolean
;
;      Do not calculate DATASUM and CHECKSUM.
;
;    no_spectral_file : in, optional, type=boolean
;
;      Do not process any spectral version of the input file.
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
;   2020-07-06 : MGL. Split from red::fitscube_export.
; 
;-
pro red::fitscube_cleanheader, filename $
                               , help = help $
                               , keywords = keywords $
                               , no_spectral_file = no_spectral_file $
                               , no_checksum = no_checksum $
                               , release_date = RELEASE $
                               , release_comment = RELEASEC
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  ;; Make prpara
  red_make_prpara, prpara, keywords     
  red_make_prpara, prpara, no_spectral_file      
  red_make_prpara, prpara, no_checksum 
  red_make_prpara, prpara, release_date      
  red_make_prpara, prpara, release_comment      
  
  ;; The keywords keyword will be checked against this list.
  red_append, approved_keywords, 'AUTHOR'
  red_append, help_keywords, 'Who designed the observation'

  red_append, approved_keywords, 'CAMPAIGN'
  red_append, help_keywords, 'Coordinated campaign name/number, including instance number'

  red_append, approved_keywords, 'DATATAGS'
  red_append, help_keywords, 'Any additional search terms that do not fit in other keywords'

  red_append, approved_keywords, 'OBS_MODE'
  red_append, help_keywords, 'A string uniquely identifying the mode of operation'

  red_append, approved_keywords, 'OBSERVER'
  red_append, help_keywords, 'Who acquired the data'

  red_append, approved_keywords, 'PLANNER'
  red_append, help_keywords, 'Observation planner(s)'

  red_append, approved_keywords, 'PROJECT'
  red_append, help_keywords, 'Name(s) of the project(s) affiliated with the data'

  red_append, approved_keywords, 'REQUESTR'
  red_append, help_keywords, 'Who requested this particular observation'

  red_append, approved_keywords, 'SETTINGS'
  red_append, help_keywords, 'Other settings - numerical values "parameter1=n, parameter2=m"'

  red_append, approved_keywords, 'TELCONFG'
  red_append, help_keywords, 'Telescope configuration'

  
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

  ;; This will be used for a new DATE header.
  new_DATE = red_timestamp(/utc,/iso)

  indir  = file_dirname(filename)+'/'
  infile = file_basename(filename)

;  hdr = red_readhead(filename)

  ;; Read header and fix PRSTEP info 
  red_fitscube_correct_prstep, filename, header = hdr

  ;; Any spectral file to copy?
  if keyword_set(no_spectral_file) then begin
    do_spectral = 0
  endif else begin
    spfile = red_strreplace(infile, '_im.fits', '_sp.fits')
    do_spectral = spfile ne infile and file_test(indir + spfile)
  endelse

  ;; Delete header keywords that should not be there
  red_fitsdelkeyword, hdr, 'STATE'

  ;; Add keywords after this one:
  anchor = 'TIMESYS'

  ;; Update the DATE keyword
  red_fitsaddkeyword, anchor = anchor, hdr, 'DATE', new_DATE

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

  
  ;; After this, this FITS file should be SOLARNET compliant
  red_fitsaddkeyword, hdr, 'SOLARNET', 1, 'SOLARNET compliant file'

  ;; Add header info about this step
  prstep = 'HEADER-CORRECTION'
  self -> headerinfo_addstep, hdr $
                              , prstep = prstep $
                              , prpara = prpara $
                              , prproc = inam

  if do_spectral then begin
    ;; Change some keywords in the spectral file
;    sphdr = headfits(indir+spfile)
    ;; Read header and fix PRSTEP info 
    red_fitscube_correct_prstep, indir+spfile, header = sphdr
    ;; Delete header keywords that should not be there
    red_fitsdelkeyword, sphdr, 'STATE'
    ;; After this, this FITS file should be SOLARNET compliant
    red_fitsaddkeyword, sphdr, 'SOLARNET', 1, 'SOLARNET compliant file'
    ;; Add keywords after this one:
    spanchor = 'TIMESYS'
    self -> headerinfo_addstep, sphdr $
                                , prstep = prstep $
                                , prpara = prpara $
                                , prproc = inam
    red_fitsaddkeyword, anchor = spanchor, sphdr, 'DATE', new_DATE
  endif


;  if n_elements(smooth_width) gt 0 then begin
;    ;; Blur the frames in the fitscube(s)
;    naxis = fxpar(hdr, 'NAXIS*')
;    Nframes = round(product(naxis[2:*]))
;
;    ;; The statistics will not be right after this!
;    
;    for iframe = 0, Nframes-1 do begin
;
;      red_progressbar, iframe, Nframes, /predict, 'Blurring frames'
;
;      red_fitscube_getframe, outdir+outfile, frame, iframe = iframe
;      red_fitscube_addframe, outdir+outfile, smooth(frame, smooth_width), iframe = iframe
;
;    endfor                      ; iframe
;    
;  endif

  if ~keyword_set(no_checksum) then begin

    ;; Don't modify the file(s) after adding datasum and checksum!
    
    ;; Checksums
    datasum = red_fitscube_datasum(filename)

    red_fitsaddkeyword, anchor = anchor, hdr, 'DATASUM', datasum
    fits_add_checksum, hdr
    red_fitscube_newheader, filename, hdr

    if do_spectral then begin
      red_fitsaddkeyword, anchor = spanchor, sphdr, 'DATASUM', datasum
      fits_add_checksum, sphdr
    endif
  endif

  if do_spectral then begin
    ;; Write the new header to the spectral file
    red_fitscube_newheader, indir+spfile, sphdr
  endif

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



a -> fitscube_cleanheader, filename  $
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

h = headfits(filename)
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

d = red_readdata(filename)
test = fits_test_checksum(h, d, ERRMSG = errmsg)
print, 'Test checksums: ', test
if test eq -1 then print, errmsg
end
