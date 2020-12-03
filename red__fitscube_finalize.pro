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
;    do_spectral_file : in, optional, type=boolean
;
;      Do process the spectral version of the input file (if there is
;      one). This should not be needed as crispex uses only the main
;      HDU data, and detects if it is different than in the regular
;      file. 
;
;    keywords : in, optional, type=struct
;
;      Keywords and values to be added to the header. Keywords are
;      upcased and checked against a list of SOLARNET-approved
;      keywords.
;
;    no_write : in, optional, type=boolean
;
;      Don't write the new header(s) (and don't compute checksums).
;      Combine this with keywords header, old_header, spectral_header,
;      and old_spectral_header to check what changes would be done.
;    
;    header : out, optional, type=strarr
;
;      The new header.    
;    
;    old_header : out, optional, type=strarr
;    
;      The old header.    
;    
;    spectral_header : out, optional, type=strarr
;    
;      The new header of the spectral file.    
;    
;    old_spectral_header : out, optional, type=strarr
;
;      The old header of the spectral file.    
;
;    no_checksum : in, optional, type=boolean
;
;      Do not calculate DATASUM and CHECKSUM.
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
;   2020-07-13 : MGL. New keywords no_write, header, old_header,
;                spectral_header, old_spectral_header.
; 
;   2020-07-15 : MGL. Remove keyword no_spectral_file, add keyword
;                do_spectral_file. 
; 
;-
pro red::fitscube_finalize, filename $
                            , do_spectral_file = do_spectral_file $
                            , header = hdr $
                            , old_header = oldhdr $
                            , spectral_header = sphdr $
                            , old_spectral_header = oldsphdr $
                            , help = help $
                            , keywords = keywords $
                            , no_checksum = no_checksum $
                            , no_write = no_write $
                            , release_date = RELEASE $
                            , release_comment = RELEASEC
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  ;; Make prpara
  red_make_prpara, prpara, keywords     
  red_make_prpara, prpara, do_spectral_file      
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
    print, 'You can add/change them by specifying them with a struct like this:'
    print, 'a -> fitscube_finalize, filename, keywords={REQUESTR:"Boss P. Erson", PROJECT:"Spots at disk center"}'
    print
    return
  endif

  ;; This will be used for a new DATE header.
  new_DATE = red_timestamp(/utc,/iso)

  indir  = file_dirname(filename)+'/'
  infile = file_basename(filename)

  ;; To compare with new header later?
  oldhdr = red_readhead(filename)
  
  ;; Read header and fix PRSTEP info 
  red_fitscube_correct_prstep, filename, header = hdr, /nowrite

  ;; Any spectral file to copy?
  if keyword_set(do_spectral_file) then begin
    spfile = red_strreplace(infile, '_im.fits', '_sp.fits')
    do_spectral = spfile ne infile and file_test(indir + spfile)
  endif else begin
    do_spectral = 0
  end

  ;; Delete header keywords that should not be there
  red_fitsdelkeyword, hdr, 'STATE'

  ;; Add keywords after this one:
  red_fitsaddkeyword, anchor = 'EXTNAME', hdr, 'TIMESYS', "UTC"
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
    print, inam + ' : No release date given. Mark as "Not proprietary"? [Yn]'
    read, s
    if s eq '' then s = 'Y'
    if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
      RELEASE = ''
      RELEASEC = 'Data is not proprietary.'
    endif
  endif else begin
    ;; Check that RELEASE is a proper date
    if ~(strmatch(RELEASE,'[12][90][0-9][0-9]-[012][0-9]-[0123][0-9]') $
         || RELEASE eq '') then begin
      print, inam + ' : Release date is not a proper ISO date: '+RELEASE
      return
    endif
  endelse
  if n_elements(RELEASE) gt 0 then begin
    ;; Proprietary data
    if n_elements(RELEASEC) eq 0 then begin
      releasec = 'Proprietary data, see RELEASE keyword for release date.'
    endif
    red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASE',  RELEASE
    red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASEC', RELEASEC
  endif
  


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
  red_fitsaddkeyword, anchor = 'DATE', hdr, 'SOLARNET', 1, 'SOLARNET compliant file'

  ;; Add header info about this step
  prstep = 'HEADER-CORRECTION'
  self -> headerinfo_addstep, hdr $
                              , prstep = prstep $
                              , prpara = prpara $
                              , prproc = inam

  if do_spectral then begin
    ;; Change some keywords in the spectral file
    oldsphdr = headfits(indir+spfile)
    ;; Read header and fix PRSTEP info 
    red_fitscube_correct_prstep, indir+spfile, header = sphdr, /nowrite
    ;; Delete header keywords that should not be there
    red_fitsdelkeyword, sphdr, 'STATE'
    ;; After this, this FITS file should be SOLARNET compliant
    red_fitsaddkeyword, anchor = 'DATE', sphdr, 'SOLARNET', 1, 'SOLARNET compliant file'
    ;; Add keywords after this one:
    red_fitsaddkeyword, anchor = 'EXTNAME', sphdr, 'TIMESYS', "UTC"
    spanchor = 'TIMESYS'
    self -> headerinfo_addstep, sphdr $
                                , prstep = prstep $
                                , prpara = prpara $
                                , prproc = inam
    red_fitsaddkeyword, anchor = spanchor, sphdr, 'DATE', new_DATE
  endif

  if ~keyword_set(no_write) then begin
    red_fitscube_newheader, filename, hdr
    red_fitscube_statistics, filename, /write
    if do_spectral then begin
      ;; Write the new header to the spectral file
      red_fitscube_newheader, indir+spfile, sphdr
    endif
  endif


  if ~keyword_set(no_checksum) and ~keyword_set(no_write) then begin
    red_fitscube_checksums, filename
    if do_spectral then red_fitscube_checksums, spfilename
    ;; Don't modify the file(s) after adding datasum and checksum!
  endif

end

;; Code for testing. Test only with smaller cubes because we have to
;; read the entire cube before calling fits_test_checksum below.

if 1 then begin
  cd, '/scratch/mats/2016.09.19/CRISP-aftersummer/'
  filename = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=2-8_stokes_corrected_im.fits'
  a = crispred(/dev)
endif else begin
  cd, '/scratch/mats/2016.09.19/CHROMIS-jan19/'
  filename = 'cubes_nb/nb_3950_2016-09-19T09:28:36_scans=69-75_corrected_im.fits'
  filename = 'cubes_nb/nb_3950_2016-09-19T10:42:01_scans=0-4_corrected_cmapcorr_im.fits'
  a = chromisred(/dev)
endelse

a -> fitscube_finalize, filename, no_write = 0 $
                        , header = hdr $
                        , old_header = oldhdr $
                        , spectral_header = sphdr $
                        , old_spectral_header = oldsphdr  $
                        , keywords = { observer:'Some person' $
                                       , requestr:'Boss person' $
                                     } $
                        , release_comment = 'These are test data, never to be released' $
                        , release_date = '2999-12-31'




end

a -> fitscube_finalize, filename  $
                        , keywords = { observer:'Some Person' $
                                       , requestr:'Boss Person' $
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
