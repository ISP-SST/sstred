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
;    feature : in, optional, type=string
;
;      The value of FEATURE keyword, names of solar features in the FOV.
;
;    observer  : in, optional, type=string
;
;      The value of the OBSERVER keyword, names of observers.
;
;    point_id : in, optional, type=string
;
;      The value of POINT_ID keyword, typically date and timestamp,
;      sometimes with '_mosaic' or '_grouped' endings.
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
;   2020-07-15 : MGL. Recursively finalize also the matching WB cube.
;
;   2022-03-08 : OA. New keywords feature, observer, point_id
;
;   2022-04-05 : OA. Moved part of the code to 'red__fitscube_header_finalize'.
; 
;-
pro red::fitscube_finalize, filename $
                            , do_spectral_file = do_spectral_file $
                            , header = hdr $
                            , help = help $
                            , keywords = keywords $
                            , no_checksum = no_checksum $
                            , no_wbcube = no_wbcube $
                            , no_write = no_write $
                            , old_header = oldhdr $
                            , old_spectral_header = oldsphdr $
                            , old_wb_header = oldwbhdr $
                            , release_comment =  release_comment $
                            , release_date = release_date $
                            , spectral_header = sphdr $
                            , wb_header = wbhdr $
                            , feature = feature $
                            , observer = observer $
                            , silent = silent $
                            , point_id = point_id
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; We could look for PRPROCn eq 'red::fitscube_finalize' in the
  ;; header and if we find it show the PRREFn date and ask if we want
  ;; to finalize again.  
  
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

  if keyword_set(point_id) then pp = point_id
  self -> fitscube_header_finalize, hdr $
                       , keywords = keywords $
                       , no_checksum = no_checksum $
                       , coordinates = wcs $
                       , release_date = release_date $
                       , release_comment = release_comment $
                       , feature = feature $
                       , observer = observer $
                       , silent = silent $
                       , point_id = pp    

  ;; Any spectral file to copy?
  if keyword_set(do_spectral_file) then begin
    spfile = red_strreplace(infile, '_im.fits', '_sp.fits')
    do_spectral = spfile ne infile and file_test(indir + spfile)
  endif else begin
    do_spectral = 0
  end
  
;specfile:
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
    ;if do_spectral then red_fitscube_checksums, indir+spfile
    ;; Don't modify the file(s) after adding datasum and checksum!
  endif
  
  if ~keyword_set(no_wbcube) then begin
    ;; Get wb file name from the last PRREF keyword with an align
    ;; reference. (If this is a NB cube!)
    indx = where(strmatch(hdr, '*Align reference:*'), Nwhere)
    if Nwhere gt 0 then begin
      key = strtrim((strsplit(hdr[indx[-1]], '=', /extract))[0], 2)
      wcfile = red_strreplace(fxpar(hdr, key), 'Align reference: ', '')
      if file_dirname(wcfile) eq '.' then wcfile = 'cubes_wb/'+wcfile
      ;; Recursively call fitscube_finalize with the wb cube name and
      ;; many of the origial keywords.
      self -> fitscube_finalize, wcfile $
                                 , header = wbhdr $
                                 , keywords = keywords $
                                 , no_checksum = no_checksum $
                                 , /no_wbcube $
                                 , no_write = no_write $
                                 , old_header = oldwbhdr $
                                 , release_comment = release_comment $
                                 , release_date = release_date $
                                 , feature = feature $
                                 , observer = observer $
                                 , silent = silent $
                                 , point_id = point_id
    endif
  endif

end


a = crispred(/dev)

nname = 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=2,3_stokes_corrected_im.fits'


h = headfits(nname)
hgrep, h, 'Align reference:'

stop

a -> fitscube_finalize, nname

h = headfits(nname)
hgrep, h, 'Align reference:'



stop

wname = 'cubes_wb/wb_6302_2016-09-19T09:30:20_scans=0-43_corrected_im.fits'

h1 = headfits(wname)

a -> fitscube_finalize, wname

h2 = headfits(wname)

stop

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
