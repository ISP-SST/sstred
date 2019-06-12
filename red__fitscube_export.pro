; docformat = 'rst'

;+
; Make an adjusted version of a fitscube that is suitable for
; exporting to an SVO.
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
;    keywords: in, optional, type=struct
;
;      Keywords and values to be added to the header. Keywords are
;      checked agains a list of SOLARNET-approved keywords.
; 
;    outdir: in, optional, type=string, default="Get from filename"
;   
;      The directory in which to store the new file.
;
;    overwrite: in, optional, type=boolean
;
;      Set this to overwrite an existing file.
;
; :History:
; 
;   2019-06-12 : MGL. First version.
; 
; 
; 
; 
;-
pro red::fitscube_export, filename $
                          , keywords = keywords $
                          , outdir = outdir $
                          , overwrite = overwrite $
                          , release_date = RELEASE $
                          , release_comment = RELEASEC 
  
  inam = red_subprogram(/low, calling = inam1)

  indir  = file_dirname(filename)+'/'
  infile = file_basename(filename)

  ;; The keywords keyword will be checked against this list.
  approved_keywords = ['AUTHOR'      $ ; Who designed the observation
                       , 'DATATAGS'  $ ; Any additional search terms that do not fit in other keywords.
                       , 'CAMPAIGN'  $ ; Coordinated campaign name/number, including instance number
                       , 'OBSERVER'  $ ; Who acquired the data.
                       , 'OBS_MODE'  $ ; A string uniquely identifying the mode of operation.
                       , 'PLANNER'   $ ; Observation planner(s).
                       , 'PROJECT'   $ ; Name(s) of the project(s) affiliated with the data
                       , 'REQUESTR'  $ ; Who requested this particular observation.
                       , 'SETTINGS'  $ ; Other settings - numerical values “parameter1=n, parameter2=m”.
                       , 'TELCONFG'  $ ; Telescope configuration.
                      ]

  
  if n_elements(outdir) eq 0 then outdir = indir

  outfile = red_strreplace(infile, '_im.fits', '_export_im.fits')

  if file_same(filename, outdir+outfile) then begin
    print, inam + ' : This operation would overwrite the input file.'
    return
  endif
  
  if file_test(outdir+'/'+outfile) and ~keyword_set(overwrite) then begin
    print, inam + ' : The output file already exists: '
    spawn, 'ls -l '+outdir+'/'+outfile
    return
  endif
  
  hdr = red_readhead(filename)
  anchor = 'DATE'

  ;; Delete header keywords that should not be there
  red_fitsdelkeyword, hdr, 'STATE'
  
  
  for itag = 0, n_tags(keywords)-1 do begin
    ;; Add FITS header keywords as requested
    if max(strmatch(approved_keywords, (tag_names(keywords))[itag])) then begin
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
  red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASE', RELEASE
  red_fitsaddkeyword, anchor = anchor, hdr, 'RELEASEC', RELEASEC



  ;; Add FITS keywords with info available in the file but not as
  ;; headers. 

  ;; WCS info: DATE-BEG and DATE-END take care of the temporal
  ;; coordinates, WAVEMIN and WAVEMAX of the spectral coordinates. But
  ;; we need to specify the pointing as keywords.
  red_fitscube_getwcs, filename, coordinates = coordinates   
  red_fitsaddkeyword, anchor = anchor, hdr, 'SVO_HPLN' $
                      , mean(coordinates.hpln), 'Average HPLN coordinate'
  red_fitsaddkeyword, anchor = anchor, hdr, 'SVO_HPLT' $
                      , mean(coordinates.hplt), 'Average HPLT coordinate'
  ;; We should also indicate whether the cube is polarimetric
  red_fitsaddkeyword, anchor = anchor, hdr, 'SVO_STOK' $
                      , boolean(fxpar(hdr,'NAXIS4') gt 1), 'Is this a Stokes cube?'
  ;; Check in SVO if data from other telescopes have kwywords for this!

  file_copy, filename, outdir+outfile
  fxhmodify, outdir+outfile, new_header = hdr

  ;; Copy also the spectral file? Or just offer to serve it together
  ;; with the search result?
  
end

a = crispred(/dev)

filename = '/scratch/mats/2016.09.19/CRISP-aftersummer/cubes_nb/nb_6302_2016-09-19T09:30:20_scans=12-16_stokes_corrected_im.fits'

a -> fitscube_export, filename $
                      , keywords = { observer:'Some person' $
                                     ,  requestr:'Boss person' $
                                   } $
                      , release_date = '2012-09-19'

end
