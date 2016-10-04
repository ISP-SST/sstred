; docformat = 'rst'

;+
; Construct filenames (or search strings matching file names) for
; CHROMIS files.
; 
; :Categories:
;
;    CHROMIS pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; :Returns:
; 
;    Filenames or search strings.
; 
; :Params:
; 
;    datatype : in, type=string
; 
;       The type of file. One of 'dark', 'flat', 'sflat', 'gain', or
;       'pinh'. Should implement also momfbd output (of both file
;       types). Offset files or the current equivalent info. Restored
;       data cubes with scans and/or movies. Any other files that are
;       relevant to the pipeline.
;
;    states : in, type="array of structs"
;
;       State info used to generate the file names.
; 
; :Keywords:
; 
;    no_dir : in, optional, type=boolean
; 
;       If set, generate file names without directory. 
; 
;    no_fits : in, optional, type=boolean
; 
;       If set, generate file names without .fits extension. 
; 
;    raw : in, optional, type=boolean
; 
;       If set, generate names for raw data of the specified type.
; 
;    wild_* : in, optional, type="boolean or string"
;   
;       If a string, use this as part of the search string. Otherwise,
;       if true, construct the proper type of wild card.
; 
; :History:
; 
;    2016-09-29 : MGL. First version.
; 
; 
;-
function chromis::filenames, datatype, states $
                             , wild_detector = wild_detector $         ; camXXVII
                             , wild_camera = wild_camera $             ; Chromis-N
                             , wild_fullstate = wild_fullstate $       ; 12.00ms_G10.00_3934_3934_+0
                             , wild_scannumber = wild_scannumber $     ; Five digits
                             , wild_framenumber = wild_framenumber $   ; Seven digits
                             , wild_tuning = wild_tuning $             ; 3934_+0
                             , wild_prefilter = wild_prefilter $       ; Four digits
                             , wild_exposure = wild_exposure $         ; 12.00ms
                             , wild_gain = wild_gain $                 ; G10.00
                             , wild_cam_settings = wild_cam_settings $ ; 12.00ms_G10.00
                             , no_fits = no_fits $
                             , no_dir = no_dir $
                             , raw = raw $
                             , timestamp = timestamp

  Nstates = n_elements(states)
  filenames = strarr(Nstates)
 
  if keyword_set(raw) then begin
     if ~keyword_set(wild_prefilter) and ~keyword_set(wild_tuning) then begin
        print, 'chromis::filenames : Can only generate raw data searchstrings with'
        print, '                     wildcards for prefilter and tuning at this point.'
        return, 0
     endif
  endif


  ;; Shorthand to make several fields wild:
  if keyword_set(wild_fullstate) then begin
     wild_exposure = 1
     wild_gain = 1
     wild_prefilter = 1
     wild_tuning = 1
  endif
  if keyword_set(wild_cam_settings) then begin
     wild_exposure = 1
     wild_gain = 1
  endif

  ;; Search strings, the wild parts.
  detector_searchstring    = 'cam[XVI]*'
  camera_searchstring      = 'Chromis-?' 
  scannumber_searchstring  = strjoin(replicate('[0-9]', 5))
  framenumber_searchstring = strjoin(replicate('[0-9]', 7))
  prefilter_searchstring   = strjoin(replicate('[0-9]', 4))
  exposure_searchstring    = '[0-9]*\.[0-9]*ms' 
  gain_searchstring        = 'G[0-9]*\.[0-9]*'
  tuning_searchstring      = strjoin(replicate('[0-9]', 4)) + '_[+-][0-9]*'
  
  ;; If wild, set to search strings (possibly the one entered with the
  ;; wild_* keyword).
  if keyword_set(wild_detector) then begin
     if strlowcase(size(wild_detector, /tname)) eq 'string' then begin
        detector = wild_detector
     endif else begin
        detector = detector_searchstring
     endelse
  endif
  if keyword_set(wild_camera) then begin
     if strlowcase(size(wild_camera, /tname)) eq 'string' then begin
        camera = wild_camera
     endif else begin
        camera = camera_searchstring
     endelse
  endif
  if keyword_set(wild_scannumber) then begin
     if strlowcase(size(wild_scannumber, /tname)) eq 'string' then begin
        scannumber = wild_scannumber
     endif else begin
        scannumber = scannumber_searchstring
     endelse
  endif
  if keyword_set(wild_framenumber) then begin
     if strlowcase(size(wild_framenumber, /tname)) eq 'string' then begin
        framenumber = wild_framenumber
     endif else begin
        framenumber = framenumber_searchstring
     endelse
  endif
  if keyword_set(wild_prefilter) then begin
     if strlowcase(size(wild_prefilter, /tname)) eq 'string' then begin
        prefilter = wild_prefilter
     endif else begin
        prefilter = prefilter_searchstring
     endelse
  endif
  if keyword_set(wild_exposure) then begin
     if strlowcase(size(wild_exposure, /tname)) eq 'string' then begin
        exposure = wild_exposure
     endif else begin
        exposure = exposure_searchstring 
     endelse
  endif
  if keyword_set(wild_gain) then begin
     if strlowcase(size(wild_gain, /tname)) eq 'string' then begin
        gain = wild_gain
     endif else begin
        gain = gain_searchstring
     endelse
  endif
  if keyword_set(wild_tuning) then begin
     if strlowcase(size(wild_tuning, /tname)) eq 'string' then begin
        tuning = wild_tuning
     endif else begin
        tuning = tuning_searchstring 
     endelse
  endif


  ;; Now generate filenames/searchstrings for the given states.
  for istate = 0, Nstates-1 do begin

     ;; If not wild, set to value from states. 
     if ~keyword_set(wild_detector)    then detector    = states[istate].detector
     if ~keyword_set(wild_camera)      then camera      = states[istate].camera
     if ~keyword_set(wild_scannumber)  then scannumber  = string(states[istate].scannumber, format = '(i05)')
     if ~keyword_set(wild_framenumber) then framenumber = string(states[istate].framenumber, format = '(i07)')
     if ~keyword_set(wild_prefilter)   then prefilter   = states[istate].prefilter
     if ~keyword_set(wild_tuning)      then tuning      = states[istate].tuning
     ;; Get exposure and detector gain from cam_settings because there
     ;; they are already correctly formatted strings.
     cam_settings = strsplit(states[istate].cam_settings, '_', /extract)
     if ~keyword_set(wild_exposure)    then exposure    = cam_settings[0]
     if ~keyword_set(wild_gain)        then gain        = cam_settings[1]

     undefine, tag_list
     dir = ''
     ext = ''

     if keyword_set(raw) then begin

        ;; Raw data

        ;; If wild_prefilter/wild_tuning are only booleans
        if strlowcase(size(wild_prefilter, /tname)) ne 'string' then prefilter = 'w*'
        if strlowcase(size(wild_tuning, /tname)) ne 'string' then tuning = 'hrz*'

        red_append, tag_list, 'sst'
        red_append, tag_list, detector

        case strlowcase(strmid(datatype, 0, 4)) of
           
           'dark' : begin
              dirs = *self.dark_dir
              red_append, tag_list, scannumber
              red_append, tag_list, framenumber
           end

           'flat' : begin
              dirs = *self.flat_dir
              red_append, tag_list, scannumber
              red_append, tag_list, framenumber
              red_append, tag_list, prefilter
              red_append, tag_list, tuning
           end

           'pinh' : begin
              dirs = *self.pinh_dirs
              red_append, tag_list, scannumber
              red_append, tag_list, framenumber
              red_append, tag_list, prefilter
              red_append, tag_list, tuning
           end

           'scie' : begin
              dirs = *self.data_dirs
              red_append, tag_list, scannumber
              red_append, tag_list, framenumber
              red_append, tag_list, prefilter
              red_append, tag_list, tuning
           end

        endcase
        
        if ~keyword_set(no_dir) then begin
           ;; Can we include a directory?
           if n_elements(dirs) ge 1 then begin
              ;; Only if dark_dir is specified
              if n_elements(dirs) eq 1 then begin
                 ;; If only one, use it (possibly with a different timestamp)
                 if n_elements(timestamp) eq 1 then begin
                    dir = file_dirname(dirs[0]) + '/' + timestamp + '/'
                 endif else begin
                    dir = dirs[0] + '/'
                 endelse
              endif else begin
                 tdirs = file_basename(dirs)  
                 bdirs = file_dirname(dirs)
                 ubdirs = bdirs[uniq(bdirs, sort(bdirs))]
                 if n_elements(ubdirs) eq 1 then begin
                    ;; If all dark_dir are the same except for the
                    ;; timestamp, we can do something:
                    if n_elements(timestamp) eq 1 then begin
                       ;; If a timestamp is given in keyword, use it
                       dir = ubdirs[0] + '/' + timestamp + '/'
                    endif else begin
                       dir = ubdirs[0] + '/[0-2][0-9]:[0-5][0-9]:[0-5][0-9]/'
                    endelse
                 endif
              endelse                         ; eq 1
           endif                              ; ge 1
           if dir ne '' then dir += 'Chromis-?/' ; We could set the correct dir based on detector!
        endif

        ext = '.fits'
        
     endif else begin

        ;; Not raw data

        case strlowcase(strmid(datatype, 0, 4)) of
           
           'dark' : begin
              dir = self.out_dir+'/darks/'
              red_append, tag_list, detector
              red_append, tag_list, exposure
              red_append, tag_list, gain
              ext = '.dark'
              if ~keyword_set(no_fits) then ext += '.fits'
           end

           'flat' : begin
              dir = self.out_dir + '/flats/' 
              red_append, tag_list, detector
              red_append, tag_list, exposure
              red_append, tag_list, gain
              red_append, tag_list, prefilter
              if states[istate].is_wb eq 0 and tuning ne '' then $
                 red_append, tag_list, tuning
              ext = '.flat'
              if ~keyword_set(no_fits) then ext += '.fits'
           end

           'sfla' : begin
              dir = self.out_dir + '/flats/' 
              red_append, tag_list, detector
              red_append, tag_list, exposure
              red_append, tag_list, gain
              red_append, tag_list, prefilter
              if states[istate].is_wb eq 0 and tuning ne '' then $
                 red_append, tag_list, tuning
              red_append, tag_list, 'summed'
              ext = '.flat'
              if ~keyword_set(no_fits) then ext += '.fits'
           end

           'gain' : begin
              dir = self.out_dir + '/gaintables/' 
              red_append, tag_list, detector
              red_append, tag_list, exposure
              red_append, tag_list, gain
              red_append, tag_list, prefilter
              if states[istate].is_wb eq 0 and tuning ne '' then $
                 red_append, tag_list, tuning
              ext = '.gain'
              if ~keyword_set(no_fits) then ext += '.fits'
           end

           'pinh' :  begin
              dir = self.out_dir+'/pinhs/'
              red_append, tag_list, detector
              red_append, tag_list, prefilter
              if states[istate].is_wb eq 0 and tuning ne '' then $
                 ;; Should be replaced with the NB tuning info!
                 red_append, tag_list, tuning
              ext = '.pinh'
              if ~keyword_set(no_fits) then ext += '.fits'
           end

        endcase

     endelse

     if keyword_set(no_dir) then begin
        filenames[istate] = strjoin(tag_list, '_') + ext
     endif else begin
        filenames[istate] = dir + strjoin(tag_list, '_') + ext
     endelse

  endfor                        ; istate
  
  return, filenames

end 


a = chromisred()

snames = file_search('flats/camXXX*fits')
;snames = file_search('flats/camXXVII*fits')
a -> extractstates, snames, states

datanames = a -> filenames('science', states, /raw, /wild_tuning, /wild_prefilter, /wild_scan, /wild_frame)


stop
dnames = a -> filenames('dark', states)

dnames = a -> filenames('dark', states, wild_detector = 'camXVI')


fnames = a -> filenames('flat', states, /wild_tuning)
pnames = a -> filenames('pinh', states, /no_fits)
gnames = a -> filenames('gain', states, /no_fits)

end
