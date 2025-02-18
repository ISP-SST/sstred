; docformat = 'rst'

;+
; Extract states information from an array of strings (typically file
; names). 
;
; Replaces the various versions of red_getstates. Based on regular
; expressions and vectorization.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    strings : in, type=strarr
;   
;        A list of strings from which to extract the states
;        information.
;   
;    states : out, optional, type=array(struct)
;
;        An array of structs, containing (partially filled) state
;        information.
; 
; 
; :Keywords:
;
;     force : in, optional, type=boolean
;
;        Do not use cached states.
;
;     strip_settings : in, optional, type=boolean
;
;        Exclude exposure/gain information from the fullstate entries.
; 
; 
; :History:
; 
;   2014-01-22 : First version.
; 
;   2014-04-?? : MGL. Added keyword wavelength.
;
;   2015-09-03 : MGL. Added the "blue" keyword, will now return
;                something meaningful also for blue tilt filter data.
;                Bugfix in qw regular expression.
;
;   2016-05-19 : THI. Partial copy to the crisp class. Modify the
;                state structures and keywords for clarity.
;
;   2016-05-24 : MGL. Removed polarization stuff. If the strings are
;                filenames, then look in the headers. Added keywords
;                gain and exposure, return this information if wanted.
;                Make the keyword fullstate set the other keywords.
;
;   2016-05-25 : MGL. Do not assume camera gain is integer. Get some
;                info not in the header from the file names for now.
;                Move comments to where they are needed.
;
;   2016-05-27 : MGL. Get more information from the headers.
;
;   2016-05-30 : MGL. Trim some whitespace. Strip trailing dots in
;                fullstate_list when there is no tuning. Set
;                scannumber and framenumber if fullstate is set. 
;
;   2016-05-31 : MGL. Detect whether filter1 keyword exists.
;
;   2016-05-31 : JLF. Begin using red_keytab to keep track of changing
;		 SOLARNET keywords. 
;
;   2016-06-01 : MGL. Get prefilter wavelengths from header. Remove
;                the boolean keywords that select info to be returned,
;                just return all info possible. Also the basename
;                boolean keyword doesn't seem to do anything. Rewrite
;                to work directly with the struct and not intermediate
;                arrays.
;
;   2016-06-01 : THI. Added gain & exposure to fullstate. Added
;                keywords strip_wb (to exclude tuning from fullstate
;                for WB cameras) and strip_settings (to exclude gain
;                and exposure from fullstate for pinholes).  
;
;   2016-06-02 : MGL. Added progress printout. 
;
;   2016-06-03 : MGL. Filter headers silently.
;
;   2016-06-09 : MGL. Trim filter string. Build the fullstate string
;                so it never starts with an underscore.
;
;   2016-06-09 : JLF. Bugfix. cam_settings couldn't handle exposure times
;		 f.o.m 10 msec, produced ****ms_G??.??. 
;
;   2016-06-09 : MGL. Change exposure time format to handle longer
;                exposures while preserving format for short
;                exposures. 
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-09-06 : MGL. Call red_meta2head rather than
;                red_filterchromisheaders. 
;
;   2016-09-19 : MGL. Let the subprogram find out its own name. Added
;                conversion from HRE digital units to wavelength
;                tuning.
;
;    2016-09-21 : MGL. Work in meters, not Å or mÅ. Typos hzr --> hrz.
;
;    2016-09-22 : MGL. Make exception for darks when looking for
;                 tuning info.
;
;    2016-09-27 : MGL. Detect is_wb also if the camera is not already
;                 known. 
;
;    2016-09-28 : MGL. Use TEXPOSUR for cam_settings of summed files.
;                 Parse STATE correctly. FULLSTATE for darks without
;                 tuning info. Do not complain if calling from
;                 hrz_zeropoint.
;
;    2016-10-13 : MGL. Added nframes state keyword.  
;
;    2016-10-27 : MGL. Put the camera into the states if found.
;
;    2017-05-08 : MGL. Match camera name to get WB status.
; 
;    2017-06-19 : MGL. Use red_fitspar_getwavelnth.
; 
;    2017-12-01 : MGL. Cache states during IDL session. New keyword
;                 force.
;
;    2019-07-23 : OA. Renamed to chromis::exctractstates_nodb, from
;                 now on chromis::extractstates is a wrapper around db
;                 and nondb versions.
;
;    2022-11-05 : MGL. Recoded to handle new (from 2022-11-03) files
;                 where filenames and headers have proper filter and
;                 tuning info. Older files handled by
;                 chromis::extractstates_wheelhrz_nondb.
;
;-
pro chromis::extractstates_nondb, strings, states $
                                  , force = force $
                                  , polcal = polcal $
                                  , strip_settings = strip_settings
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)            

  if self.isodate lt red_dates(tag = 'CHROMIS tuning metadata') then begin
    ;; Old data with wheel and hrz instead of prefilter and tuning. 
    self -> extractstates_wheelhrz_nondb, strings, states $
       , force = force $
       , strip_settings = strip_settings
  endif

  strings = strtrim(strings,2)
  idx = where( strings ne '' )
  if min(idx) ge 0 then strings = strings[ idx ] $
  else return
  Nstrings = n_elements(strings)
  if( Nstrings eq 0 ) then return

  ;; Create array of structs to holed the state information
  states = replicate( {CHROMIS_STATE}, Nstrings )
  
  ;; Read headers and extract information. This should perhaps return
  ;; an array the length of the number of frames rather than the
  ;; number of files?

  for ifile = 0, Nstrings-1 do begin

    if keyword_set(force) then $
       cnt = 0 $
    else $
       this_cache = rdx_cacheget(strings[ifile], count = cnt)
    
    if cnt gt 0 then begin

      red_progressbar, ifile, Nstrings, 'Extract state info from cache', /predict

      states[ifile] = this_cache.state
      
    endif else begin

      red_progressbar, ifile, Nstrings, 'Extract state info from file headers', /predict

      ;; Get the header
      status = -1
      if file_test(strings[ifile]) then begin
        head = red_readhead(strings[ifile], /silent, status=status)
        states[ifile].filename = strings[ifile]
      endif else begin
        print,'file does not exist: ', strings[ifile]
        print,'state information will be incomplete!'
      endelse
      if status ne 0 then begin
        mkhdr, head, ''         ; create a dummy header
        head = red_meta2head(head, metadata = {filename:strings[ifile]})
      endif
      
      ;; Detector
      detector = fxpar(head, red_keytab('detector'), count=count)
      if count gt 0 then states[ifile].detector = strtrim(detector, 2)

      ;; Camera
      camera = fxpar(head, red_keytab('camera'), count=count)
      if count gt 0 then states[ifile].camera = strtrim(camera,2)

      ;; WB or NB?
      states[ifile].is_wb = strmatch(states[ifile].camera,'*-[DW]')             

      ;; Some numerical keywords
      naxis3 = fxpar(head, 'NAXIS3', count=hasnframes)
      if hasnframes then states[ifile].nframes = naxis3 else states[ifile].nframes = 1
      states[ifile].scannumber = fxpar(head, red_keytab('scannumber'))
      states[ifile].framenumber = fxpar(head, red_keytab('framenumber'))
      states[ifile].gain = fxpar(head, 'DETGAIN', count=hasgain) 
      states[ifile].exposure = fxpar(head, 'XPOSURE', count=hasexp)
      texposur = fxpar(head, 'TEXPOSUR', count=hastexp)
      
      ;; State
      state = fxpar(head, 'STATE', count=hasstate)
      if hasstate gt 0 then state_split = strsplit( state, '_',  /extr )
      
      if states[ifile].is_wb then begin
        ;; WB
        red_fitspar_getwavelnth, head, wavelnth = wavelnth, haswav = haswav
        if haswav ne 0 then begin
          states[ifile].pf_wavelength = float(wavelnth)
        endif
        if hasstate then begin
          if n_elements(state_split) eq 1 then begin
            ;; WB flats are not collected together with NB flats,
            ;; hence no real fpi_state.
            states[ifile].fpi_state = state
          endif else begin
            states[ifile].fpi_state = strjoin(state_split[1:2], '_')
          endelse 
          states[ifile].prefilter = state_split[0]
          states[ifile].tuning = state_split[0]+'+0'    
        endif
      endif else begin
        ;; NB
        if hasstate then begin
          states[ifile].fpi_state = strjoin(state_split[1:2], '_')
          states[ifile].prefilter = state_split[1]
          states[ifile].pf_wavelength = double(states[ifile].prefilter)*1e-10
          states[ifile].tuning = strjoin(state_split[1:2], '_')
        endif
      endelse
      
      states[ifile].tun_wavelength = double(strmid(states[ifile].tuning, 0, 4))*1d-10 $
                                     + double(strmid(states[ifile].tuning, 5))*1d-13

      ;; Camera settings
      if hastexp then begin
        ;; This is a summed file, use the single-exposure exposure
        ;; time for the camera setting.
        states[ifile].cam_settings = strtrim(string(texposur*1000 $
                                                    , format = '(f9.2)'), 2) + 'ms'
      end else if hasexp gt 0 then begin
        ;; This is not a summed file.
        states[ifile].cam_settings = strtrim(string(states[ifile].exposure*1000 $
                                                    , format = '(f9.2)'), 2) + 'ms'
      endif
      if hasgain gt 0 then begin
        if hasexp gt 0 then states[ifile].cam_settings += '_'
        states[ifile].cam_settings += 'G' + string(states[ifile].gain, format = '(f05.2)')
      endif

      if states[ifile].tuning eq '0000_+0' then states[ifile].tuning = ''
      
      ;; Remove zero padding
      split_tuning = strsplit(states[ifile].tuning, '_', /extract)
      if n_elements(split_tuning) eq 2 then begin
        states[ifile].tuning = split_tuning[0] + '_' $
                               + strmid(split_tuning[1],0,1) + strtrim(round(abs(split_tuning[1])),2)
      endif
      split_tuning = strsplit(states[ifile].fpi_state, '_', /extract)
      if n_elements(split_tuning) eq 2 then begin
        states[ifile].fpi_state = split_tuning[0] + '_' $
                                  + strmid(split_tuning[1],0,1) + strtrim(round(abs(split_tuning[1])),2)
      endif  
      
      ;; The fullstate string
      undefine, fullstate_list
      if ~keyword_set(strip_settings) then red_append, fullstate_list, states[ifile].cam_settings
      if states[ifile].prefilter ne '' then red_append, fullstate_list, states[ifile].prefilter
      if states[ifile].tuning ne '' then $             
          red_append, fullstate_list, states[ifile].tuning
      states[ifile].fullstate = strjoin(fullstate_list, '_')

      ;; Store in cache
      rdx_cache, strings[ifile], { state:states[ifile] }

    endelse
    
  endfor                        ; ifile
  
end

cd, '/scratch/mats/test_chromisnames/CHROMIS/'

a = chromisred('config.txt', /dev, /no)

dir = '/storage/Incoming/2022.11.03/CHROMIS-data/09:05:35/' ; Regular data
dir = '/storage/Incoming/2022.11.03/CHROMIS-mosaic/08:56:23/' ; Mosaic data
dirN = dir+'Chromis-N/'
dirW = dir+'Chromis-W/'
fnamesW = file_search(dirW+'*_00000_*fits', count = NfilesW)
fnamesN = file_search(dirN+'*_00000_*fits', count = NfilesN)

a -> extractstates, fnamesN, statesN, /force
a -> extractstates, fnamesW, statesW, /force

help, statesN[0]
help, statesW[0]

end

filesN = file_search('/storage/Incoming/2022.11.03/CHROMIS-data/09:05:35/Chromis-N/*_00001_*.fits')
filesW = file_search('/storage/Incoming/2022.11.03/CHROMIS-data/09:05:35/Chromis-W/*_00001*.fits')
a -> extractstates, filesN, statesN
a -> extractstates, filesW, statesW

stop

;; Test caching

dirN = '/storage/sand02/Incoming/2016.09.11/CHROMIS-flats/*/Chromis-N/'
dirW = '/storage/sand02/Incoming/2016.09.11/CHROMIS-flats/*/Chromis-W/'
fnamesW = file_search(dirW+'*fits', count = NfilesW, /force)
fnamesN = file_search(dirN+'*fits', count = NfilesN, /force)

tic
a -> extractstates, fnamesN, statesN1
toc

tic
a -> extractstates, fnamesW, statesW
toc

tic
a -> extractstates, fnamesN, statesN1
toc

tic
a -> extractstates, fnamesN, statesN2, /force
toc


tic
a -> extractstates, fnamesW, statesW
toc

stop

;; Test darks
files = file_search('darks/cam*.dark', count = Nfiles)
a -> extractstates, files, states

stop

;; Test flats
files = file_search('flats/camXXX_*.flat', count = Nfiles)
a -> extractstates, files, states


stop

;; Test narrowband
dirN = '/storage/sand02/Incoming/2016.09.11/CHROMIS-flats/*/Chromis-N/'
fnamesN = file_search(dirN+'*fits', count = NfilesN)
if NfilesN gt 0 then a -> extractstates, fnamesN, statesN

stop

;; Test wideband
dirW = '/storage/sand02/Incoming/2016.09.11/CHROMIS-flats/*/Chromis-W/'
fnamesW = file_search(dirW+'*fits', count = NfilesW)
if NfilesW gt 0 then a -> extractstates, fnamesW, statesW

end
