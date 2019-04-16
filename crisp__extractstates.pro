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
;      A list of strings from which to extract the states information.
;   
;    states : out, optional, type=array(struct)
;
;        An array of structs, containing (partially filled) state information.
; 
; 
; :Keywords:
; 
;     force : in, optional, type=boolean
;
;        Do not use cached states.
; 
;     polcal : in, optional, type=boolean
; 
;        Set this to add polcal-specific items in the states, qw and
;        lp. 
;
;     strip_wb : in, optional, type=boolean
;
;        Exclude tuning information from the fullstate entries for WB
;        cameras
; 
;     strip_settings : in, optional, type=boolean
;
;        Exclude exposure/gain information from the fullstate entries.
; 
; 
; :History:
; 
;   2017-07-28 : MGL. New version based on chromis::extractstates.
;
;   2017-07-06 : MGL. New keyword polcal. Do not sort files.
; 
;   2018-04-16 : MGL. Make old-style CRISP files a special case,
;                working more like in the old code base.
;   
;
;-
pro crisp::extractstates, strings, states $
                          , force = force $
                          , strip_wb = strip_wb $
                          , strip_settings = strip_settings $
                          , polcal = polcal

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
 
  Nstrings = n_elements(strings)
  if( Nstrings eq 0 ) then return

  ;; If the first string ends with a number, assume all the strings
  ;; are old-style CRISP file names.
;  if (strsplit(strings[0], '.', /extract, count = Nsplit))[Nsplit-1]
;  ne 'fits' then begin

  ;; Create array of structs to hold the state information
  if keyword_set(polcal) then begin
    states = replicate( {CRISP_POLCAL_STATE}, Nstrings )
  endif else begin
    states = replicate( {CRISP_STATE}, Nstrings )
  endelse
  states.nframes = 1            ; single frame by default



  print, inam + ' : Try to extract state info from cache'
  for ifile = 0, Nstrings-1 do begin

    if keyword_set(force) then $
       cnt = 0 $
    else $
       this_cache = rdx_cacheget(strings[ifile], count = cnt)
    
    if cnt gt 0 then states[ifile] = this_cache.state

  endfor                        ; ifile

  ;; Anything not already cached?
  ncindx = where(states.filename ne strings, Nnotcached)

  if Nnotcached eq 0 then return


  
  if strmatch(strmid(strings[0],strlen(strings[0])-1,1),'[0-9]') then begin


    print, inam + ' : Get un-cached state info from file names'
    ;; Most of the info is in the file names
    if keyword_set(polcal) then begin
      red_extractstates, strings[ncindx] $
                         , basename  = basename   $
                         , cam       = detector   $
                         , dwav      = dwav       $
                         , focus     = focus      $
                         , fullstate = fullstate  $
;                         , hscan     = hscan      $
                         , lambda    = lambda     $
                         , lc        = lc         $
                         , lp        = lp         $
                         , nums      = nums       $
                         , pref      = pref       $
                         , qw        = qw         $
;                         , rscan     = rscan      $
                         , scan      = scan       $
                         , wav       = wav        
      states[ncindx].qw = qw
      states[ncindx].lp = lp
    endif else begin
      red_extractstates, strings[ncindx] $
                         , basename  = basename   $
                         , cam       = detector   $
                         , dwav      = dwav       $
                         , focus     = focus      $
                         , fullstate = fullstate  $
;                         , hscan     = hscan      $
                         , lambda    = lambda     $
                         , lc        = lc         $
                         , nums      = nums       $
                         , pref      = pref       $
;                         , rscan     = rscan      $
                         , scan      = scan       $
                         , wav       = wav        

;      print, strings[ncindx]
;      print, wav, dwav, lambda
    endelse

    states[ncindx].filename         = strings[ncindx] ; File names in input strings
    states[ncindx].detector         = detector        ; "camXIX" 
    states[ncindx].framenumber      = nums            ; Frame numbers
    states[ncindx].lc               = lc              ; Liquid crystal state
    states[ncindx].nframes          = 1               ; Single frame in old CRISP files
    states[ncindx].pf_wavelength    = lambda          ; [nm] Prefilter wavelength
    states[ncindx].prefilter        = pref            ; "6302"
    states[ncindx].scannumber       = scan            ; Scan number

    
    ;; Some info is in the ANA headers
    for ifile = 0L, Nnotcached-1 do begin
      
      if ifile mod 100 eq 0 then $
         red_progressbar, ifile, Nnotcached, 'Get remaining state info from ANA headers', /predict

      anahdr = fzhead(strings[ncindx[ifile]])    ; Read the ANA header

      ;; Exposure time [s]
      tspos = strpos(anahdr, 'Ts=')
      tepos = strpos(anahdr, 'Te=')
      if tspos ne -1 and tepos ne -1 then begin
        Ts = strmid(anahdr, tspos+3, 26)
        Te = strmid(anahdr, tepos+3, 26)
        states[ncindx[ifile]].exposure = red_time2double(strmid(Te, strpos(Te, ' '))) $
                                         - red_time2double(strmid(Ts, strpos(Ts, ' ')))
      endif

      ;; Camera "Crisp-W" etc.
      states[ncindx[ifile]].camera = red_strreplace((stregex(anahdr,'(\[)(CRISP-[WDTR]+)(\])' $
                                                             , /extr, /subexp))[2] $
                                                    , 'CRISP', 'Crisp') 
      states[ncindx[ifile]].is_wb = strmatch(states[ncindx[ifile]].camera,'*-[DW]') 

      if states[ncindx[ifile]].is_wb then begin
        states[ncindx[ifile]].tun_wavelength = states[ncindx[ifile]].pf_wavelength ; [nm] Tuning wavelength
        states[ncindx[ifile]].tuning         = states[ncindx[ifile]].prefilter+'_+0'      ; "6302_+0"
      endif else begin
        states[ncindx[ifile]].tun_wavelength = dwav[ifile]*1e-10 ; [nm] Tuning wavelength
        states[ncindx[ifile]].tuning         = wav[ifile]        ; "6301_+100"
      endelse

      ;; For CHROMIS, states.fullstate is the settings and the fpi_state
      ;; (e.g. "12.00ms_G10.00_3934_3934_+0") but for (old) CRISP we
      ;; don't have an exact expsure time setting and no info about the
      ;; camera gain. On the other hand, the LC state is important and
      ;; known. So construct it like this: "6302_6302_+0_lc0'" (or
      ;; without the LC part for WB).
      states[ncindx[ifile]].fullstate = states[ncindx[ifile]].prefilter $
                                        + '_' + states[ncindx[ifile]].tuning
      if ~states[ncindx[ifile]].is_wb then $
         states[ncindx[ifile]].fullstate += '_' + 'lc' + strtrim(long(states[ncindx[ifile]].lc), 2)
    endfor                      ; ifile

    ;; Add polcal state if neccessary
    if keyword_set(polcal) then $
       states[ncindx].fullstate = 'lp' + string(lp, format = '(i03)') + '_' $
                                  + 'qw' + string(qw, format = '(i03)') + '_' $
                                  + states[ncindx].fullstate
    
    
;   states[ncindx].skip             =  ;

    ;; states.fpi_state is the prefilter and the tuning
    ;; "6302_6301_+100", NB info also for WB
    states[ncindx].fpi_state = states[ncindx].prefilter + '_' + wav

    ;; states.settings is exposure time and camera gain for CHROMIS
    ;; (e.g., "12.00ms_G10.00") but we don't know the gain for CRISP.
    ;; So use only exposure time!
    states[ncindx].cam_settings = string(states[ncindx].exposure*1e3 $
                                         , format = '(f05.2)')+'ms' 

  endif

  ;; If we get to this point, the strings are FITS files. So the
  ;; following code should be revised once CRISP is run with the new
  ;; camera system.
  
  ;; Some info from file names
  if keyword_set(polcal) then begin
    red_extractstates,strings, lc=lc $
                      , wav = wav, pref = pref $
                      , qw=qw, lp=lp
    states.qw = qw
    states.lp = lp
  endif else begin
    red_extractstates,strings, lc=lc $
                      , wav = wav, pref = pref
  endelse
  states.lc = lc
;  states.fpi_state = pref+'_'+wav
  states.fpi_state = wav

  ;; Are the strings actually names of existing files? Then look in
  ;; the headers (for some info).
  AreFiles = min(file_test(strings))

  ;; Read headers and extract information. This should perhaps return
  ;; an array the length of the number of frames rather than the
  ;; number of files?

  for ifile = 0, Nstrings-1 do begin

    if file_test(strings[ifile]) then begin
      head = red_readhead(strings[ifile], /silent)
      states[ifile].filename = strings[ifile]
    endif else begin
      mkhdr, head, ''           ; create a dummy header
      head = red_meta2head(head, metadata = {filename:strings[ifile]})
      print,'file does not exist: ', strings[ifile]
      print,'state information will be incomplete!'
    endelse

    ;; Numerical keywords
    naxis3 = fxpar(head, 'NAXIS3', count=hasnframes)
    if hasnframes then states[ifile].nframes = naxis3
;    states[ifile].gain = fxpar(head, 'GAIN', count=hasgain)
;    ;; New keyword for gain from 2016.08.30:
;    if hasgain eq 0 then states[ifile].gain = fxpar(head, 'DETGAIN', count=hasgain) 
    states[ifile].exposure = fxpar(head, 'XPOSURE', count=hasexp)
    texposur = fxpar(head, 'TEXPOSUR', count=hastexp)

    red_fitspar_getwavelnth, head, wavelnth = wavelnth, haswav = haswav
    if haswav ne 0 then begin
;      if wavelnth ne '        ' then $
      states[ifile].pf_wavelength = float(wavelnth)
    endif

    states[ifile].scannumber = fxpar(head, red_keytab('scannumber'))
    states[ifile].framenumber = fxpar(head, red_keytab('framenumber'))

;    state = fxpar(head, 'STATE', count=hasstate)
;    if hasstate gt 0 then begin
;      state_split = strsplit( state, '_',  /extr )
;      if n_elements(state_split) gt 1 then begin
;        states[ifile].tuning = state_split[1]
;      endif
;    endif

    ;; String keywords require checking
    detector = fxpar(head, red_keytab('detector'), count=count)
    if count eq 0 then begin    ; Temporary fugly hack to work on data processed before 2016-08-23
      detector = fxpar(head, red_keytab('camera'), count=count)
    endif
    if count gt 0 then states[ifile].detector = strtrim(detector, 2)
    filter = fxpar(head, red_keytab('prefilter'), count=count)
    if count gt 0 then states[ifile].prefilter = strtrim(filter, 2)
    camera = fxpar(head, red_keytab('camera'), count=count)
    ;;camera = fxpar(head, red_keytab('old_channel'), count=count)   ; Temporary fugly hack to work on data processed before 2016-08-23

    if count gt 0 then begin
      ;; The camera is given in the header
      camera = strtrim(camera,2)
      states[ifile].camera = camera
;      states[ifile].is_wb = strmatch(states[ifile].camera,'*-[DW]') 
;      if camera eq 'Chromis-W' or camera eq 'Chromis-D' then states[ifile].is_wb = 1
    endif else begin
      ;; Try to get the camera from the detector
      self->getdetectors
      indx = where(strmatch(*self.detectors,strtrim(detector, 2)),count) 
      if count eq 1 then begin
        camera = (*self.cameras)[indx[0]]
        states[ifile].camera = camera
      endif 
    endelse
    states[ifile].is_wb = strmatch(states[ifile].camera,'*-[DW]') 

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

;    if hasgain gt 0 then begin
;      if hasexp gt 0 then states[ifile].cam_settings += '_'
;      states[ifile].cam_settings += 'G' + string(states[ifile].gain, format = '(f05.2)')
;    endif

    ;; Replace the following regexp code when this info is in the
    ;; header.

    fname = file_basename(strings[ifile], '.fits')

    ;; The focus field is an f followed by a sign and at least one
    ;; digit for the amount of focus (in ?unit?). The focus added by
    ;; the AO system (in order to compensate for prefilters with
    ;; optical power). Unit?
    focus = (stregex(fname $
                     , '(\.|^)(F[+-][0-9]+)(\.|$)' $
                     , /extr, /subexp, /fold_case))[2,*]
    ;; states.focus = focus   TODO: add field to states struct (in an SST class?)

    if ~keyword_set(polcal) then begin
      ;; The CRISP tuning information consists of a four digit
      ;; wavelength (in Å) followed by an underscore, a sign (+ or -),
      ;; and at least one digit for the finetuning (in mÅ).
      state = fxpar(head, 'STATE')
      if strtrim(state,2) eq '' || strmid(state, 0, 1) eq '_' then begin
        ;; Probably a dark frame, no tuning
      endif else begin
        tuninfo = stregex(state $
                          , '([0-9][0-9][0-9][0-9])_([+-][0-9]*)' $
                          , /extract, /subexpr) 
        
        ;; Make tuning without zero-padding in the finetuning part
        split_tuning = strsplit(tuninfo[0], '_', /extract)
        states[ifile].tuning = split_tuning[0] + '_' + strmid(split_tuning[1],0,1) + strtrim(round(abs(split_tuning[1])),2)
      endelse
    endif

    if states[ifile].tuning eq '0000_+0' then states[ifile].tuning = ''

    ;; The fullstate string
    undefine, fullstate_list
    ;; No cam settings in CRISP fullstate because exposure time varies
    ;; while observing.
;    if ~keyword_set(strip_settings) then red_append, fullstate_list, states[ifile].cam_settings
    if keyword_set(polcal) then begin
      red_append, fullstate_list, 'lp'+string(round(states[ifile].lp), format = '(i03)')
      red_append, fullstate_list, 'qw'+string(round(states[ifile].qw), format = '(i03)')
;      red_append, fullstate_list, states[ifile].lp
;      red_append, fullstate_list, states[ifile].qw
    endif
    if states[ifile].prefilter ne '' then red_append, fullstate_list, states[ifile].prefilter
    if states[ifile].tuning ne '' then begin     
      if keyword_set(strip_wb) then begin
        if states[ifile].is_wb eq 0 then $
           red_append, fullstate_list, states[ifile].tuning
      endif else begin
        red_append, fullstate_list, states[ifile].tuning
      endelse
    endif
    if states[ifile].is_wb eq 0 then red_append, fullstate_list, 'lc'+strtrim(lc[ifile], 2)
    if n_elements(fullstate_list) gt 0 then states[ifile].fullstate = strjoin(fullstate_list, '_')
    
    red_progressbar, ifile, Nstrings, 'Extract state info from file headers', /predict

    ;; Some things differ between NB and WB
    if states[ifile].is_wb then begin
      states[ifile].tun_wavelength = states[ifile].pf_wavelength
    endif else begin
      tuninfo = stregex(wav[ifile] $
                        , '([0-9][0-9][0-9][0-9])_([+-][0-9]*)' $
                        , /extract, /subexpr) 
      states[ifile].tun_wavelength = total(double(tuninfo[1:2])*[1d-10, 1d-13])
    endelse
    
  endfor                        ; ifile

;     ;; Some things differ between NB and WB
;     wbindx = where(states.is_wb, Nwb, comp=nbindx, Ncomp = Nnb)
;     if Nwb gt 0 then begin
;   
;       ;; Keywords LC and FPI_STATE should specify the NB state also for
;       ;; WB data! Fix this!!
;   
;   ;    states[wbindx].lc = ''
;   ;    states[wbindx].fpi_state = pref[wbindx]+'_'+pref[wbindx]+'+0'
;       states[wbindx].tun_wavelength = states[wbindx].pf_wavelength
;     endif
;     if Nnb gt 0 then begin
;   
;       stop
;   ;    states[nbindx].lc = lc[nbindx]
;   ;    states[nbindx].fpi_state = pref[nbindx]+'_'+wav[nbindx]
;       states[nbindx].tun_wavelength = double(pref[nbindx])*1d-10 + double(wav[nbindx])*1d-13
;     endif

  ;; Store in cache
  for ifile = 0L, Nnotcached-1 do $
     rdx_cache, strings[ncindx[ifile]], { state:states[ncindx[ifile]] }
  

  return


end


a = crispred('config.txt', /dev)


;; Test polcal
files = file_search('/data/2016/2016.09/2016.09.19/Polcal/8542/12:44:29/Crisp-T/cam*', count = Nfiles)
a -> extractstates, files, pstates, /polcal

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
