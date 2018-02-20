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
; 
;-
pro crisp::extractstates, strings, states $
                          , strip_wb = strip_wb $
                          , strip_settings = strip_settings $
                          , polcal = polcal

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
 
  Nstrings = n_elements(strings)
  if( Nstrings eq 0 ) then return

  ;; Create array of structs to hold the state information
  if keyword_set(polcal) then begin
    states = replicate( {CRISP_POLCAL_STATE}, Nstrings )
  endif else begin
    states = replicate( {CRISP_STATE}, Nstrings )
  endelse
  states.nframes = 1            ; single frame by default

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
      states.filename = strings
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

    ;; The CRISP tuning information consists of a four digit
    ;; wavelength (in Å) followed by an underscore, a sign (+ or -),
    ;; and at least one digit for the finetuning (in mÅ).
    tuninfo = stregex(fxpar(head, 'STATE') $
                      , '([0-9][0-9][0-9][0-9])_([+-][0-9]*)' $
                      , /extract, /subexpr) 
    states[ifile].tuning = tuninfo[0]
    if states[ifile].tuning eq '0000_+0' then states[ifile].tuning = ''

    ;; The fullstate string
    undefine, fullstate_list
    ;; No cam settings in CRISP fullstate because exposure time varies
    ;; while observing.
;    if ~keyword_set(strip_settings) then red_append, fullstate_list, states[ifile].cam_settings
    if keyword_set(polcal) then begin
      red_append, fullstate_list, 'LP'+strtrim(long(states[ifile].lp), 2)
      red_append, fullstate_list, 'qw'+strtrim(long(states[ifile].qw), 2)
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

  endfor                        ; ifile

  ;; Some things differ between NB and WB
  wbindx = where(states.is_wb, Nwb, comp=nbindx, Ncomp = Nnb)
  if Nwb gt 0 then begin

    ;; Keywords LC and FPI_STATE should specify the NB state also for
    ;; WB data! Fix this!!

;    states[wbindx].lc = ''
;    states[wbindx].fpi_state = pref[wbindx]+'_'+pref[wbindx]+'+0'
    states[wbindx].tun_wavelength = states[wbindx].pf_wavelength
  endif
  if Nnb gt 0 then begin
;    states[nbindx].lc = lc[nbindx]
;    states[nbindx].fpi_state = pref[nbindx]+'_'+wav[nbindx]
    states[nbindx].tun_wavelength = double(pref[nbindx])*1d-10 + double(wav[nbindx])*1d-13
  endif


end


a = crispred('config.txt', /dev)


;; Test polcal
files = file_search('/storage/sand02/Incoming/2016.09.19/Polcal/8542/12:44:29/Crisp-T/cam*', count = Nfiles)
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
