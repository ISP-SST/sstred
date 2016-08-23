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
;-
pro chromis::extractstates, strings, states $
                            , strip_wb = strip_wb $
                            , strip_settings = strip_settings
  
  Nstrings = n_elements(strings)
  if( Nstrings eq 0 ) then return

  ;; Create array of structs to holed the state information
  states = replicate( {CHROMIS_STATE}, Nstrings )

  ;; Are the strings actually names of existing files? Then look in
  ;; the headers (for some info).
  AreFiles = min(file_test(strings))

  ;; Read headers and extract information. This should perhaps return
  ;; an array the length of the number of frames rather than the
  ;; number of files?
  progress_message = 'Extract state info from file headers'
  red_progressbar, 0, Nstrings, message = progress_message
  for ifile = 0, Nstrings-1 do begin

     if file_test(strings[ifile]) then begin
        head = red_readhead(strings[ifile], /silent)
        states.filename = strings
     endif else begin
        mkhdr, head, ''         ; create a dummy header
        head = red_filterchromisheaders( head, meta={filename:strings[ifile]}, /silent)
        print,'file does not exist: ', strings[ifile]
        print,'state information will be incomplete!'
     endelse

     ;; Numerical keywords
     states[ifile].gain = fxpar(head, 'GAIN', count=hasgain)
     states[ifile].exposure = fxpar(head, 'XPOSURE', count=hasexp)
     states[ifile].pf_wavelength = fxpar(head, 'WAVELNTH')
     states[ifile].scannumber = fxpar(head, red_keytab('scannumber'))
     states[ifile].framenumber = fxpar(head, red_keytab('framenumber'))

     ;; String keywords require checking
     detector = fxpar(head, red_keytab('detector'), count=count)
     if count eq 0 then begin   ; Temporary fugly hack to work on data processed before 2016-08-23
        detector = fxpar(head, red_keytab('camera'), count=count)
     endif
     if count gt 0 then states[ifile].detector = strtrim(detector, 2)
     filter = fxpar(head, red_keytab('prefilter'), count=count)
     if count gt 0 then states[ifile].prefilter = strtrim(filter, 2)
     ;camera = fxpar(head, red_keytab('camera'), count=count)
     camera = fxpar(head, red_keytab('old_channel'), count=count)   ; Temporary fugly hack to work on data processed before 2016-08-23
     if count gt 0 then begin
         camera = strtrim(camera,2)
         states[ifile].camera = camera
         if camera eq 'Chromis-W' or camera eq 'Chromis-D' then states[ifile].is_wb = 1
     endif

     if hasexp gt 0 then begin
        states[ifile].cam_settings = strtrim(string(states[ifile].exposure*1000 $
                                                    , format = '(f9.2)'), 2) + 'ms'
     endif
     if hasgain gt 0 then begin
         if hasexp gt 0 then states[ifile].cam_settings += '_'
         states[ifile].cam_settings += 'G' + string(states[ifile].gain, format = '(f05.2)')
     endif
     
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


     ;; The tuning information consists of a four digit wavelength (in
     ;; Å) followed by an underscore, a sign (+ or -), and at least
     ;; one digit for the finetuning (in mÅ). Return the wavelength
     ;; tuning information in the form tuning_finetuning, where tuning
     ;; is the approximate wavelength in Å and finetuning is the
     ;; (signed) fine tuning in mÅ. (string)
         
     states[ifile].tuning = (stregex(fname $
                                     , '(\.|^)([0-9][0-9][0-9][0-9]_[+-][0-9]+)(\.|$)' $
                                     ,  /extr, /subexp))[2,*]
     
     ;; Convert to a wavelength
     ;; Return the tuning in decimal form [Å]
     states[ifile].tun_wavelength = total(double(strsplit(states[ifile].tuning,'_', /extract))*[1d, 1d-3])

     ;; The fullstate string
     undefine, fullstate_list
     if ~keyword_set(strip_settings) then red_append, fullstate_list, states[ifile].cam_settings
     if states[ifile].prefilter ne '' then red_append, fullstate_list, states[ifile].prefilter
     if keyword_set(strip_wb) then begin
        if states[ifile].is_wb eq 0  and states[ifile].tuning ne '' then $
           red_append, fullstate_list, states[ifile].tuning
     endif else begin
        if states[ifile].tuning ne '' then  $
           red_append, fullstate_list, states[ifile].tuning
     endelse
     states[ifile].fullstate = strjoin(fullstate_list, '_')

     red_progressbar, ifile, Nstrings, message = progress_message
     
  endfor                        ; ifile
  red_progressbar, /finished, message = progress_message

end
