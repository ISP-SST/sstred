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
        head = red_filterchromisheaders( head, meta={filename:strings[ifile]} )
        print,'file does not exist: ', strings[ifile]
        print,'state information will be incomplete!'
     endelse

     ;; Numerical keywords
     states[ifile].gain = fxpar(head, 'GAIN')
     states[ifile].exposure = fxpar(head, 'XPOSURE')
     states[ifile].pf_wavelength = fxpar(head, 'WAVELNTH', count=count)

     ;; String keywords require checking
     camtag = fxpar(head, red_keytab('cam'), count=count)
     if count gt 0 then states[ifile].camtag = strtrim(fxpar(head, red_keytab('cam')), 2)
     filter = fxpar(head, red_keytab('prefilter'), count=count)
     if count gt 0 then states[ifile].prefilter = filter

     ;; These keywords are temporary, change when we change in
     ;; red_filterchromisheaders.
     states[ifile].scannumber = fxpar(head, red_keytab('scannumber'))
     states[ifile].framenumber = fxpar(head, red_keytab('framenumber'))

     
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
 ;    if states[ifile].tuning ne '' then begin
 ;       states[ifile].fullstate = states[ifile].prefilter + '.' + states[ifile].tuning
 ;    endif else begin
 ;       states[ifile].fullstate = states[ifile].prefilter
 ;    endelse
 ;
 ;
     if ~keyword_set(strip_settings) then begin
        states[ifile].fullstate = string(states[ifile].exposure*1000, format = '(f4.2)') + 'ms' $
                                + '_' + 'G' + string(states[ifile].gain, format = '(f05.2)')
     endif
     if states[ifile].prefilter ne '' then states[ifile].fullstate += '_' + states[ifile].prefilter
     if keyword_set(strip_wb) then begin
        if states[ifile].camtag eq self->getcamtag('Chromis-N') $
           and states[ifile].tuning ne '' then states[ifile].fullstate += '_' + states[ifile].tuning
     endif else begin
        if states[ifile].tuning ne '' then states[ifile].fullstate += '_' + states[ifile].tuning
     endelse

     red_progressbar, ifile, Nstrings, message = progress_message
     
  endfor                        ; ifile
  red_progressbar, /finished, message = progress_message

end
