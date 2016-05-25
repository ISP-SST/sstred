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
; 
; :Keywords:
; 
;     states : out, optional, type=array(struct)
;
;        An array of structs, containing (partially filled) state information.
; 
;     cam : in, optional, type=boolean
;
;        Return the camera tag. (string)
; 
;     scannumber : in, optional, type=boolean
;
;        Return the scan-number. (int)
; 
;     prefilter : in, optional, type=boolean
;
;        Return the prefilter. (string)
; 
;     pf_wavelength : in, optional, type=boolean
;
;        Return the prefilter wavelength. (float)
; 
;     tuning : in, optional, type=boolean
;
;        Return the wavelength tuning information in the form
;        tuning_finetuning, where tuning is the approximate wavelength
;        in Å and finetuning is the (signed) fine tuning in mÅ. (string)
; 
;     tun_wavelength : in, optional, type=boolean
;
;        Return the tuning in decimal form [Å]. (double)
; 
;     lc : in, optional, type=boolean
;
;        Return the LC state. (string)
;
;     framenumber : in, optional, type=boolean
;
;        Return the frame-number. (long)
; 
;     fullstate : in, optional, type=boolean
;
;        Return the fullstate field.
; 
;     gain : in, optional, type=boolean
;
;        Return the camera gain. 
; 
;     exposure : in, optional, type=boolean
;
;        Return the exposure time. 
; 
;     focus : in, optional, type=boolean
;
;        The focus added by the AO system (in order to compensate for
;        prefilters with optical power). Unit?
; 
;     basename : in, optional, type=boolean
;
;        Set this to remove directory information from the beginning
;        of the strings.
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
;
; 
;-
pro chromis::extractstates, strings $
                        , states $
                        , cam = cam $
                        , scannumber = scannumber $
                        , prefilter = prefilter $
                        , pf_wavelength = pf_wavelength $
                        , tuning = tuning $
                        , tun_wavelength = tun_wavelength $
                        , framenumber = framenumber $
                        , fullstate = fullstate $
                        , focus = focus $
                        , basename = basename


    if keyword_set(basename) then begin
        strlist = file_basename(strings)
    endif else strlist = strings

    if fullstate then begin
       gain = 1
       exposure = 1
       tuning = 1
       prefilter = 1
    endif
    
    nt = n_elements(strings)

    if( nt eq 0 ) then return

    ;; Are the strings actually names of existing files? Then look in
    ;; the headers (for some info).
    AreFiles = min(file_test(strings))

    if AreFiles then begin
       
       ;; The scan number is the only field that is exactly five digits
       ;; long:
       if keyword_set(scannumber) then $
          scan_list = strarr(nt)

       ;; The focus field is an f folowed by a sign and at least one digit
       ;; for the amount of focus (in ?unit?):
       if keyword_set(focus) then $
          focus_list = strarr(nt) 

       ;; The frame number is the last field iff it consists entirely of
       ;; digits. The third subexpression of the regular expression matches
       ;; only the end of the string because that's where it is if it is
       ;; present. We do not know the length of the frame number field so
       ;; if the third subexpression were allowed to match a dot we would
       ;; get false matches with the scan and prefilter fields.
       if keyword_set(framenumber) then $
          num_list = strarr(nt)

       ;; The camera name consists of the string 'cam' followed by a roman
       ;; number.
       if keyword_set(cam) then $
          cam_list = strarr(nt) 

       ;; The prefilter is the only field that is exactly four digits
       if keyword_set(prefilter) or keyword_set(wavelength) then $
          prefilter_list = strarr(nt) 

       ;; The tuning information consists of a four digit wavelength (in Å)
       ;; followed by an underscore, a sign (+ or -), and at least one
       ;; digit for the finetuning (in mÅ).
       if keyword_set(tuning) or keyword_set(dwav) then $
          tuning_list = strarr(nt) 
 
       ;; The gain information is a floating point number in the
       ;; header but is really an integer.
       if keyword_set(gain) then $
          gain_list = strarr(nt)

       ;; The exposure time is a floating point number of seconds.
       if keyword_set(exposure) then $
          exposure_list = strarr(nt)

       ;; Read headers and extract information.
       ;; This should perhaps return an array the length of the number
       ;; of frames rather than the number of files?
       for ifile = 0, nt-1 do begin

          head = red_readhead(strings[ifile])
          fname = file_basename(strings[ifile], '.fits')

          if keyword_set(gain) then gain_list[ifile] = fxpar(head, 'GAIN')

          if keyword_set(exposure) then exposure_list[ifile] = fxpar(head, 'XPOSURE')

          ;; Replace the following regexp expressions when this info
          ;; is in the header.

          if keyword_set(scannumber) then $
             scan_list[ifile] = (stregex(fname $
                                         , '(\.|^)([0-9]{5})(\.|$)' $
                                         , /extr, /subexp))[2,*]

          if keyword_set(focus) then $
             focus_list[ifile] = (stregex(fname $
                                          , '(\.|^)(F[+-][0-9]+)(\.|$)' $
                                          , /extr, /subexp, /fold_case))[2,*]

     
          if keyword_set(framenumber) then $
             num_list[ifile] = (stregex(fname $
                                        , '(\.)([0-9]+)($)' $
                                        , /extr, /subexp))[2,*]

          if keyword_set(cam) then $
             cam_list[ifile] = (stregex(fname  $
                                        , '(\.|^)(cam[IVX]+)(\.|$)' $
                                        , /extr, /subexp))[2,*]

          if keyword_set(prefilter) or keyword_set(wavelength) then $
             prefilter_list[ifile] = (stregex(fname $
                                              , '(\.|^)([0-9]{4})(\.|$)'  $
                                              , /extr, /subexp))[2,*]

          if keyword_set(tuning) or keyword_set(dwav) then $
             tuning_list[ifile] = (stregex(fname $
                                           , '(\.|^)([0-9][0-9][0-9][0-9]_[+-][0-9]+)(\.|$)' $
                                           ,  /extr, /subexp))[2,*]


       endfor                   ; ifile

    endif else begin            ; AreFiles

       ;; Quantities extracted directly from the strings ----

       ;; The regular expressions below consist of three subexpressions
       ;; enclosed in parentheses. The first subexpression matches the
       ;; beginning of the line or a dot (the field separator). The second
       ;; subexpression matches the field we want to extract. The third
       ;; subexpression matches a dot or the end of the string. (With one
       ;; exception, see the frame number field below.) When calling
       ;; stregex with such a regular expression and the flag /subexpr, we
       ;; get a 2D string array, the second column of which is the field we
       ;; want to extract.

       ;; In the fields identified by a text label, like focus, lc, qw,
       ;; etc, we match with /fold_case. So if the case changes in the file
       ;; names we are OK. If another field is defined, with the same label
       ;; except for case, then this will fail. So don't do that!

       ;; The scan number is the only field that is exactly five digits
       ;; long:
       if keyword_set(scannumber) then $
          scan_list = reform((stregex(strlist,'(\.|^)([0-9]{5})(\.|$)', /extr, /subexp))[2,*])

       ;; The focus field is an f folowed by a sign and at least one digit
       ;; for the amount of focus (in ?unit?):
       if keyword_set(focus) then $
          focus_list = reform((stregex(strlist,'(\.|^)(F[+-][0-9]+)(\.|$)', /extr, /subexp, /fold_case))[2,*])

       ;; The frame number is the last field iff it consists entirely of
       ;; digits. The third subexpression of the regular expression matches
       ;; only the end of the string because that's where it is if it is
       ;; present. We do not know the length of the frame number field so
       ;; if the third subexpression were allowed to match a dot we would
       ;; get false matches with the scan and prefilter fields.
       if keyword_set(framenumber) then $
          num_list = reform((stregex(strlist,'(\.)([0-9]+)($)', /extr, /subexp))[2,*])

       ;; The camera name consists of the string 'cam' followed by a roman
       ;; number.
       if keyword_set(cam) then $
          cam_list = reform((stregex(strlist,'(\.|^)(cam[IVX]+)(\.|$)', /extr, /subexp))[2,*])

       ;; The prefilter is the only field that is exactly four digits
       if keyword_set(prefilter) or keyword_set(wavelength) then $
          prefilter_list = reform((stregex(strlist,'(\.|^)([0-9]{4})(\.|$)', /extr, /subexp))[2,*])

       ;; The tuning information consists of a four digit wavelength (in Å)
       ;; followed by an underscore, a sign (+ or -), and at least one
       ;; digit for the finetuning (in mÅ).
       if keyword_set(tuning) or keyword_set(dwav) then $
          tuning_list = reform((stregex(strlist,'(\.|^)([0-9][0-9][0-9][0-9]_[+-][0-9]+)(\.|$)', /extr, /subexp))[2,*])

    endelse                     ; AreFiles

    ;; Quantities calculated from the extracted quantities ----

    ;; The tuning as a single double precision number (in Å)
    if keyword_set(tun_wavelength) then begin
        tun_wavelength_list = dblarr(nt)
        dfac = [1d, 1d-3]
        ;; In IDL v.8 this could be done without a loop using strsplit
        for ii = 0L, nt -1 do tun_wavelength_list[ii] = total(double(strsplit(tuning_list[ii],'_', /extract))*dfac)
    endif

    if keyword_set(fullstate) then $
        fullstate_list = strjoin(transpose([[prefilter_list], [tuning_list]]), '.')

    if keyword_set(pf_wavelength) then $
        pf_wavelength_list = float(prefilter_list)*1e-10


    ;; Create array with state information
    states = replicate( {CHROMIS_STATE}, nt )
    
    states.camtag = cam_list
    states.filename = strings
    if keyword_set(gain) then states.gain = gain_list
    if keyword_set(exposure) then states.exposure = exposure_list
    if keyword_set(scannumber) then states.scannumber = scan_list
    if keyword_set(framenumber) then states.framenumber = num_list
    if keyword_set(tuning) then states.tuning = tuning_list
    if keyword_set(prefilter) then states.prefilter = prefilter_list
    if keyword_set(pf_wavelength) then states.pf_wavelength = pf_wavelength_list
    if keyword_set(dwav) then states.tun_wavelength = tun_wavelength_list
    if keyword_set(lc) then states.lc = lc_list
    if keyword_set(fullstate) then states.fullstate = fullstate_list
    ;if keyword_set(focus) then states.focus = focus_list   TODO: add field to state struct (in an SST class?)
  
end
