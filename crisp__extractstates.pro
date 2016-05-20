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
; :author:
; 
;     Mats Löfdahl, ISP
; 
; 
; :returns:
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
;     polcal : in, optional, type=boolean
;
;        Return the polcal state information (linear polarizer
;        and quarter-wave plate angles). (string)
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
; :history:
; 
;   2014-01-22 : First version.
; 
;   2014-04-?? : MGL. Added keyword wavelength.
;
;   2015-09-03 : MGL. Added the "blue" keyword, will now return
;                something meaningful also for blue tilt filter data.
;                Bugfix in qw regular expression.
;
;   2016-05-19 : THI. Partial copy to the crisp class. Modify the state structures
;                and keywords for clarity.
;
; 
;-
pro crisp::extractstates, strings $
                        , states $
                        , cam = cam $
                        , scannumber = scannumber $
                        , prefilter = prefilter $
                        , pf_wavelength = pf_wavelength $
                        , tuning = tuning $
                        , tun_wavelength = tun_wavelength $
                        , lc = lc $
                        , framenumber = framenumber $
                        , fullstate = fullstate $
                        , polcal = polcal $
                        , focus = focus $
                        , basename = basename


    if keyword_set(basename) then begin
        strlist = file_basename(strings)
    endif else strlist = strings

    nt = n_elements(strings)

    if( nt eq 0 ) then return

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
    if keyword_set(prefilter) or keyword_set(fullstate) or keyword_set(wavelength) then $
        prefilter_list = reform((stregex(strlist,'(\.|^)([0-9]{4})(\.|$)', /extr, /subexp))[2,*])

    ;; The tuning information consists of a four digit wavelength (in Å)
    ;; followed by an underscore, a sign (+ or -), and at least one
    ;; digit for the finetuning (in mÅ).
    if keyword_set(tuning) or keyword_set(dwav) or keyword_set(fullstate) then $
        tuning_list = reform((stregex(strlist,'(\.|^)([0-9][0-9][0-9][0-9]_[+-][0-9]+)(\.|$)', /extr, /subexp))[2,*])

    ;; The LC state is the string 'LC' followed by a single digit
    if keyword_set(lc) or keyword_set(fullstate) then $
        lc_list = reform((stregex(strlist,'(\.|^)(LC[0-9])(\.|$)', /extr, /subexp, /fold_case))[2,*])

    ;; For polcal
    if keyword_set(polcal) then begin
        ;; The linear polarizer state
        lp_list = reform((stregex(strlist,'(\.|^)(LP[0-3][0-9]{2})(\.|$)', /extr, /subexp, /fold_case))[2,*])
        ;; The quarter wave plate state
        qw_list = reform((stregex(strlist,'(\.|^)(QW[0-3][0-9]{2})(\.|$)', /extr, /subexp, /fold_case))[2,*])
    endif

    ;; Quantities calculated from the extracted quantities ----

    ;; The tuning as a single double precision number (in Å)
    if keyword_set(tun_wavelength) then begin
        tun_wavelength_list = dblarr(nt)
        dfac = [1d, 1d-3]
        ;; In IDL v.8 this could be done without a loop using strsplit
        for ii = 0L, nt -1 do tun_wavelength_list[ii] = total(double(strsplit(tuning_list[ii],'_', /extract))*dfac)
    endif

    if keyword_set(fullstate) then $
        fullstate_list = strjoin(transpose([[prefilter_list], [tuning_list], [lc_list]]), '.')

    if keyword_set(pf_wavelength) then $
        pf_wavelength_list = float(prefilter_list)*1e-10


    ;; Create array with state information

    if keyword_set(polcal) then begin
        states = replicate( {CRISP_POLCAL_STATE}, nt )
        states.lp = lp_list
        states.qw = qw_list
    endif else begin
        states = replicate( {CRISP_STATE}, nt )
    endelse
    
    states.camtag = cam_list
    states.filename = strings
    if keyword_set(scannumber) then states.scannumber = scan_list
    if keyword_set(framenumber) then states.framenumber = num_list
    if keyword_set(tuning) then states.tuning = tuning_list
    if keyword_set(prefilter) then states.prefilter = prefilter_list
    if keyword_set(pf_wavelength) then states.pf_wavelength = pf_wavelength_list
    if keyword_set(dwav) then states.tun_wavelength = tun_wavelength_list
    if keyword_set(lc) then states.lc = lc_list
    if keyword_set(fullstate) then states.fullstate = fullstate_list
    ;if keyword_set(focus) then states.focus = focus_list   TODO: add field to state struct (in an SST class?)



  ;; For polarimetry
;   if arg_present(pstates) then pstates = { state:strjoin(transpose([[lp], [qw], [pref], [wav], [lc]]), '.') $ 
;                                            , lp:lp $
;                                            , qw:qw $
;                                            , pref:pref $
;                                            , wav:wav $
;                                            , lc:lc $
;                                            , nums:nums $
;                                            , files:strlist $
;                                            , star:bytarr(nt) $ ; Empty field, see red_flagtuning.
;                                          }
;   
;   ;; For polarimetry out (?)
;   if arg_present(pstates_out) then pstates_out = { state:strjoin(transpose([[lps], [qws], [pref], [lcs]]), '.') $ 
;                                                    , lp:float(strmid(lp, 2, 3)) $
;                                                    , qw:float(strmid(qw, 2, 3)) $
;                                                    , pref:pref $
;                                                    , wav:wav $
;                                                    , lc:lc $
;                                                    , lcs:strmid(lc, 2, 1) $
;                                                    , lps:lp $
;                                                    , qws:qw $
;                                                    , cam:cam $
;                                                  }
  
end
