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
;     states : out, optional, type=struct
;
;        A struct containing some of the state information as string
;        arrays. 
; 
;     fullstate : out, optional, type=strarr
;
;        The states in the form pref.wav.lc (see those keywords).
; 
;     nums : out, optional, type=strarr
;
;        The frame numbers.
; 
;     wav : out, optional, type=strarr
;
;        The wavelength tuning information in the form
;        tuning_finetuning, where tuning is the approximate wavelength
;        in Å and finetuning is the (signed) fine tuning in mÅ.
; 
;     lc : out, optional, type=strarr
;
;        The LC state.
; 
;     pref : out, optional, type=strarr
;
;        The prefilter wavelength in Å.
; 
;     scan : out, optional, type=strarr
;
;        The scan number.
; 
;     dwav : out, optional, type=strarr
;
;        The tuning in decimal form [Å].
; 
;     lambda : out, optional, type=fltarr
;
;        The prefilter wavelength in decimal form [m].
; 
;     cam : out, optional, type=strarr
;
;        The camera name.
; 
;     rscan : out, optional, type=strarr
;
;        ?? Can be removed when removed from prepmomfbd and prepmomfbd2.
; 
;     hscan : out, optional, type=strarr
;
;        ?? Can be removed when removed from prepmomfbd and prepmomfbd2.
;  
;     focus : out, optional, type=strarr
;
;        The focus added by the AO system (in order to compensate for
;        prefilters with optical power). Unit?
; 
;     basename : in, optional, type=boolean
;
;        Set this to remove directory information from the beginning
;        of the strings.
;
;     blue : in, optional, type=boolean
; 
;        Set this if it's data from a blue tower or tilt filter
;        camera. Will change what is returned in the fullstate and wav
;        keywords. Will not return anything in many keywords that are
;        only relevant for CRISP data.  
; 
; 
; :History:
; 
;   2014-01-22 : First version.
; 
;   2014-04-?? : MGL. Added keyword lambda.
;
;   2015-09-03 : MGL. Added the "blue" keyword, will now return
;                something meaningful also for blue tilt filter data.
;                Bugfix in qw regular expression.
;
;   2018-06-11 : MGL. Add conversion from CHROMIS wheel status to
;                prefilter. New keywords is_wb and is_pd.
;
; 
;-
pro red_extractstates, strings $
                       , is_wb = is_wb $
                       , is_pd = is_pd $
                       , states = states $
                       , pstates = pstates $
                       , pstates_out = pstates_out $
                       , fullstate = fullstate $
                       , nums = nums $
                       , wav = wav $
                       , qw = qw $
                       , lc = lc $
                       , lp = lp $
                       , cam = cam $
                       , pref = pref $
                       , scan = scan $
                       , dwav = dwav $
                       , lambda = lambda $
                       , rscan = rscan $
                       , hscan = hscan  $
                       , focus = focus $
                       , blue = blue $
                       , basename = basename

  nt = n_elements(strings)

  if keyword_set(basename) then begin
    strlist = file_basename(strings)
  endif else strlist = strings

  if nt eq 0 then return

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
  if arg_present(scan) or arg_present(rscan)  or $
     arg_present(hscan) or arg_present(states) then $
        scan = reform((stregex(strlist,'(_|\.|^)([0-9]{5})(_|\.|$)', /extr, /subexp))[2,*])

  ;; The focus field is an f folowed by a sign and at least one digit
  ;; for the amount of focus (in ?unit?):
  if arg_present(focus) then $
     focus = reform((stregex(strlist,'(_|\.|^)(F[+-][0-9]+)(_|\.|$)', /extr, /subexp, /fold_case))[2,*])

  ;; The frame number is the last field iff it consists entirely of
  ;; digits. The third subexpression of the regular expression matches
  ;; only the end of the string because that's where it is if it is
  ;; present. We do not know the length of the frame number field so
  ;; if the third subexpression were allowed to match a dot we would
  ;; get false matches with the scan and pref fields.
  if arg_present(nums) or arg_present(states) or arg_present(pstates) then $
     nums = reform((stregex(strlist,'(_|\.)([0-9]+)($)', /extr, /subexp))[2,*])

  ;; The camera name consists of the string 'cam' followed by a roman
  ;; number.
  if arg_present(cam) or arg_present(pstates_out) then $
     cam = reform((stregex(file_basename(strlist),'(_|\.|^)(cam[IVX]+)(_|\.|$)', /extr, /subexp))[2,*])
  
  ;; The prefilter is the only field that is exactly four digits
  if arg_present(pref) or arg_present(fullstate) or arg_present(lambda) or $
     arg_present(states) or arg_present(pstates) or arg_present(pstates_out) or $
     arg_present(is_wb) or arg_present(is_pd) $
  then begin

    pref = reform((stregex(strlist,'(_|\.|^)([0-9]{4})(_|\.|$)', /extr, /subexp))[2,*])
    is_wb = bytarr(n_elements(pref))
    is_pd = bytarr(n_elements(pref))

    ;; Chromis files have "wheel0000i" instead of prefilter
    indx = where(pref eq '', Nwhere)
    for iwhere = 0, Nwhere-1 do begin
      camera = stregex(strings[indx[iwhere]], 'Chromis-[WDN]', /extract)
      wheel = stregex(strlist[indx[iwhere]], 'wheel0000[0-9]', /extract)
      if camera ne '' and wheel ne '' then begin
        if camera eq 'Chromis-N' then begin
          ;; Chromis-N
          case wheel of
            'wheel00001' : pref[iwhere] = '3925' ; Ca II K blue wing
            'wheel00002' : pref[iwhere] = '3934' ; Ca II K core
            'wheel00003' : pref[iwhere] = '3969' ; Ca II H core
            'wheel00004' : pref[iwhere] = '3978' ; Ca II H red wing
            'wheel00005' : pref[iwhere] = '3999' ; Ca II H continuum
            'wheel00006' : pref[iwhere] = '4862' ; H-beta core
            else :
          endcase
        endif else begin
          ;; Chromis-W and Chromis-D
          case wheel of
            'wheel00006' : pref[iwhere] = '4846' ; H-beta continuum
            else         : pref[iwhere] = '3950' ; Ca II HK wideband
          endcase
        endelse
      endif
      is_wb[iwhere] = ~strmatch(camera, '*-N')
      is_pd[iwhere] =  strmatch(camera, '*-D')
    endfor                      ; iwhere
  endif
    
  ;; The tuning information consists of a four digit wavelength (in Å)
  ;; followed by an underscore, a sign (+ or -), and at least one
  ;; digit for the finetuning (in mÅ).
  if arg_present(wav) or arg_present(dwav) or arg_present(fullstate) or $
     arg_present(states) or arg_present(pstates) then $
        wav = reform((stregex(strlist,'(_|\.|^)([0-9][0-9][0-9][0-9]_[+-][0-9]+)(_|\.|$)', /extr, /subexp))[2,*])

  ;; The LC state is the string 'LC' followed by a single digit
  if arg_present(lc) or arg_present(fullstate) or arg_present(states) or $
     arg_present(pstates) or arg_present(pstates_out) then $
        lc = long(red_strreplace(reform( $
             (stregex(strlist,'(_|\.|^)(LC[0-9])(_|\.|$)', /extr, /subexp, /fold_case))[2,*]) $
                                 , 'lc',''))
  
  ;; For polcal, the linear polarizer state
  if arg_present(lp) or arg_present(pstates) or arg_present(pstates_out) then begin
    strlst = red_strreplace(file_basename(strlist),'LP','lp')
    lp = float(red_strreplace(reform( $
         (stregex(strlst,'(_|\.|^)(LP[0-3][0-9]{2})(_|\.|$)', /extr, /subexp, /fold_case))[2,*]) $
                              , 'lp', ''))
  endif
  
  ;; For polcal, the quarter wave plate state
  if arg_present(qw) or arg_present(pstates) or arg_present(pstates_out) then $
     qw = float(red_strreplace(reform( $
          (stregex(strlist,'(_|\.|^)(QW[0-3][0-9]{2})(_|\.|$)', /extr, /subexp, /fold_case))[2,*]) $
                               , 'qw', ''))


  ;; Quantities calculated from the extracted quantities ----

  ;; The tuning as a single double precision number (in Å)
  if arg_present(dwav) or arg_present(states) then begin
    dwav = dblarr(nt)
    dfac = [1d, 1d-3]
    ;; In IDL v.8 this could be done without a loop using strsplit
    for ii = 0L, nt -1 do dwav[ii] = total(double(strsplit(wav[ii],'_', /extract))*dfac)
  endif

  if arg_present(hscan) or arg_present(rscan) or arg_present(states) then begin
    rscan = strarr(nt)
    hscan = strarr(nt)
    for ii = 0L, nt -1 do begin
      rscan[ii] = red_decode_scan(scan[ii], hscan=hs)
      hscan[ii] = hs   
    endfor
  endif
  
  if arg_present(fullstate) or arg_present(states) then $
     fullstate = strjoin(transpose([[pref], [wav], [lc]]), '.')

  if arg_present(lambda) then lambda = float(pref)*1e-10

  ;; Create structures with collections of states

  if arg_present(states) then states = { files:strlist $
                                         , nums:nums $
                                         , wav:wav $
                                         , lc:lc $
                                         , pref:pref $
                                         , state:fullstate $
                                         , star:bytarr(nt) $
                                         , scan:scan $
                                         , rscan:rscan $
                                         , hscan:hscan $
                                         , dwav:dwav $
                                       }
  
  ;; For polarimetry
  if arg_present(pstates) then pstates = { state:strjoin(transpose([[lp], [qw], [pref], [wav], [lc]]), '.') $ 
                                           , lp:lp $
                                           , qw:qw $
                                           , pref:pref $
                                           , wav:wav $
                                           , lc:lc $
                                           , nums:nums $
                                           , files:strlist $
                                           , star:bytarr(nt) $ ; Empty field, see red_flagtuning.
                                         }
  
  ;; For polarimetry out (?)
  if arg_present(pstates_out) then pstates_out = { state:strjoin(transpose([[lps], [qws], [pref], [lcs]]), '.') $ 
                                                   , lp:lp $
                                                   , qw:qw $
                                                   , pref:pref $
                                                   , wav:wav $
                                                   , lc:lc $
                                                   , lcs:strmid(lc, 2, 1) $
                                                   , lps:lp $
                                                   , qws:qw $
                                                   , cam:cam $
                                                 }
  
end
