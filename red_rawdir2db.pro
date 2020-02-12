; docformat = 'rst'

;+
; Write metadata for raw datasets into the database.
;
; A dataset is defined by its timestamp directory.
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
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
;
;   date : in, optional, type=str
;
;     The date of observations.
;
;   dir : in, optional, type=str, default=all
;   
;     A top directory to start search for data files.
;
;   instrument : in, optional, type=str
;
;     An instrument to include.
;
;   all : in, optional, type=boolean
;
;     Set to include in search all subdirectories.
; 
; :History:
; 
;   2019-01-31 : MGL. First version.
;
;   2019-05-28 : OA. Second version
;
;   2019-06-06 : OA. Added check for incomplete scans (to be rejected)
;
;   2020-02-04 : OA. Added calls to different procedures for different
;                instruments.
;-
pro red_rawdir2db, all = all $
                   , date = date $
                   , dir = dir $
                   , instrument = instrument $
                   , debug = debug

  inam = red_subprogram(/low, calling = inam1)
  if n_elements(date) eq 0 and n_elements(dir) eq 0 then begin
    print, inam+' : Please provide at least one of the dir or date keywords.'
    retall
  endif
  
  if n_elements(date) eq 1 then begin
    ;; If the date keyword is given, then make a top search directory
    ;; based on it.
    red_currentsite,  site = site, search_dirs = search_dirs, date = date
    ;; Then append the dir keyword if given.
    if n_elements(dir) gt 0 then search_dirs += dir
  endif else begin
    ;; If no date is given, then check whether the dir keyword is an
    ;; absolute path or a relative path.
    if strmid(dir, 0, 1) eq '/' then begin
      ;; If an absolute path, then use it.
      search_dirs = dir
    endif else begin
      ;; If a relative path, then append it to the site's
      ;; date-less search dir.
      red_currentsite, site = site, search_dirs = search_dirs
      search_dirs += dir
    endelse
  endelse 
  
  search_dirs = red_strreplace(search_dirs + '/', '//', '/')
  
  found_dirs = file_search(search_dirs+'*', count = Nfound)

  if Nfound eq 0 then begin
    print, inam + ' : No matches when searching for '+search_dirs+'*'
    retall
  endif

  if n_elements(instrument) gt 0 then begin
    is_chromis = strmatch(found_dirs, '*CHROMIS*')
    indx_chromis = where(is_chromis, Nchromis $
                         , complement = indx_crisp, ncomplement = Ncrisp)
    case instrument of
      'CRISP' : begin
        if Ncrisp ne 0 then found_dirs = found_dirs[indx_crisp] else begin
          print, inam + ' : No data for instrument='+instrument+' in '+search_dirs+'*'
          retall
        endelse
      end
      'CHROMIS' : begin
        if Nchromis ne 0 then found_dirs = found_dirs[indx_chromis] else begin
          print, inam + ' : No data for instrument='+instrument+' in '+search_dirs+'*'
          retall
        endelse
      end
      else : begin
        print, inam + ' : Unknown instrument: '+instrument
        retall
      end
    endcase
  endif
  
  if ~keyword_set(all) then begin
    tmp = red_select_subset(found_dirs $
                            , default = '*' $
                            , qstring = 'Select directories' $
                            , indx = sindx $
                            , count = Nselected)
    if Nselected eq 0 then begin
      print, inam + ' : No directories selected.'
      retall
    endif
    dirs = found_dirs[sindx]
  endif else dirs = found_dirs
  Ndirs = n_elements(dirs)
  
  ;; Now go as deep in the selected directories as needed to find the
  ;; camera directories.
  for idir = 0, Ndirs-1 do begin
    tmpdirs = dirs[idir]
    while ~strmatch((strsplit(tmpdirs[0], '/', /extract))[-1], '*-[NWDTR]') do $
       tmpdirs = file_search(tmpdirs+'/*', count = Nfound)

    if Nfound gt 0 then red_append, camdirs, tmpdirs
    
  endfor                        ; ifound

  Ncamdirs = n_elements(camdirs)
  if Ncamdirs eq 0 then begin
    print, inam + ' : No camera directories found in '
    print, dirs
    retall
  endif
 
  for idir = 0, Ncamdirs-1 do begin

    splitdir = (strsplit(camdirs[idir],'/',/extract))
    camera = splitdir[-1]
    date = stregex(camdirs[idir], '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]', /extract)
    timestamp = stregex(camdirs[idir], '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]', /extract)
    isodate = red_strreplace(date, '.', '-', n = 2)
    date_obs = isodate+' '+timestamp ; 'T'
    
;    is_wb = strmatch(camera, '*-[WD]')
    
    ;; Get the DATATYPE from the path.
    a = stregex(camdirs[idir], '/CHROMIS-([a-zA-Z]*)/', /extract, /sub)
 ;   if n_elements(a) eq 2 then begin
    if a[0] ne '' then begin
      
      instrume = 'CHROMIS'

      datatype = a[1]
      case datatype of
        'darks'    : datatype = 'darks'
        'flats'    : datatype = 'flats'
        'pinholes' : datatype = 'pinholes'
        'data'     : datatype = 'science' 
        else : datatype = ' (unknown)'
      endcase

      ;; Array with NB prefilters to be indexed with the "wheel" number
      ;; minus 1. Same index to use for lambda_ref, du_ref, convfac.
      nbprefs = ['3925', '3934', '3969', '3978', '3999', '4862']
      Nnbprefs = n_elements(nbprefs)
      lambda_ref = dblarr(Nnbprefs)
      du_ref     = dblarr(Nnbprefs)
      convfac    = dblarr(Nnbprefs)
      for wheel = 1, Nnbprefs do begin
        zfile = './info/hrz_zeropoint_' + nbprefs[wheel-1] + '.fz'
                                ;refinfo = f0(zfile)
        refinfo = dblarr(3)
        hh = bytarr(512)
        openr,11,zfile
        readu,11,hh
        readu,11,refinfo
        close,11
        lambda_ref[wheel-1] = refinfo[0]
        du_ref[wheel-1]     = refinfo[1]
        convfac[wheel-1]    = refinfo[2]
      endfor                    ; wheel
      

      ;; Get the hrz conversion info
      ;a = stregex(camdirs[idir], '/([12][0-9][0-9][0-9][-.][01][0-9][-.][0-3][0-9])/', /extract, /sub)
;      red_download_linedefs, isodate, strjoin(splitdir[0:-2],'/'), './'
      ;chromis_hrz_zeropoint, './'
      
    endif else begin
      ;; CRISP directory structure is not strict
      instrume = 'CRISP'
      datatype = '(unknown)'
      if stregex(camdirs[idir],'Darks',/boolean) then datatype = 'darks'
      if stregex(camdirs[idir],'Flats',/boolean) then datatype = 'flats'
      if stregex(camdirs[idir],'Pinholes',/boolean) then datatype = 'pinholes'
      if stregex(camdirs[idir],'Polcal',/boolean) then datatype = 'polcal'
      if stregex(camdirs[idir],'Science',/boolean) then datatype = 'science'
      
    endelse
    
    ;; Now search for the actual files
    files = file_search(camdirs[idir]+'/*', count = Nfiles)
    
    if Nfiles eq 0 then continue
    
    ;; Make struct for info from header, that should go into the
    ;; database. Some info from the path is set right away.

    strct  = {BITPIX           : 0           $  
              , CADAVG         : 0.          $
              , CAMERA         : camera      $ 
              , DATATYPE       : datatype    $
              , DATE           : ''          $
              , DATE_BEGS      : ''          $
              , DATE_OBS       : date_obs    $
              , DETECTOR       : ''          $
              , DETFIRM        : '0.0'       $
              , DETGAIN        : 0.          $
 ;             , DETMODEL       : ''          $
              , DETOFFS        : 0L          $
 ;             , FILENAME       : ''          $
              , FILTER1        : '0'          $
              , FIRST_FRAME    : 0L           $
              , INSTRUME       : instrume    $
 ;             , IS_WB          : is_wb       $
 ;             , NAXIS          : 0L          $
              , NAXIS1         : 0L          $          
              , NAXIS2         : 0L          $
              , NAXIS3         : 0L          $
 ;             , OBJECT         : ''          $
 ;             , OBSERVER       : ''          $
 ;             , OBS_HDU        : 0B          $
 ;             , ORIGIN         : ''          $
              , SCANNUM        : 0L          $
              , SOLARNET       : 0.          $
              , STATE          : ''          $
              , WHEEL          : 0           $
              , HRZ            : 0L           $
 ;             , TELESCOP       : ''          $
              , WAVEBAND       : ''          $
              , WAVELNTH       : 0.          $
              , WAVEMAX        : 0.          $
              , WAVEMIN        : 0.          $
              , WAVEUNIT       : 0           $
              , XPOSURE        : 0.          $
              , DETTEMP        : 0.          $
              , DIR_TEMPLATE   : ''          $
              , FNM_TEMPLATE   : ''          $
              , LC_STATE       : '-1'          $ ; fake lc_state (crisp darks)
              , QW_STATE       : ''          $
              , LP_STATE       : ''          $
              , FOCUS          : 888         $
              , FRAME_MAX      : ''          $
              , FRAME_MIN      : ''          $
              , FRAME_MEDIAN   : ''          $
              , FRAME_STDDEV   : ''          $
              , NB_PREFS       : ''          $
              , LAMBDA_REF     : ''          $
              , DU_REF         : ''          $
              , CONVFAC        : ''          $
             }

    ;; Make array of such structs and set the values known without
    ;; reading the headers
    dbinfo = replicate(strct, Nfiles)
    ;dbinfo.FILENAME = file_basename(files)

    ;; Fill in other values in the array of structs from the header:
    for ifile = 0, Nfiles-1 do begin

      red_progressbar, ifile, Nfiles, /predict, 'Get DB values for '+camdirs[idir]
      
      h = red_readhead(files[ifile], date_beg=date_beg, framenumbers=framenumbers)

      if ifile eq 0 then begin
        ;; Set the keywords that absolutely do not change from file to
        ;; file within the directory.
        value = fxpar(h, 'DETECTOR', count=cnt) & if cnt eq 1 then dbinfo.DETECTOR = value 
        value = fxpar(h, 'DETFIRM' , count=cnt) & if cnt eq 1 then dbinfo.DETFIRM  = value 
;        value = fxpar(h, 'DETMODEL', count=cnt) & if cnt eq 1 then dbinfo.DETMODEL = value 
;        value = fxpar(h, 'OBJECT'  , count=cnt) & if cnt eq 1 then dbinfo.OBJECT   = value 
;        value = fxpar(h, 'OBSERVER', count=cnt) & if cnt eq 1 then dbinfo.OBSERVER = value 
;        value = fxpar(h, 'OBS_HDU' , count=cnt) & if cnt eq 1 then dbinfo.OBS_HDU  = value 
;        value = fxpar(h, 'ORIGIN'  , count=cnt) & if cnt eq 1 then dbinfo.ORIGIN   = value 
;        value = fxpar(h, 'TELESCOP', count=cnt) & if cnt eq 1 then dbinfo.TELESCOP = value
        if instrume eq 'CRISP' then begin
          crisp_fnm_gen, files[ifile], fnm_gen = fnm_gen, dir_gen = dir_gen
          dbinfo.DIR_TEMPLATE = dir_gen
          dbinfo.FNM_TEMPLATE = fnm_gen
        endif else begin ; CHROMIS
          chromis_fnm_gen, files[ifile], fnm_gen = fnm_gen, dir_gen = dir_gen
          dbinfo.DIR_TEMPLATE = dir_gen
          dbinfo.FNM_TEMPLATE = fnm_gen
          dbinfo.NB_PREFS = json_serialize(nbprefs)
          dbinfo.LAMBDA_REF = json_serialize(lambda_ref)
          dbinfo.CONVFAC = json_serialize(convfac)
          dbinfo.DU_REF = json_serialize(du_ref)
        endelse        
      endif

      ;; We need to take some extra care with the STATE keyword for
      ;; CHROMIS data.
      state = strtrim(fxpar(h, 'STATE'))      
      
      if instrume eq 'CHROMIS' then begin
        
        dbinfo[ifile].FIRST_FRAME = framenumbers[0] ;red_collapserange(framenumbers, ld='', rd='')
        dbinfo[ifile].DATE_BEGS = json_serialize(date_beg)
        
        ;STATISTIC
        dat = red_readdata(files[ifile])
        Nframes = n_elements(framenumbers)
        mn = intarr(Nframes)
        mx = intarr(Nframes)
        med = intarr(Nframes)
        std = fltarr(Nframes)
        for ifrm=0,Nframes-1 do begin
          mx[ifrm] = max(dat[*,*,ifrm], min=mnm)
          mn[ifrm] = mnm
          med[ifrm] = median(dat[*,*,ifrm])
          std[ifrm] = stddev(dat[*,*,ifrm])
        endfor
        dbinfo[ifile].FRAME_MAX = json_serialize(mx)
        dbinfo[ifile].FRAME_MIN = json_serialize(mn)
        dbinfo[ifile].FRAME_MEDIAN = json_serialize(med)
        dbinfo[ifile].FRAME_STDDEV = json_serialize(std)
          
        ;; For Chromis data we may have to convert prefilter and tuning info in
        ;; wheel+hrz form to something readable.
        ;; Note that the STATE keyword is the NB state, also for WB data.
        case 1 of
          strmatch(state, 'wheel[0-9]*_hrz[0-9]*') : begin
            ;; State matches the wheel+hrz pattern
            fpi_info = stregex(state, 'wheel([0-9]*)_hrz([0-9]*)', /extract, /sub)
            ;; NB prefilter
            wheel = long(fpi_info[1])
            nbpref = nbprefs[wheel-1]
            ;; Tuning in digital units
            du = long(fpi_info[2])
            dbinfo[ifile].WHEEL = wheel
            dbinfo[ifile].HRZ = du
            ;; Tuning in [m]
            dlambda = convfac[wheel-1] * (du-du_ref[wheel-1]) 
            lambda_ref_string = string(round(lambda_ref[wheel-1]*1d10), format = '(i04)')
            tuning_string = strtrim(round(dlambda*1d13), 2)
            if strmid(tuning_string, 0, 1) ne '-' then tuning_string = '+'+tuning_string
            state = nbpref + '_' + lambda_ref_string + '_' + tuning_string
          end
          strmatch(state, '[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-][0-9]*') : begin
            ;; State is in the "prefilter_line_[+-]tuning" format.
            ;; Don't change it!
          end
          strmatch(state,'wheel[0-9]*') : begin
            if ~strmatch(state,'hrz[0-9]*') then begin
              ; WB flats
              pos=stregex(state,'wheel([0-9]*)')
              wheel =  long(strmid(state,pos+5,5))
              dbinfo[ifile].WHEEL = wheel
              nbpref = nbprefs[wheel-1]              
              state = nbpref + '_' + nbpref + '_+000' ;fake state
            endif
          end
          else : begin
            ;; No state specified, could be darks. Or WB flats for
            ;; which no simultaneous NB data were collected. 
            state = ''
          end
        endcase

      endif else begin ; only two instruments at the moment     
        ;; CRISP has one frame per file at the present
        if state eq '__+0000' then state='' ; for some reason for Crisp-W darks state is '__+0000'
        if state ne '' then begin
          ; for Crisp WB flats red_readhead construct state without tuning
          ; information and we need it for correct database, 
          ; so let's get it from a filename
          st = stregex(files[ifile], '([0-9][0-9][0-9][0-9])_([+-][0-9]*)', /extract,/subexpr)
          state = st[0]
          dbinfo[ifile].WHEEL = fix(st[1]) ; in fact reference line
          dbinfo[ifile].HRZ = fix(st[2]) ; in fact tuning
        endif
          dbinfo[ifile].FIRST_FRAME = framenumbers
          dbinfo[ifile].DATE_BEGS = date_beg
          fnm = files[ifile]
          pos = stregex(fnm, 'lc[0-9]')
          if pos ne -1 then dbinfo[ifile].LC_STATE = strmid(fnm, pos+2, 1)
          pos = stregex(fnm, 'LP[0-9][0-9][0-9]')
          if pos ne -1 then dbinfo[ifile].LP_STATE = strmid(fnm, pos+2, 3)
          pos = stregex(fnm, 'qw[0-9][0-9][0-9]')
          if pos ne -1 then dbinfo[ifile].QW_STATE = strmid(fnm, pos+2, 3)
          pos = stregex(fnm, 'f[+-][0-9][0-9][0-9]')
          if pos ne -1 then dbinfo[ifile].FOCUS = strmid(fnm, pos+1, 4)        
        ;STATISTIC
        dat = red_readdata(files[ifile])
        dbinfo[ifile].FRAME_MAX = string(max(dat,min=mn))
        dbinfo[ifile].FRAME_MIN = string(mn)
        dbinfo[ifile].FRAME_MEDIAN = string(median(dat))
        dbinfo[ifile].FRAME_STDDEV = string(stddev(dat))
        
      endelse

      dbinfo[ifile].STATE = state

      ;; From the header
      value = fxpar(h, 'BITPIX'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].BITPIX   = value
      value = fxpar(h, 'CADENCE' , count=cnt) & if cnt eq 1 then dbinfo[ifile].CADAVG   = value ; Use the SOLARNET CADAVG keyword
      value = fxpar(h, 'DATE'    , count=cnt) & if cnt eq 1 then dbinfo[ifile].DATE     = value 
      gain = float(fxpar(h, 'DETGAIN' , count=cnt)) & if cnt eq 1 then dbinfo[ifile].DETGAIN  = gain 
      value = fxpar(h, 'DETOFFS' , count=cnt) & if cnt eq 1 then dbinfo[ifile].DETOFFS  = value 
      value = fxpar(h, 'FILTER1' , count=cnt) & if cnt eq 1 then dbinfo[ifile].FILTER1  = value        
;      value = fxpar(h, 'NAXIS'   , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS    = value 
      value = fxpar(h, 'NAXIS1'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS1   = value           
      value = fxpar(h, 'NAXIS2'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS2   = value 
      value = fxpar(h, 'NAXIS3'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS3   = value 
      value = fxpar(h, 'SCANNUM' , count=cnt) & if cnt eq 1 then dbinfo[ifile].SCANNUM  = value 
      value = fxpar(h, 'SOLARNET', count=cnt) & if cnt eq 1 then dbinfo[ifile].SOLARNET = value 
      value = fxpar(h, 'WAVEBAND', count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEBAND = value 
      value = fxpar(h, 'WAVELNTH', count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVELNTH = value 
      value = fxpar(h, 'WAVEMAX' , count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEMAX  = value 
      value = fxpar(h, 'WAVEMIN' , count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEMIN  = value 
      value = fxpar(h, 'WAVEUNIT', count=cnt) & if cnt eq 1 then dbinfo[ifile].WAVEUNIT = value 
      xpos = float(fxpar(h, 'XPOSURE' , count=cnt)) & if cnt eq 1 then dbinfo[ifile].XPOSURE  = xpos
      value = fxpar(h, 'DETTEMP' , count=cnt) & if cnt eq 1 then dbinfo[ifile].DETTEMP  = value

      ;; We need to find unique exposure/gain pairs for CHROMIS
      ;; (For CHROMIS they can change during a scan)
      if instrume eq 'CHROMIS' then begin
        if ifile eq 0 then begin
          cf = fltarr(1,2)
          cf[0,0] = xpos
          cf[0,1] = gain
          red_append,config,cf
        endif else begin
          sz = size(config)
          cnt = 0
          for ll=0, sz[1]-1 do begin
            if config[ll,0] eq xpos and config[ll,1] eq gain then cnt += 1
          endfor
          if cnt eq 0 then begin
            cf[0,0] = xpos
            cf[0,1] = gain
            red_append, config,cf
          endif
        endelse
      endif

    endfor                      ; ifile

    ;; Check for complete scans only
    scans = dbinfo[uniq(dbinfo.scannum, sort(dbinfo.scannum))].scannum
    Nscans = n_elements(scans)
    f_scan = lonarr(Nscans)
    for iscan = 0L, Nscans-1 do $
       f_scan[iscan] = n_elements(where(dbinfo.scannum eq scans[iscan]))
    mask = replicate(1B, Nfiles)
    for iscan = 1L, Nscans-1 do begin
      if f_scan[iscan]-f_scan[0] lt 0 then begin
        print, inam + ' : WARNING : ' + dbinfo[0].camera + ': Incomplete scan nr ' + strtrim(scans[iscan], 2)
        print, inam + '             only ' + strtrim(f_scan[iscan], 2) + ' of ' $
               + strtrim(f_scan(0), 2) + ' files.  Skipping it'
        mask[where(dbinfo.scannum EQ scans[iscan])] = 0
      endif
    endfor                  ; iscan
    idx = where(mask)
    dbinfo = dbinfo[idx]

    ;; Send the array of structs as input to a command that knows what
    ;; info should go into what table, and writes it there.   
    if keyword_set(debug) then save, dbinfo, filename = date + '_' + timestamp + '_' + camera + '.sav'
    if instrume eq 'CHROMIS' then $
       chromis_rawfile2db, dbinfo, config=config, debug=debug $
    else $
       crisp_rawfile2db, dbinfo, debug=debug    

  endfor                        ; icam

end

