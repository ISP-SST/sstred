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
;   timestamp, in, type=string
; 
;     The timestamp-directory to process.
; 
; 
; :Keywords:
; 
;   cams : in, optional, type=strarr, default=all
;   
;     The cameras to include.
; 
; 
; :History:
; 
;   2019-01-31 : MGL. First version. 
; 
;-
pro red_rawdir2db, all = all $
                   , date = date $
                   , dir = dir $
                   , instrument = instrument

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
    timestamp = splitdir[-2]
    isodate = red_strreplace(splitdir[-4], '.', '-', n = 2)

    date_obs = isodate+'T'+timestamp
    
    is_wb = strmatch(camera, '*-[WD]')
    
    ;; Get the DATATYPE from the path.
    a = stregex(camdirs[idir], '/CHROMIS-([a-zA-Z]*)/', /extract, /sub)
    if n_elements(a) eq 2 then begin
      
      instrume = 'CHROMIS'

      datatype = a[1]
      case datatype of
        'data'     : datatype = 'science' 
        else : datatype += ' (unknown)'
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
        refinfo = f0(zfile)
        lambda_ref[wheel-1] = refinfo[0]
        du_ref[wheel-1]     = refinfo[1]
        convfac[wheel-1]    = refinfo[2]
      endfor                    ; wheel

      ;; Get the hrz conversion info
      a = stregex(camdirs[idir], '/([12][0-9][0-9][0-9][-.][01][0-9][-.][0-3][0-9])/', /extract, /sub)
      red_download_linedefs, isodate, strjoin(splitdir[0:-2],'/'), './'
      chromis_hrz_zeropoint, './'
      
    endif else begin
      
      instrume = 'CRISP'
      
      datatype = splitdir[-3]
      case datatype of
        'Darks'    : datatype = 'darks'    
        'Flats'    : datatype = 'flats'    
        'Pinholes' : datatype = 'pinholes' 
        'Polcal'   : datatype = 'polcal'   
        'Science'  : datatype = 'science'  
        else : datatype += ' (unknown)'
      endcase 
      
    endelse
    
    ;; Now search for the actual files
    files = file_search(camdirs[idir]+'/*', count = Nfiles)
    
    if Nfiles eq 0 then continue
    

    
    ;; Make struct for info from header, that should go into the
    ;; database. Some info from the path is set right away.

    strct  = {BITPIX           : 0B          $  
              , CADAVG         : 0.          $
              , CAMERA         : camera      $ 
              , DATATYPE       : datatype    $
              , DATE           : ''          $
              , DATE_BEGS      : ''          $
              , DATE_OBS       : date_obs    $
              , DETECTOR       : ''          $
              , DETFIRM        : ''          $
              , DETGAIN        : 0.          $
              , DETMODEL       : ''          $
              , DETOFFS        : 0L          $
              , FILENAME       : ''          $
              , FILTER1        : ''          $
              , FRAMENUMS      : ''          $
              , INSTRUME       : instrume    $
              , IS_WB          : is_wb       $
              , NAXIS          : 0B          $
              , NAXIS1         : 0L          $          
              , NAXIS2         : 0L          $
              , NAXIS3         : 0L          $
              , OBJECT         : ''          $
              , OBSERVER       : ''          $
              , OBS_HDU        : 0B          $
              , ORIGIN         : ''          $
              , SCANNUM        : 0L          $
              , SOLARNET       : 0.          $
              , STATE          : ''          $
              , TELESCOP       : ''          $
              , WAVEBAND       : ''          $
              , WAVELNTH       : 0.          $
              , WAVEMAX        : 0.          $
              , WAVEMIN        : 0.          $
              , WAVEUNIT       : 0           $
              , XPOSURE        : 0.          $
             }

    ;; Make array of such structs and set the values known without
    ;; reading the headers
    dbinfo = replicate(strct, Nfiles)
    dbinfo.FILENAME = file_basename(files)

    ;; Fill in other values in the array of structs from the header:
    for ifile = 0, Nfiles-1 do begin

      red_progressbar, ifile, Nfiles, /predict, 'Get DB values for '+camdirs[idir]
      
      h = red_readhead(files[ifile], date_beg=date_beg, framenumbers=framenumbers)

      if ifile eq 0 then begin
        ;; Set the keywords that absolutely do not change from file to
        ;; file within the directory.
        value = fxpar(h, 'DETECTOR', count=cnt) & if cnt eq 1 then dbinfo.DETECTOR = value 
        value = fxpar(h, 'DETFIRM' , count=cnt) & if cnt eq 1 then dbinfo.DETFIRM  = value 
        value = fxpar(h, 'DETMODEL', count=cnt) & if cnt eq 1 then dbinfo.DETMODEL = value 
        value = fxpar(h, 'OBJECT'  , count=cnt) & if cnt eq 1 then dbinfo.OBJECT   = value 
        value = fxpar(h, 'OBSERVER', count=cnt) & if cnt eq 1 then dbinfo.OBSERVER = value 
        value = fxpar(h, 'OBS_HDU' , count=cnt) & if cnt eq 1 then dbinfo.OBS_HDU  = value 
        value = fxpar(h, 'ORIGIN'  , count=cnt) & if cnt eq 1 then dbinfo.ORIGIN   = value 
        value = fxpar(h, 'TELESCOP', count=cnt) & if cnt eq 1 then dbinfo.TELESCOP = value         
      endif
      
      ;; These have one value per frame so we need to turn them into strings:
      dbinfo[ifile].FRAMENUMS = red_collapserange(framenumbers, ld='', rd='')
      dbinfo[ifile].DATE_BEGS = json_serialize(date_beg)

      ;; We need to take some extra care with the STATE keyword for
      ;; CHROMIS data.
      state = fxpar(h, 'STATE')
      if instrume eq 'CHROMIS' then begin
          
        ;; For Chromis data we may have to convert prefilter and tuning info in
        ;; wheel+hrz form to something readable.
        ;; Note that the STATE keyword is the NB state, also for WB data.

        case 1 of
          strmatch(state, 'wheel[0-9]*_hrz[0-9]*') : begin
            ;; State matches the wheel+hrz pattern
            fpi_info = stregex(state, 'wheel([0-9]*)_hrz([0-9]*)', /extract, /sub)
            
            ;; State matches the wheel+hrz pattern
            
            ;; NB prefilter
            wheel = long(fpi_info[1])
            nbpref = nbprefs[wheel-1]
            
            ;; Tuning in digital units
            du = fpi_info[2]
            
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
          else : begin
            ;; No state specified, could be darks. Or WB flats for
            ;; which no simultaneous NB data were collected. 
            state = ''
          end
        endcase

      endif                     ; CHROMIS
      dbinfo[ifile].STATE = state

      ;; From the header
      value = fxpar(h, 'BITPIX'  , count=cnt) & if cnt eq 1 then dbinfo[ifile].BITPIX   = value
      value = fxpar(h, 'CADENCE' , count=cnt) & if cnt eq 1 then dbinfo[ifile].CADAVG   = value ; Use the SOLARNET CADAVG keyword
      value = fxpar(h, 'DATE'    , count=cnt) & if cnt eq 1 then dbinfo[ifile].DATE     = value 
      value = fxpar(h, 'DETGAIN' , count=cnt) & if cnt eq 1 then dbinfo[ifile].DETGAIN  = value 
      value = fxpar(h, 'DETOFFS' , count=cnt) & if cnt eq 1 then dbinfo[ifile].DETOFFS  = value 
      value = fxpar(h, 'FILTER1' , count=cnt) & if cnt eq 1 then dbinfo[ifile].FILTER1  = value        
      value = fxpar(h, 'NAXIS'   , count=cnt) & if cnt eq 1 then dbinfo[ifile].NAXIS    = value 
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
      value = fxpar(h, 'XPOSURE' , count=cnt) & if cnt eq 1 then dbinfo[ifile].XPOSURE  = value 

    endfor                      ; ifile
    
    ;; Send the array of structs as input to a command that knows what
    ;; info should go into what table, and writes it there.
    stop
    
    ;;red_rawfile2db, dbinfo
    
  endfor                        ; icam

  stop
end



;; a = chromisred(/dev)

;; a -> rawdir2db, '09:28:36'

;end                             
