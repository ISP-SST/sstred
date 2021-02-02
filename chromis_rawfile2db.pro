; docformat = 'rst'

;+
; Write metadata from an array of structures made in red_rawdir2db
; procedure into the database. 
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
;   dbinfo, in, type=structarr
; 
;     An array of structs, one per file, with metadata info.
; 
; 
; :Keywords:
; 
;   debug, in, optional, type=boolean
;   
;     Set to have verbose output from red_mysql_cmd
; 
; 
; :History:
; 
;   2018-07-02 : MGL. First version.
;
;   2019-05-28 : OA. Second version.
;
;   2019-10-07 : OA. Optimized sql queries.
;
;   2020-02-04 : OA. Split red_rawfile2db to crisp_rawfile2d &
;                chromis_rawfile2d. Added loop for different
;                exposure/gain settings during a scan which are
;                transfered in config array.
;                               
;-
pro chromis_rawfile2db, dbinfo, config, debug=debug
  
  inam = red_subprogram(/low, calling = inam1)

  Nfiles =  n_elements(dbinfo)
  if Nfiles eq 0 then begin
    print, inam, 'There is nothing to put in the database. The array of structures is empty.' 
    return
  endif
  
  ;; Make sure database is open
  red_mysql_check, handle

  ;; In rawdir2db this procedure is called in cameras loop, so there
  ;; is only one unique instrument in the 'dbinfo' as well as only one
  ;; 'DATE_OBS' and one 'CAMERA'
  instrument = dbinfo[0].INSTRUME
  camera =  dbinfo[0].CAMERA
  is_wb = strmatch(camera, '*-[WD]')
  datatype = dbinfo[0].DATATYPE
  tab = string(9B)  

  d = strsplit(dbinfo[0].DATE_OBS,' ',/extract)
  if instrument eq 'CHROMIS' then begin      
    query = 'SELECT id, filter1 FROM calibrations WHERE date = "' + d[0] + '";'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    if nl eq 1 then begin
      ;; insert data in calibrations table for the date
      nbprefs = json_parse(dbinfo[0].NB_PREFS,/TOARRAY)
      lambda_ref = json_parse(dbinfo[0].LAMBDA_REF,/TOARRAY)
      convfac = json_parse(dbinfo[0].CONVFAC,/TOARRAY)
      du_ref = json_parse(dbinfo[0].DU_REF,/TOARRAY)
      nel = n_elements(nbprefs)
      query = 'INSERT INTO calibrations (filter1, convfac, du_ref, date, lambda_ref) VALUES '
      for jj=0,nel-1  do $
        query += '(' + nbprefs[jj] + ', ' + strtrim(string(convfac[jj]),2) + ', ' + $
        strtrim(string(du_ref[jj]),2) + ', "' + d[0] + '", ' + strtrim(string(lambda_ref[jj]),2) + '),'
      ll = strlen(query)
      query = strmid(query, 0, ll-1) + ';'
      red_mysql_cmd, handle, query, ans, nl, debug=debug ;insert
      query = 'SELECT id, filter1 FROM calibrations WHERE date = "' + d[0] + '";'
      red_mysql_cmd, handle, query, ans, nl, debug=debug ; get id, filter1
    endif
    calib = strarr(nl-1,2)
    for jj = 1, nl -1 do $
      calib[jj-1,*] = strsplit(ans[jj],tab,/extract,/preserve_null) ; we need calibration_id for bursts table
  endif

  print,'---------------------------------------------------------------------------------'
  print, inam+' : Write DB info for ' + dbinfo[0].DATE_OBS + '  ' + dbinfo[0].CAMERA

;    observer=''
;    if dbinfo[0].OBSERVER eq '' and keyword_set(description) then begin
;       print, 'Enter observer name for ', instrument, '  ', dbinfo[0].DATE_OBS, ' dataset.'
;       read, observer
;    endif else begin 
;        observer = dbinfo[0].OBSERVER
;    endelse
;    desc=''
;    if keyword_set(description) then begin
;      print, 'Enter description for ', instrument, '  ', dbinfo[0].DATE_OBS, ' dataset.'
;      read, desc
;    endif

  
  query = 'INSERT INTO datasets (date_obs, instrument, data_type, solarnet) VALUES ("'+ $ ;observer, description,
          dbinfo[0].DATE_OBS + '", "'+ instrument + '", "' + $                            ;observer + '", "' + desc + '", "' + 
          datatype + '", ' + strtrim(string(dbinfo[0].solarnet),2) + ')' + $
          'ON DUPLICATE KEY UPDATE ' + $ ;observer = ' +  VALUES(observer), description = VALUES(description), 
          'data_type = VALUES(data_type), solarnet = VALUES(solarnet);'
  query += 'SELECT LAST_INSERT_ID();'
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  set_id=ans[1]
  if set_id eq 0 then begin
    query='SELECT id FROM datasets WHERE date_obs = "' + dbinfo[0].DATE_OBS + '" AND instrument = "' + instrument + '";'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    set_id=ans[1]
  endif
            

  ;; collect IDs for data from different tables to be put in configs table
  query = "SELECT id from  dirname_templates WHERE template = '" + dbinfo[0].DIR_TEMPLATE + "';"
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  if nl eq 1 then begin
    query = "INSERT INTO dirname_templates (template) VALUES('" + dbinfo[0].DIR_TEMPLATE + "');"
    query += 'SELECT LAST_INSERT_ID();'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    dir_templt_id=ans[1]
  endif else dir_templt_id=ans[1]

  query = "SELECT id from  filename_templates WHERE template = '" + dbinfo[0].FNM_TEMPLATE + "';"
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  if nl eq 1 then begin
    query = "INSERT INTO filename_templates (template) VALUES('" + dbinfo[0].FNM_TEMPLATE + "');"
    query += 'SELECT LAST_INSERT_ID();'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    fnm_templt_id=ans[1]
  endif else fnm_templt_id=ans[1]        

  query = 'SELECT id, detfirm FROM detectors WHERE name = "' + dbinfo[0].DETECTOR + '";'
  detectors_id = '0'
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  if nl ne 1 then begin
    ;; Insert new row in the table for the detector with new firmware if it doesn't exist
    ;; Copy column values from the last row for the detector    
    dets = strsplit(ans,tab,/extract,/preserve_null)
    f_indx = where(dets[1] eq dbinfo[0].DETFIRM)          
    if f_indx eq -1 then begin
      id = ans[0,nl-1] ; last id with the current detector
      query = 'INSERT INTO detectors (manufacturer, model, xsize, ysize, pixel_pitch, serial_number, name, detfirm) ' + $
        'SELECT manufacturer, model, xsize, ysize, pixel_pitch, serial_number, name, "' + dbinfo[0].DETFIRM + $
        '" FROM detectors WHERE id = ' + id + ';'
      query += 'SELECT LAST_INSERT_ID();'
      red_mysql_cmd, handle, query, ans, nl, debug=debug
      detectors_id = ans[1]
    endif $ ;; take id from the detectors table for the detector with current firmware
      else detectors_id = dets[f_indx,0]            
  endif
  
  scan0 = where(dbinfo.scannum eq 0)  
  ;; Insert information for prefilter if it's not there.
  ;; We use one set of filters at the moment and don't have to
  ;; distinguish similar ones. May be it will change...
  filt_in = uniq(dbinfo[scan0].FILTER1, sort(dbinfo[scan0].FILTER1))
  ;; if it is darks then skip
  if dbinfo[filt_in[0]].WAVEMAX ne 0. then begin
    query = 'INSERT INTO filters (filter1, wavemax, wavemin, wavelnth, waveunit) VALUES ' ;
    for ifilt = 0, n_elements(filt_in)-1 do begin
      fin = filt_in[ifilt]
      query += '(' + dbinfo[fin].FILTER1 + $
               ', ' + strtrim(string(dbinfo[fin].WAVEMAX),2) + ', ' + strtrim(string(dbinfo[fin].WAVEMIN),2) + ', ' + $
               strtrim(string(dbinfo[fin].WAVELNTH),2) + ', ' + strtrim(string(dbinfo[fin].WAVEUNIT),2) + '),'    
    endfor
    ll = strlen(query)
    query = strmid(query, 0, ll-1)
    query += ' ON DUPLICATE KEY UPDATE wavemax = VALUES(wavemax), wavemin = VALUES(wavemin), wavelnth = VALUES(wavelnth);'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
  endif


  Nconf = size(config)

  for iconf = 0,Nconf[1]-1 do begin
    ;; Insert or update information in configs table 
    config_id = '0'
    xpos = config[iconf,0]
    gain = config[iconf,1]
    query = 'INSERT INTO configs (sets_id, camera, detgain, xposure, naxis1, naxis2, naxis3, bitpix, cadavg, ' + $
            'detoffs, detectors_id, dir_templt_id, fnm_templt_id) VALUES (' + set_id + ', "' + $
            dbinfo[0].CAMERA + '", ' + strtrim(string(gain),2) + ', ' + strtrim(string(xpos),2)+ ', ' + $
            strtrim(string(dbinfo[0].NAXIS1),2) + ', ' + strtrim(string(dbinfo[0].NAXIS2),2) + ', ' + $
            strtrim(string(dbinfo[0].NAXIS3),2) + ', ' + strtrim(string(dbinfo[0].BITPIX),2) + ', ' + $
            strtrim(string(dbinfo[0].CADAVG),2) + ', ' + strtrim(string(dbinfo[0].DETOFFS),2) + ', ' + $
            strtrim(string(detectors_id),2) + ', ' + strtrim(string(dir_templt_id),2) + ', ' + $
            strtrim(string(fnm_templt_id),2) + $
            ') ON DUPLICATE KEY UPDATE detgain = VALUES(detgain), xposure = VALUES(xposure), ' + $
            'naxis1 = VALUES(naxis1), naxis2 = VALUES(naxis2), naxis3 = VALUES(naxis3), bitpix = VALUES(bitpix), ' + $ 
            'cadavg = VALUES(cadavg), detoffs = VALUES(detoffs), detectors_id = VALUES(detectors_id), ' + $
            'dir_templt_id = VALUES(dir_templt_id), fnm_templt_id = VALUES(fnm_templt_id);'
    ;query += 'SELECT LAST_INSERT_ID();'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    ;config_id=ans[1]
    ;if fix(config_id) > 100 then stop
                                ;if config_id eq '0' then begin

    ;;For some reason there were occasions when
    ;;LAST_INSERT_ID()returned weird values for config_id. So,
    ;;let's keep to query it explicitly
      query = 'SELECT id FROM configs WHERE camera = "' + dbinfo[0].CAMERA + $
              '" AND sets_id = ' + set_id + ' AND ABS(xposure - ' + strtrim(string(config[iconf,0]),2) + $
              ')<0.0001 AND ABS(detgain - ' + strtrim(string(config[iconf,1]),2) + ')<0.0001;'
      red_mysql_cmd, handle, query, ans, nl, debug=debug
      config_id=ans[1]
    ;endif

    ;; we need to check only exposure and gain as the camera is the same
    files_in = where(dbinfo.xposure eq xpos and dbinfo.detgain eq gain)
    Nfiles = n_elements(files_in)
    for ifl = 0, Nfiles-1 do begin
      red_progressbar, ifl, Nfiles, /predict, strtrim(string(Nfiles),2) + ' bursts.'        
      ifile = files_in[ifl]

      st = strsplit(dbinfo[ifile].STATE,'_',/extract)
      in_cal = where(strmatch(calib[*,1],st[0]) eq 1)
      ;; We need calibration_id for burst table
      if in_cal ne -1 then calibration_id = calib[in_cal,0] else calibration_id = '0'

      ;; change date format to be complient with mariaDB
      ss = strsplit(dbinfo[ifile].DATE, 'T',/extract)
      dt = ss[0]
      zz = strsplit(ss[1],'.',/extract)
      dt += ' ' + zz[0]
      ;;We put 'wheel' and 'hrz' values in 'line' and 'tuning' columns 
      ;;of the bursts table to reconstruct filenames and FP state
      ;;after.
      query = 'BEGIN;'
      query += 'INSERT INTO bursts(config_id, filter1, line, tuning, calibration_id, scannum, first_frame, dettemp, date) VALUES(' + $
               strtrim(string(config_id),2) + ', ' + dbinfo[ifile].FILTER1 + ', ' + strtrim(string(dbinfo[ifile].WHEEL),2) + ', ' + $
               strtrim(string(dbinfo[ifile].HRZ),2) + ', ' + calibration_id + ', ' + $
               strtrim(string(dbinfo[ifile].SCANNUM), 2) + ', ' + strtrim(string(dbinfo[ifile].FIRST_FRAME),2) + ', ' + $
               strtrim(string(dbinfo[ifile].DETTEMP),2) + ', "' + dt + '") ON DUPLICATE KEY UPDATE ' + $
               'filter1 = VALUES(filter1), line = VALUES(line), tuning = VALUES(tuning), calibration_id = VALUES(calibration_id), ' + $
               'scannum = VALUES(scannum), dettemp = VALUES(dettemp), date = VALUES(date);'
      query += 'SELECT id FROM bursts WHERE first_frame = ' + strtrim(string(dbinfo[ifile].FIRST_FRAME),2) + $
               ' AND config_id = ' + config_id + ' INTO @burst_id;'     

      Nframe = dbinfo[ifile].NAXIS3
      date_begs = json_parse(dbinfo[ifile].DATE_BEGS,/TOARRAY)
      mx = json_parse(dbinfo[ifile].FRAME_MAX,/TOARRAY)
      mn = json_parse(dbinfo[ifile].FRAME_MIN,/TOARRAY)
      median = json_parse(dbinfo[ifile].FRAME_MEDIAN,/TOARRAY)
      stddev = json_parse(dbinfo[ifile].FRAME_STDDEV,/TOARRAY)   

      query += 'INSERT INTO chromis_frames (burst_id, framenum, date_begs, beg_frac, max, min, median, stddev) VALUES '
      for iframe = 0, Nframe-1 do begin
        frame_num = dbinfo[ifile].FIRST_FRAME + iframe ;
        ss = strsplit(date_begs[iframe],'T',/extract)
        zz = strsplit(ss[1],'.',/extract)
        if n_elements(zz) eq 1 then begin ; a rare ocasion of datastamp without fraction of second
          beg_time = zz
          beg_frac = '0.'
        endif else begin
          beg_time = zz[0]
          beg_frac = '.'+zz[1]
        endelse
        query += '( @burst_id, ' + strtrim(string(frame_num),2) + ', "' + beg_time + '", ' + $
                 beg_frac + ', ' + strtrim(string(mx[iframe]),2) + ', ' + strtrim(string(mn[iframe]),2) + $
                 ', ' + strtrim(string(median[iframe]),2) + ', ' + strtrim(string(stddev[iframe]),2) + '),'
      endfor                    ; iframe
      ll = strlen(query)
      query = strmid(query, 0, ll-1)
      query += ' ON DUPLICATE KEY UPDATE framenum = VALUES(framenum), date_begs = VALUES(date_begs), ' + $
               'beg_frac = VALUES(beg_frac), max = VALUES(max), min = VALUES(min), median = VALUES(median), stddev = VALUES(stddev);'
      ;; 'BEGIN; ... COMMIT;' -- this is to make bulc query safer
      query += 'COMMIT;'
      red_mysql_cmd, handle, query, ans, nl, debug=debug

    endfor                      ; ifile
  endfor                        ; iconf
    

  print,'-----------------------------------------------------------------'
  free_lun,handle
end
