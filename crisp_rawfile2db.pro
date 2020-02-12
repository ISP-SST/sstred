; docformat = 'rst'

;+
; Write metadata from an array of structures made in red_rawdir2db procedure into the database.
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
;   2019-10-07 : OA. Optimized sql queries
;
;   2020-02-04 : OA. Split red_rawfile2db to crisp_rawfile2d & chromis_rawfile2d
; 
;-
pro crisp_rawfile2db, dbinfo, config=config, debug=debug
  
  inam = red_subprogram(/low, calling = inam1)

  Nfiles =  n_elements(dbinfo)
  if Nfiles eq 0 then begin
    print, inam, 'There is nothing to put in the database. The array of structures is empty.' 
    return
  endif
  
  ;; Make sure database is open
  red_mysql_check, handle

  ; In rawdir2db this procedure is called in cameras loop, so there is only one unique instrument in the 'dbinfo'
  ; as well as only one 'DATE_OBS' and one 'CAMERA'
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
  query = 'INSERT INTO filters (filter1, wavemax, wavemin, wavelnth, waveunit) VALUES ' ;
  ;; if it is darks then skip
  if dbinfo[filt_in[0]].filter1 ne '0' then begin
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

  
  calibration_id='0'            ; no calibration data for CRISP

  ;; Insert or update information in configs table 
  config_id = '0'
  query = 'INSERT INTO configs (sets_id, camera, detgain, xposure, naxis1, naxis2, naxis3, bitpix, cadavg, ' + $
          'detoffs, detectors_id, dir_templt_id, fnm_templt_id) VALUES (' + set_id + ', "' + $
          dbinfo[0].CAMERA + '", ' + strtrim(string(dbinfo[0].DETGAIN),2) + ', ' + $
          strtrim(string(dbinfo[0].XPOSURE),2)+ ', ' + $ 
          strtrim(string(dbinfo[0].NAXIS1),2) + ', ' + strtrim(string(dbinfo[0].NAXIS2),2) + ', ' + $
          strtrim(string(dbinfo[0].NAXIS3),2) + ', ' + strtrim(string(dbinfo[0].BITPIX),2) + ', ' + $
          strtrim(string(dbinfo[0].CADAVG),2) + ', ' + strtrim(string(dbinfo[0].DETOFFS),2) + ', ' + $
          strtrim(string(detectors_id),2) + ', ' + strtrim(string(dir_templt_id),2) + ', ' + strtrim(string(fnm_templt_id),2) + $
          ') ON DUPLICATE KEY UPDATE detgain = VALUES(detgain), xposure = VALUES(xposure), ' + $
          'naxis1 = VALUES(naxis1), naxis2 = VALUES(naxis2), naxis3 = VALUES(naxis3), bitpix = VALUES(bitpix), ' + $ 
          'cadavg = VALUES(cadavg), detoffs = VALUES(detoffs), detectors_id = VALUES(detectors_id), ' + $
          'dir_templt_id = VALUES(dir_templt_id), fnm_templt_id = VALUES(fnm_templt_id);'
  query += 'SELECT LAST_INSERT_ID();'
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  config_id=ans[1]
  if config_id eq '0' then begin
    query = 'SELECT id FROM configs WHERE camera = "' + dbinfo[0].CAMERA + $
            '" AND sets_id = ' + set_id + ';'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    config_id=ans[1]
  endif
  
                                ;find number of scans
  scans = uniq(dbinfo.scannum, sort(dbinfo.scannum))
  Nscans= n_elements(scans)
                                ;find number of states (wavelength positions in a scan)
  states = uniq(dbinfo.state, sort(dbinfo.state))
  Nstates = n_elements(states)
                                ; for darks we have to use prefilters instead of tuning states
  if Nstates eq 1 then if dbinfo[states].state eq '' then begin
    states = uniq(dbinfo.filter1, sort(dbinfo.filter1))
    Nstates = n_elements(states)
  endif

  for iscan=0, Nscans-1 do begin      
    scannum = dbinfo[scans[iscan]].SCANNUM
    red_progressbar, iscan, Nscans, /predict, strtrim(string(Nscans),2) + ' scans.'

    query = 'BEGIN;'
    for istate=0, Nstates-1 do begin           
      state = dbinfo[states[istate]].state
      pref = dbinfo[states[istate]].filter1
      if state eq '' then $
                                ; for darks we have to use prefilters instead of tuning states
         frms_indx = where(dbinfo.scannum eq scannum and dbinfo.filter1 eq pref) $
      else $          
         frms_indx = where(dbinfo.scannum eq scannum and dbinfo.STATE eq state)
                                ;for CRISP we treat files with the same tuning state as frames in burst
      Nframes = n_elements(frms_indx)
      fst_frm_num = min(dbinfo[frms_indx].first_frame, inn)
      fst_frm_indx = frms_indx[inn]

      line = '0'
      tuning = '0'        
      if state ne '' then begin         
        st = strsplit(state,'_',/extract)                
        line = st[0]
        tuning = st[1]
      endif
      if datatype eq 'flats' and is_wb then line = pref

                                ; change date format to be complient with mariaDB
      ss = strsplit(dbinfo[fst_frm_indx].DATE, 'T',/extract)
      dt = ss[0]
      if n_elements(ss) gt 1 then begin
        zz = strsplit(ss[1],'.',/extract)
        dt += ' ' + zz[0]
      endif
      query += 'INSERT INTO bursts(config_id, line, tuning, scannum, first_frame, dettemp, date, calibration_id, filter1) ' + $
               'VALUES(' + config_id + ', ' + line + ', ' + tuning + ', ' + strtrim(string(scannum), 2) + ', ' + $
               strtrim(string(dbinfo[fst_frm_indx].FIRST_FRAME),2) + ', ' + strtrim(string(dbinfo[fst_frm_indx].DETTEMP),2) + $
               ', "' + dt + '", ' + calibration_id + ', ' + pref + ') ON DUPLICATE KEY UPDATE ' + $
               'line = VALUES(line), tuning = VALUES(tuning), scannum = VALUES(scannum), first_frame = VALUES(first_frame), ' + $
               'dettemp = VALUES(dettemp), date = VALUES(date), filter1 = VALUES(filter1);'
      query += 'SELECT id FROM bursts WHERE first_frame = ' + strtrim(string(dbinfo[fst_frm_indx].FIRST_FRAME),2) + $
               ' AND config_id = ' + config_id + ' INTO @burst_id;'

      if datatype ne 'polcal' then begin

        query += 'INSERT INTO crisp_frames (burst_id, framenum, date_begs, beg_frac, max, min, median, stddev, lc_state) VALUES '
        for ifrm=0, Nframes -1 do begin
          frm_indx = frms_indx[ifrm]
          ss = strsplit(dbinfo[frm_indx].DATE_BEGS,'T',/extract)
          zz = strsplit(ss[1],'.',/extract)
          if n_elements(zz) eq 1 then begin ; a rare ocasion of datastamp without fraction of second
            beg_time = zz
            beg_frac = '0.'
          endif else begin
            beg_time = zz[0]
            beg_frac = '.'+zz[1]
          endelse
          ;;here first frame is a number of crisp file in the scan sequence              
          query += '( @burst_id, ' + strtrim(string(dbinfo[frm_indx].FIRST_FRAME),2) + ', "' + beg_time + '", ' + beg_frac +  $
                   ', ' + strtrim(string(dbinfo[frm_indx].FRAME_MAX),2) + ', ' + strtrim(string(dbinfo[frm_indx].FRAME_MIN),2) + ', ' + $
                   strtrim(string(dbinfo[frm_indx].FRAME_MEDIAN),2) + ', ' + strtrim(string(dbinfo[frm_indx].FRAME_STDDEV),2) + ', ' + $
                   dbinfo[frm_indx].LC_STATE + '),'
        endfor                  ;ifrm
        ll = strlen(query)
        query = strmid(query, 0, ll-1)
        query += ' ON DUPLICATE KEY UPDATE framenum = VALUES(framenum), date_begs = VALUES(date_begs), ' + $
                 'beg_frac = VALUES(beg_frac), max = VALUES(max), min = VALUES(min), median = VALUES(median), ' + $
                 'stddev = VALUES(stddev), lc_state = VALUES(lc_state);'          
      endif else begin

        ;; THIS IS FOR CRISP POLCAL DATA
        query += 'INSERT INTO crisp_polcal_frames (burst_id, framenum, date_begs, beg_frac, max, min, median, stddev, ' + $
                 'lc_state, qw_state, lp_state) VALUES '
        for ifrm=0, Nframes -1 do begin
          frm_indx = frms_indx[ifrm]
          ss = strsplit(dbinfo[frm_indx].DATE_BEGS,'T',/extract)
          zz = strsplit(ss[1],'.',/extract)
          if n_elements(zz) eq 1 then begin ; a rare ocasion of datastamp without fraction of second
            beg_time = zz
            beg_frac = '0.'
          endif else begin
            beg_time = zz[0]
            beg_frac = '.'+zz[1]
          endelse
          ;;here first frame is a number of crisp file in the scan sequence
          query += '(@burst_id, ' + strtrim(string(dbinfo[frm_indx].FIRST_FRAME),2) + ', "' + beg_time + '", ' + beg_frac + $
                   ', ' + strtrim(string(dbinfo[frm_indx].FRAME_MAX),2) + ', ' + strtrim(string(dbinfo[frm_indx].FRAME_MIN),2) + ', ' + $
                   strtrim(string(dbinfo[frm_indx].FRAME_MEDIAN),2) + ', ' + strtrim(string(dbinfo[frm_indx].FRAME_STDDEV),2) + ', ' + $
                   dbinfo[frm_indx].LC_STATE + ', ' + dbinfo[frm_indx].QW_STATE + ', ' + dbinfo[frm_indx].LP_STATE + '),'
        endfor                  ;ifrm
        ll = strlen(query)
        query = strmid(query, 0, ll-1)
        query += ' ON DUPLICATE KEY UPDATE framenum = VALUES(framenum), date_begs = VALUES(date_begs), ' + $
                 'beg_frac = VALUES(beg_frac), max = VALUES(max), min = VALUES(min), median = VALUES(median), ' + $
                 'stddev = VALUES(stddev), lc_state = VALUES(lc_state), qw_state = VALUES(qw_state), lp_state = VALUES(lp_state); '
      endelse                   ;polcal

    endfor                      ; istate

    ;; 'BEGIN; ... COMMIT;' -- this is to make bulk query safer
    query += 'COMMIT;'      
    red_mysql_cmd, handle, query, ans, nl, debug=debug
  endfor                        ; iscan   

  print,'-----------------------------------------------------------------'
  free_lun,handle
end
