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
;   2019-10-07 : OA. Optimized sql queries
;
;   2020-02-04 : OA. Split red_rawfile2db to crisp_rawfile2d &
;                chromis_rawfile2d
;
;   2020-02-17 : OA. Added loop for different prefilter/filename
;                template pairs which are transfered in config array.
;
;-
pro crisp2_rawfile2db, dbinfo, config, debug=debug
  
  inam = red_subprogram(/low, calling = inam1)

  Nfiles =  n_elements(dbinfo)
  if Nfiles eq 0 then begin
    print, inam, 'There is nothing to put in the database. The array of structures is empty.' 
    return
  endif

  ;; Make sure database is open
  year = strmid(dbinfo[0].DATE_OBS,0,4)
  if year ge 2023 then db = 'sst_db_' + year else db = 'sst_db'  
  red_mysql_check, handle, database=db

  ;; In rawdir2db this procedure is called in cameras loop, so there
  ;; is only one unique instrument in the 'dbinfo' as well as only one
  ;; 'DATE_OBS' and one 'CAMERA'
  instrument = dbinfo[0].INSTRUME
  camera =  dbinfo[0].CAMERA
  is_wb = strmatch(camera, '*-[WD]')
  datatype = dbinfo[0].DATATYPE
  tab = string(9B)  

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

  query='SELECT id FROM datasets WHERE date_obs = "' + dbinfo[0].DATE_OBS + '" AND instrument = "CRISP";'
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  if nl eq 1 then begin
    query = 'INSERT INTO datasets (date_obs, instrument, data_type, solarnet) VALUES ("'+ $ ;observer, description,
            dbinfo[0].DATE_OBS + '", "CRISP", "' + $                            ;observer + '", "' + desc + '", "' + 
            datatype + '", ' + strtrim(string(dbinfo[0].solarnet),2) + ');'
    query += 'SELECT LAST_INSERT_ID();'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
  endif
  set_id=ans[1]            

  ;; collect IDs for data from different tables to be put in configs table
  query = "SELECT id from  dirname_templates WHERE template = '" + dbinfo[0].DIR_TEMPLATE + "';"
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  if nl eq 1 then begin
    query = "INSERT INTO dirname_templates (template) VALUES('" + dbinfo[0].DIR_TEMPLATE + "');"
    query += 'SELECT LAST_INSERT_ID();'
    red_mysql_cmd, handle, query, ans, nl, debug=debug    
  endif
  dir_templt_id=ans[1]

  query = 'SELECT id, detfirm FROM detectors WHERE name = "' + dbinfo[0].DETECTOR + '";'
  detectors_id = '0'
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  if nl ne 1 then begin
    ;; Insert new row in the table for the detector with new firmware if it doesn't exist
    ;; Copy column values from the last row for the detector
    dets = strarr(nl-1,2)    
    for jj = 1, nl -1 do $
      dets[jj-1,*] = strsplit(ans[jj],tab,/extract,/preserve_null)
    f_indx = where(dets[*,1] eq dbinfo[0].DETFIRM)          
    if f_indx eq -1 then begin
      id = dets[-1,0] ; last id with the current detector
      query = 'INSERT INTO detectors (manufacturer, model, xsize, ysize, pixel_pitch, serial_number, name, detfirm) ' + $
        'SELECT manufacturer, model, xsize, ysize, pixel_pitch, serial_number, name, "' + dbinfo[0].DETFIRM + $
        '" FROM detectors WHERE id = ' + id + ';'
      query += 'SELECT LAST_INSERT_ID();'
      red_mysql_cmd, handle, query, ans, nl, debug=debug
      detectors_id = ans[1]
    endif $ ;; take id from the detectors table for the detector with current firmware
      else detectors_id = dets[f_indx,0]            
  endif else stop
  
  scan0 = where(dbinfo.scannum eq 0)  
  ;; Insert information for prefilter if it's not there.
  ;; We use one set of filters at the moment and don't have to
  ;; distinguish similar ones. May be it will change...
  
  ;; filt_in = uniq(dbinfo[scan0].FILTER1, sort(dbinfo[scan0].FILTER1))
  ;; query = 'INSERT INTO filters (filter1, wavemax, wavemin, wavelnth, waveunit) VALUES ' ;
  ;; ;; if it is darks then skip
  ;; if dbinfo[filt_in[0]].wavemax ne 0. then begin
  ;;   for ifilt = 0, n_elements(filt_in)-1 do begin
  ;;     fin = filt_in[ifilt]
  ;;     query += '(' + dbinfo[fin].FILTER1 + $
  ;;              ', ' + strtrim(string(dbinfo[fin].WAVEMAX),2) + ', ' + strtrim(string(dbinfo[fin].WAVEMIN),2) + ', ' + $
  ;;              strtrim(string(dbinfo[fin].WAVELNTH),2) + ', ' + strtrim(string(dbinfo[fin].WAVEUNIT),2) + '),'    
  ;;   endfor
  ;;   ll = strlen(query)
  ;;   query = strmid(query, 0, ll-1)
  ;;   query += ' ON DUPLICATE KEY UPDATE wavemax = VALUES(wavemax), wavemin = VALUES(wavemin), wavelnth = VALUES(wavelnth);'
  ;;   red_mysql_cmd, handle, query, ans, nl, debug=debug
  ;; endif

  query = "SELECT id from  filename_templates WHERE template = '" + dbinfo[0].FNM_TEMPLATE + "';"
  red_mysql_cmd, handle, query, ans, nl, debug=debug
  if nl eq 1 then begin
    query = "INSERT INTO filename_templates (template) VALUES('" + dbinfo[0].FNM_TEMPLATE + "');"
    query += 'SELECT LAST_INSERT_ID();'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    fnm_templt_id=ans[1]
  endif else fnm_templt_id=ans[1]
  
  calibration_id='0'            ; no calibration data for CRISP

  Nconf = size(config)
  for iconf = 0, Nconf[1]-1 do begin

    ;; Insert or update information in configs table
    xpos = config[iconf,0]
    gain = config[iconf,1]
    query = 'SELECT id FROM configs WHERE camera = "' + dbinfo[0].CAMERA + $
            '" AND sets_id = ' + set_id +  $
            ' AND ABS(xposure - ' + strtrim(xpos,2) + $
            ')<0.0001 AND ABS(detgain - ' + strtrim(gain,2) + ')<0.0001;'
    
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    if nl eq 1 then begin
      query = 'INSERT INTO configs (sets_id, camera, detgain, xposure, naxis1, naxis2, naxis3, bitpix, cadavg, ' + $
              'detoffs, detectors_id, dir_templt_id, fnm_templt_id) VALUES (' + set_id + ', "' + $
              dbinfo[0].CAMERA + '", ' + strtrim(string(gain),2) + ', ' + strtrim(string(xpos),2) + ', ' + $ 
              strtrim(string(dbinfo[0].NAXIS1),2) + ', ' + strtrim(string(dbinfo[0].NAXIS2),2) + ', ' + $
              strtrim(string(dbinfo[0].NAXIS3),2) + ', ' + strtrim(string(dbinfo[0].BITPIX),2) + ', ' + $
              strtrim(string(dbinfo[0].CADAVG),2) + ', ' + strtrim(string(dbinfo[0].DETOFFS),2) + ', ' + $
              strtrim(string(detectors_id),2) + ', ' + strtrim(string(dir_templt_id),2) + ', ' + strtrim(string(fnm_templt_id),2) + ');'
      query += 'SELECT LAST_INSERT_ID();'        
      red_mysql_cmd, handle, query, ans, nl, debug=debug
    endif    
    config_id=ans[1]
    
    files_in = where(dbinfo.xposure eq xpos and dbinfo.detgain eq gain)
    Nfiles = n_elements(files_in)    

    for ifl = 0, Nfiles-1 do begin
      ;; CRISP2 files contain frames with same tuning and
      ;; polarization state. So, frames consequent in time are
      ;; spreaded over different files. Still, let's call them bursts.
      red_progressbar, ifl, Nfiles, /predict, strtrim(string(Nfiles),2) + ' bursts.'        
      ifile = files_in[ifl]

      state = dbinfo[ifile].STATE
      lc_state = dbinfo[ifile].LC_STATE      
      pref = dbinfo[ifile].FILTER1
      line = '0'
      tuning = '0'        
      if state ne '' then begin         
        st = strsplit(state,'_',/extract)                
        line = st[0]
        tuning = st[1]
      endif
      if datatype eq 'flats' and is_wb then line = pref

      ;; change date format to be complient with mariaDB
      ss = strsplit(dbinfo[ifile].DATE, 'T',/extract)
      dt = ss[0]
      if n_elements(ss) gt 1 then begin
        zz = strsplit(ss[1],'.',/extract)
        dt += ' ' + zz[0]
      endif
      
      query = 'BEGIN;'
      query += 'INSERT INTO bursts(config_id, line, tuning, scannum, first_frame, dettemp, date, calibration_id, filter1, mosaic_pos) ' + $
                 'VALUES(' + config_id + ', ' + line + ', ' + tuning + ', ' + strtrim(string(dbinfo[ifile].SCANNUM), 2) + ', ' + $
                 strtrim(string(dbinfo[ifile].FIRST_FRAME),2) + ', ' + strtrim(string(dbinfo[ifile].DETTEMP),2) + $
                 ', "' + dt + '", ' + calibration_id + ', ' + pref  + ', ' + dbinfo[ifile].MOSAIC_POS + ') ON DUPLICATE KEY UPDATE ' + $
                 'line = VALUES(line), tuning = VALUES(tuning), scannum = VALUES(scannum), first_frame = VALUES(first_frame), ' + $
                 'dettemp = VALUES(dettemp), date = VALUES(date), filter1 = VALUES(filter1), mosaic_pos = VALUES(mosaic_pos);'
      query += 'SELECT id FROM bursts WHERE first_frame = ' + strtrim(string(dbinfo[ifile].FIRST_FRAME),2) + $
                 ' AND config_id = ' + config_id + ' INTO @burst_id;'

      Nframes = dbinfo[ifile].NAXIS3
      framenumbers = json_parse(dbinfo[ifile].FRAMENUMBERS,/TOARRAY)
      date_begs = json_parse(dbinfo[ifile].DATE_BEGS,/TOARRAY)
      mx = json_parse(dbinfo[ifile].FRAME_MAX,/TOARRAY)
      mn = json_parse(dbinfo[ifile].FRAME_MIN,/TOARRAY)
      median = json_parse(dbinfo[ifile].FRAME_MEDIAN,/TOARRAY)
      stddev = json_parse(dbinfo[ifile].FRAME_STDDEV,/TOARRAY)
        
      if datatype ne 'polcal' then begin

        query += 'INSERT INTO crisp_frames (burst_id, framenum, date_begs, beg_frac, max, min, median, stddev, lc_state) VALUES '
        for ifrm=0, Nframes -1 do begin            

          frame_num = framenumbers[ifrm]          
          if n_elements(mx) ne 1 then begin
            maxx = mx[ifrm]
            minn = mn[ifrm]
            std = stddev[ifrm]
            mdn = median[ifrm]
          endif else begin
            maxx = 0
            minn = 0
            std = 0
            mdn = 0
          endelse

          ss = strsplit(date_begs[ifrm],'T',/extract)
          zz = strsplit(ss[1],'.',/extract)
          if n_elements(zz) eq 1 then begin ; a rare ocasion of datastamp without fraction of second
            beg_time = zz
            beg_frac = '0.'
          endif else begin
            beg_time = zz[0]
            beg_frac = '.'+zz[1]
          endelse
          ;;here first frame is a number of crisp file in the scan sequence              
          query += '( @burst_id, ' + strtrim(frame_num,2) + ', "' + beg_time + '", ' + beg_frac +  $
                     ', ' + strtrim(maxx,2) + ', ' + strtrim(minn,2) + ', ' + $
                     strtrim(mdn,2) + ', ' + strtrim(std,2) + ', ' + $
                     lc_state + '),'
        endfor                ;ifrm

        ll = strlen(query)
        query = strmid(query, 0, ll-1)
        query += ' ON DUPLICATE KEY UPDATE framenum = VALUES(framenum), date_begs = VALUES(date_begs), ' + $
                   'beg_frac = VALUES(beg_frac), max = VALUES(max), min = VALUES(min), median = VALUES(median), ' + $
                   'stddev = VALUES(stddev), lc_state = VALUES(lc_state);'          

      endif else begin

         ;; THIS IS FOR CRISP POLCAL DATA
        qw_state = strtrim(fix(dbinfo[ifile].QW_STATE),2)
        lp_state = strtrim(fix(dbinfo[ifile].LP_STATE),2)
        query += 'INSERT INTO crisp_polcal_frames (burst_id, framenum, date_begs, beg_frac, max, min, median, stddev, ' + $
                   'lc_state, qw_state, lp_state) VALUES '
        for ifrm=0, Nframes -1 do begin            

          frame_num = framenumbers[ifrm]
          if n_elements(mx) ne 1 then begin
            maxx = mx[ifrm]
            minn = mn[ifrm]
            std = stddev[ifrm]
            mdn = median[ifrm]
          endif else begin
            maxx = 0
            minn = 0
            std = 0
            mdn = 0
          endelse
          
          ss = strsplit(date_begs[ifrm],'T',/extract)
          zz = strsplit(ss[1],'.',/extract)
          if n_elements(zz) eq 1 then begin ; a rare ocasion of datastamp without fraction of second
            beg_time = zz
            beg_frac = '0.'
          endif else begin
            beg_time = zz[0]
            beg_frac = '.'+zz[1]
          endelse          
          ;;here first frame is a number of crisp file in the scan sequence
          query += '(@burst_id, ' + strtrim(frame_num,2) + ', "' + beg_time + '", ' + beg_frac + $
                     ', ' + strtrim(maxx,2) + ', ' + strtrim(minn,2) + ', ' + $
                     strtrim(mdn,2) + ', ' + strtrim(std,2) + ', ' + $
                     lc_state + ', ' + qw_state + ', ' + lp_state + '),'
        endfor                ;ifrm
        ll = strlen(query)
        query = strmid(query, 0, ll-1)
        query += ' ON DUPLICATE KEY UPDATE framenum = VALUES(framenum), date_begs = VALUES(date_begs), ' + $
                   'beg_frac = VALUES(beg_frac), max = VALUES(max), min = VALUES(min), median = VALUES(median), ' + $
                   'stddev = VALUES(stddev), lc_state = VALUES(lc_state), qw_state = VALUES(qw_state), lp_state = VALUES(lp_state); '
      endelse                 ;polcal

      ;; 'BEGIN; ... COMMIT;' -- this is to make bulk query safer
      query += 'COMMIT;'      
      red_mysql_cmd, handle, query, ans, nl, debug=debug
    endfor                      ; ifile

  endfor ; iconf

  print,'-----------------------------------------------------------------'
  free_lun,handle
end
