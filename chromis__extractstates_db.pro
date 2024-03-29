; docformat = 'rst'

;+
; Extract states information from the database
;
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;     Oleksii Andriienko, Institute for Solar Physics
; 
; 
; :Returns:
;
; 
; :Params:
; 
;    datasets : in, type=strarr
;   
;      A list of datasets ('YYYY-MM-DD HH:MM:SS) for which to
;      extract the states information. 
;   
;    states : out, type=array(struct)
;
;        An array of structs, containing (partially filled) state
;        information. 
; 
; :History:
; 
;   2019-07-10 : OA. Created. (Derived from chromis__extractstates and
;                red_readhead_db)
;
;   2019-07-11 : MGL. Call like extractstates.
;
;   2020-03-11 : OA. Added check for unique cameras when in the
;                procedure call filenames(strings) are used instead of
;                datasets. Added keyword 'cam'.
;
;-
pro chromis::extractstates_db, strings, states, datasets = datasets, cam = cam
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)          
  
  use_strings = n_elements(datasets) eq 0
  
  ;; clear states and try to get them from the cache
  undefine,states 
  if use_strings then begin
    Nstr = n_elements(strings)
    strings = red_strreplace(strings,'//','/',n=100) 
    ; cache can be slow comparing to the database at least when BeeGFS is fatigue
    if ~keyword_set(force) then begin      
      for ifile=0,Nstr-1 do begin
        this_cache = rdx_cacheget(strings[ifile], count = cnt)   
        if cnt gt 0 then begin
          red_progressbar, ifile, Nstr, 'Extract state info from cache', /predict
          red_append,states,this_cache.state
        endif
      endfor
      if n_elements(states) eq Nstr then return else undefine,states
    endif
    timestamps = stregex(strings,'[0-2][0-9]:[0-5][0-9]:[0-5][0-9]', /extract)
    timestamp = timestamps[uniq(timestamps,sort(timestamps))]
    datasets = self.isodate + ' ' + timestamp

    if keyword_set(cam) then ucams = [cam] else begin
      cms = strarr(Nstr)
      for ii = 0, Nstr-1 do begin
        dr = file_dirname(strings[ii])
        cm = (strsplit(dr,'/',/extract))[-1]
        if strmatch(cm,'*nostate*') then cm = (strsplit(cm,'_',/extract))[0]
        cms[ii] = cm
      endfor
      ucams = cms[uniq(cms,sort(cms))]
    endelse
  endif

  instrument = 'CHROMIS'

  year = strmid(self.isodate,0,4)
  if year ge 2023 then db = 'sst_db_' + year else db = 'sst_db'
  red_mysql_check, handle, database=db

  Nsets = n_elements(datasets)
  set_msg = string('Get DB values for ', strtrim(string(Nsets),2), ' datasets.')

  split_dir = strsplit(self.root_dir,'/',/extract)
  dir_root = '/' + strjoin(split_dir[0:-2],'/') + '/'

  for iset = 0, Nsets-1 do begin
    state = {CHROMIS_STATE}
    red_progressbar, iset, Nsets, /predict, set_msg     
    split_set = strsplit(datasets[iset], ' ', /extract)
    date_time = datasets[iset]
    
    ;; we need Y,...,sec to generate the directory name
    split_date = strsplit(split_set[0],'-',/extract)
    Y = fix(split_date[0])
    M = fix(split_date[1])
    D = fix(split_date[2])

    split_time = strsplit(split_set[1],':',/extract)
    hr = fix(split_time[0])
    min = fix(split_time[1])
    sec = fix(split_time[2])
    
    query='SELECT * FROM datasets WHERE date_obs = "' + date_time + '" AND instrument = "' + instrument + '";'
    red_mysql_cmd,handle,query,ans,nl
    if nl eq 1 then begin
      print, inam, ': There is no entry for ', date_time, ' dataset in the database.'
      print, "Let's run red_rawdir2db first."
      flat_dirs = *self.flat_dir
      dark_dirs = *self.dark_dir
      data_dirs = *self.data_dirs
      pinh_dirs = *self.pinh_dirs
      dirs = [flat_dirs, dark_dirs, pinh_dirs, data_dirs]
      in = where(strmatch(dirs,'*'+split_set[1]+'*'))
      if n_elements(in) eq 1 then begin
        if in eq -1 then begin
          print, inam, ': There are no directories for ', +split_set[1], ' timestamp.'
          free_lun,handle
          return
        endif
        dir = dirs[in]
      endif else begin          ; a rather accident case of several directories with same timestamp
        print,'For some reason there are several directories with same timestamp. Please, choose one or exit.'
        for ll = 0, n_elements(in)-1 do $
          print,'[',strtrim(ll,2),']   ',dirs[in[ll]]
        print,'[',strtrim(ll,2),']   exit'
        read,choice
        if choice eq ll then exit
        dir = dirs[in[choice]]
      endelse
      red_rawdir2db,dir=dir,/all
      red_mysql_cmd,handle,query,ans,nl
    endif
    ;;parse a result of the query
    tab = string(9B)
    set = strsplit(ans[1],tab,/extract,/preserve_null)
    set_id = set[0]
    datatype = set[2]

    query = 'SELECT * FROM configs WHERE sets_id = ' + set_id + ';'
    red_mysql_cmd,handle,query,conf_ans,nl,debug=debug
    if nl eq 1 then begin
      print, inam, ': There is no entry in configs table for ' + set_id + ' sets_id.\r'
      print,'Check the database integrity.'
      return
    endif
    Nconfs = nl-1                ; number of configs in the dataset

    for iconf=1,Nconfs do begin
      conf = strsplit(conf_ans[iconf],tab,/extract,/preserve_null)
      config_id = conf[0]

      camera = conf[2]          ; need this exact variable name to generate dirname
      ;; skip database query if we don't need an information
      if use_strings then $
        if where(strmatch(ucams,'*'+camera+'*')) eq -1 then continue

      state.gain = float(conf[3])
      state.exposure = float(conf[4])
      state.nframes = fix(conf[7]) ; NAXIS3 
      state.camera = camera
      state.is_wb = strmatch(camera,'*-[DW]')
      state.cam_settings = strtrim(string(state.exposure*1000, format = '(f9.2)'), 2) + 'ms'$
              + '_G' + string(state.gain, format = '(f05.2)')                                      

      det_id = conf[11]
      ; get detector information
      query = 'SELECT manufacturer, model, serial_number, name, detfirm FROM detectors WHERE id = ' + det_id + ';'
      red_mysql_cmd,handle,query,det_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, ': There is no entry in detectors table for detector_id = ' + det_id
        print,'Check the database integrity.'
        return
      endif
      dets = strsplit(det_ans[1],tab,/extract,/preserve_null)
      detector = dets[3]        ; need this exact variable name to generate filename
      state.detector = detector

      ;; get dirname generating string
      tmplt_id = conf[12]
      query = 'SELECT template FROM dirname_templates WHERE id = ' + tmplt_id + ';'
      red_mysql_cmd,handle,query,tmplt_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, ': There is no entry in dirname_templates table for fnm_templt_id = ' + tmplt_id
        print,'Check the database integrity.'
        return
      endif
      dir_gen = tmplt_ans[1]

      ;; get filename generating string
      tmplt_id = conf[13]
      query = 'SELECT template FROM filename_templates WHERE id = ' + tmplt_id + ';'
      red_mysql_cmd,handle,query,tmplt_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, 'There is no entry in filename_templates table for fnm_templt_id = ' + tmplt_id
        print,'Check the database integrity.'
        return
      endif
      fnm_gen = tmplt_ans[1]

      ;; get infromation from burst table                       
      query = 'SELECT * FROM bursts WHERE config_id = ' + config_id + ';' 
      red_mysql_cmd,handle,query,burst_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, 'There is no entry in bursts table for the configuration ' + date_time[iset] + ' ' + conf[2]
        print,'Check the database integrity.'
        return
      endif

      Nbursts = nl-1
      brst_msg = string('Get DB values for ', strtrim(string(Nbursts),2), ' bursts.')
      st = replicate( {CHROMIS_STATE}, Nbursts )
      
      for ifile = 0, Nbursts-1 do begin
        red_progressbar, ifile, Nbursts, /predict, brst_msg

        burst = strsplit(burst_ans[ifile+1],tab,/extract,/preserve_null)         
        scannum = burst[3] ; need this for filename generation        
        first_frame = burst[6]  ; need this for filename generation
        state.scannumber = long(scannum)
        state.framenumber = long(first_frame)
        cal_id = burst[8]
        filter1 = fix(burst[9]) ; need this for filename generation
        mosaic_pos = burst[10]  ; need this for filename generation
        wheel = fix(burst[7])   ; need this for filename generation
        hrz = long(burst[2])    ; need this for filename generation

        state.fullstate = state.cam_settings 
        
        if datatype ne 'darks' then begin 
          if year lt 2023 then begin
            ;; get calibration data for CHROMIS to convert wheel*_hrz* into line+tuning
            query = 'SELECT filter1, convfac, du_ref, lambda_ref FROM calibrations WHERE id = ' + cal_id + ';'
            red_mysql_cmd,handle,query,calib_ans,nl,debug=debug
            if nl eq 1 then begin
               print,inam, ': There is no entry in calibrations table for id ' + cal_id
               print,'Check the database integrity.'
               return
            endif
            calib = strsplit(calib_ans[1],tab,/extract,/preserve_null)
            nbpref = calib[0]
            ;; convert 'wheel_hrz' to 'line_+tuning'
            convfac = double(calib[1])
            du_ref = double(calib[2])
            lambda_ref = double(calib[3]) 
            dlambda = convfac * (hrz-du_ref)
            lambda_ref_string = string(round(lambda_ref*1d10), format = '(I04)')
            tuning_string = strtrim(string(round(dlambda*1d13), format = '(I+)'),2)
            tuning = lambda_ref_string + '_' + tuning_string
            state.tun_wavelength = lambda_ref + dlambda
          endif else begin
            nbpref = strtrim(wheel, 2)
            tuning = strtrim(wheel, 2) + '_' + strtrim(string(hrz, format = '(I+)'),2)
            state.tun_wavelength = double(nbpref)*1d-10 + double(tuning)*1d-10
          endelse        
          
          if ~state.is_wb then begin            
            state.prefilter = nbpref
            ;; get information about prefilter 
            query = 'SELECT * FROM filters WHERE filter1 = ' + nbpref + ';'
            red_mysql_cmd,handle,query,filt_ans,nl,debug=debug
            if nl eq 1 then begin
              print,inam, ': There is no entry in filters table for prefilter = ' +  burst[9]
              print, 'Check the database integrity.'
              return
            endif
            filt = strsplit(filt_ans[1],tab,/extract,/preserve_null)
            waveunit = fix(filt[5])
            state.pf_wavelength = float(filt[4])*10.^waveunit

            state.tuning = tuning
            state.fpi_state = state.tuning ;nbpref + '_' + state.tuning                        
            state.fullstate += '_' + nbpref + '_' + tuning
          endif else begin
            state.prefilter = burst[9]
            state.tuning = burst[9] + '_+0'
            if datatype eq 'flats' then $        ; WB flats
               state.fpi_state = burst[9] + '_+0' $ ;fake state
            else $
               state.fpi_state = tuning           
            state.pf_wavelength = double(filter1)*1d-10
            state.tun_wavelength = state.pf_wavelength
            state.fullstate += '_' + burst[9] + '_' + burst[9] + '_+0'
          endelse
        endif else begin
          state.tun_wavelength = 0
          state.tuning = ''
          state.fpi_state = ''
        endelse
        
        ; do we need this?
        ;;if states[ifile].tuning eq '0000_+0' then states[ifile].tuning = ''

        ;; generate dirname (required Y...s, datatype, camera)
        if datatype eq 'science' then datatype='data'
        v=execute(dir_gen)
        if ~v then begin
          print,inam,': Failed to generate dirname'
          return
        endif        
        
        ;; generate filename (required detector, scannum, first_frame, wheel, hrz variables)
        v=execute(fnm_gen)
        if ~v then begin
          print,inam,': Failed to generate filename'
          return
        endif
        zz=strpos(fnm,'hrz')
        if strmid(fnm,zz+3,3) eq '***' then stop

        state.filename = dir_root + dir + fnm
        st[ifile] = state
      endfor                    ; bursts

      ;; append states for different cams (configs)
      red_append, states, st 
      
    endfor ; cams (configs)
  endfor   ;states

  if use_strings then begin
    if strmatch(strings[0],'*CHROMIS/data/*') then begin
      ;;We have to generate filenames as link_data routine does.
      Nst = n_elements(states)
      files=strarr(Nst)       
      for ifile=0,Nst-1 do begin
        if year lt 2023 or strmatch(strings[0], '*nostate*') then begin
          if states[ifile].is_wb then begin
            files[ifile] =  states[ifile].camera + '/' +  states[ifile].detector $
                             + '_' + string(states[ifile].scannumber, format = '(i05)') $
                             + '_' + strtrim(states[ifile].prefilter, 2) $
                             + '_' + string(states[ifile].framenumber, format = '(i07)') $
                             + '.fits'      
          endif else begin 
            files[ifile] = states[ifile].camera + '/' + states[ifile].detector $
                            + '_' + string(states[ifile].scannumber, format = '(i05)') $
                            + '_' + strtrim(states[ifile].fullstate, 2) $
                            + '_' + string(states[ifile].framenumber, format = '(i07)')  $
                            + '.fits'
          endelse
        endif else $
            files[ifile] = states[ifile].camera + '/' + file_basename(states[ifile].filename)
      endfor
      ;; Filter the found states with respect to the file names in
      ;; strings
      Ncams=n_elements(ucams)
      for icam=0,Ncams-1 do begin
        ss = where(strmatch(strings,'*'+ucams[icam]+'*'), count)
        if count eq 0 then continue
        fs = ucams[icam] + '/' + file_basename(strings[ss])
        match2, fs, files, suba ;, subb 
        stt = states[suba]
        stt.filename = strings[ss]
        red_append, sts,stt
      endfor
      states = sts
    endif else begin      
      ;; Filter the found states with respect to the file names in strings
      match2, strings, states.filename, suba ;, subb
      if min(suba) eq -1 then begin
        print,'There are no matches in the database for provided filenames.'
        print,'Something is wrong.'
        stop
      endif
      states = states[suba]
    endelse
    ;; Store in cache
    Nstr = n_elements(strings)
    for ifile=0,Nstr-1 do $
      rdx_cache, strings[ifile], { state:states[ifile] }
  endif 

  free_lun,handle  
end
