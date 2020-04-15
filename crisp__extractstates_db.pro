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
;      A list of datasets ('YYYY-MM-DD HH:MM:SS) for which to extract the states information.
;   
;    states : out, type=array(struct)
;
;        An array of structs, containing (partially filled) state information.
; 
; :History:
; 
;   2019-07-10 : OA. Created. (Derived from crisp__extractstates and red_readhead_db)
;   
;   2019-07-11 : MGL. Call like extractstates.
;
;   2020-03-11 : OA. Added check for unique cameras when in the procedure call
;                filenames(strings) are used instead of datasets. Added
;                keyword 'cam'.
;
;-
pro crisp::extractstates_db, strings, states, datasets = datasets, cam = cam

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  use_strings = n_elements(datasets) eq 0

  ;; cache is slow comparing to the database at least when BeeGFS is fatigue
  ;; let's always set /force till further examination
  ;force = 1B
  
  ;;Clear states and try to get them from cache.
  undefine,states
  if use_strings then begin
    Nstr = n_elements(strings)
    ;; if ~keyword_set(force) then begin
    ;;   for ifile=0,Nstr-1 do begin
    ;;     this_cache = rdx_cacheget(strings[ifile], count = cnt)   
    ;;     if cnt gt 0 then begin
    ;;       red_progressbar, ifile, Nstr, 'Extract state info from cache', /predict
    ;;       red_append,states,this_cache.state
    ;;     endif
    ;;   endfor
    ;;   if n_elements(states) eq Nstr then return else undefine,states
    ;; endif
    timestamps = stregex(strings,'[0-2][0-9]:[0-5][0-9]:[0-5][0-9]', /extract)
    utimestamps = timestamps[uniq(timestamps,sort(timestamps))]
    datasets = self.isodate + ' ' + utimestamps
    
    if keyword_set(cam) then ucams = [cam] else begin
      cms = strarr(Nstr)
      for ii = 0, Nstr-1 do begin
        dr = file_dirname(strings[ii])
        cms[ii] = (strsplit(dr,'/',/extract))[-1]
      endfor
      ucams = cms[uniq(cms,sort(cms))]
    endelse
  endif
  
  instrument = 'CRISP'

  red_mysql_check, handle

  Nsets = n_elements(datasets)
  set_msg = string('Get metadata from the database for ', strtrim(string(Nsets),2), ' datasets.')

  split_dir = strsplit(self.root_dir,'/',/extract)
  dir_root = '/' + strjoin(split_dir[0:-2],'/') + '/'

  for iset = 0, Nsets-1 do begin

    red_progressbar, iset, Nsets, /predict, set_msg     
    split_set = strsplit(datasets[iset], ' ', /extract)
    timestamp = split_set[1]
    date_time = datasets[iset]
    
    ; we need Y,...,sec to generate the directory name
    split_date = strsplit(split_set[0],'-',/extract)
    Y = fix(split_date[0])
    M = fix(split_date[1])
    D = fix(split_date[2])

    split_time = strsplit(timestamp,':',/extract)
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
      pinh_dirs = *self.pinh_dirs
      polcal_dirs = *self.polcal_dir     
      data_dirs = *self.data_dirs
      dirs = [flat_dirs, dark_dirs, pinh_dirs, polcal_dirs, data_dirs]
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
    ;parse a result of the query
    tab = string(9B)
    set = strsplit(ans[1],tab,/extract,/preserve_null)
    set_id = set[0]
    datatype = set[2]
    case datatype of
      'darks': datatype='Darks'
      'flats': datatype='Flats'
      'pinholes': datatype='Pinholes'
      'polcal': datatype='Polcal'
      'science': datatype='Science'
    endcase
    if datatype eq 'Polcal' then begin
      state = {CRISP_POLCAL_STATE}
    endif else begin
      state = {CRISP_STATE}
    endelse

    query = 'SELECT * FROM configs WHERE sets_id = ' + set_id + ';'
    red_mysql_cmd,handle,query,conf_ans,nl,debug=debug
    if nl eq 1 then begin
      print, inam, ': There is no entry in configs table for ' + instrument + ' ' + date_time[iset] + ' dataset.\r'
      print,'Check the database integrity.'
      return
    endif
    Ncams = nl-1                ; number of cameras (configs) in the dataset
    ;cams=strarr(Ncams)

    for icam=1,Ncams do begin
      conf = strsplit(conf_ans[icam],tab,/extract,/preserve_null)
      config_id = conf[0]

      camera = conf[2]          ; need this exact variable name to generate dirname
      ;; skip database query if we don't need an information
      if use_strings then $
        if where(strmatch(ucams,'*'+camera+'*') eq 1) eq -1 then continue
      
      state.exposure = float(conf[4])
      state.camera = camera
      ;cams[icam-1]=camera
      state.is_wb = strmatch(camera,'*-[DW]')
      state.cam_settings = strtrim(string(state.exposure*1000, format = '(f9.2)'), 2) + 'ms'      

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

      ; get dirname generating string
      tmplt_id = conf[12]
      query = 'SELECT template FROM dirname_templates WHERE id = ' + tmplt_id + ';'
      red_mysql_cmd,handle,query,tmplt_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, ': There is no entry in dirname_templates table for fnm_templt_id = ' + tmplt_id
        print,'Check the database integrity.'
        return
      endif
      dir_gen = tmplt_ans[1]

      ; get filename generating string
      tmplt_id = conf[13]
      query = 'SELECT template FROM filename_templates WHERE id = ' + tmplt_id + ';'
      red_mysql_cmd,handle,query,tmplt_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, 'There is no entry in filename_templates table for fnm_templt_id = ' + tmplt_id
        print,'Check the database integrity.'
        return
      endif
      fnm_gen = tmplt_ans[1]

      ; get infromation from burst table                       
      query = 'SELECT * FROM bursts WHERE config_id = ' + config_id + ';' 
      red_mysql_cmd,handle,query,burst_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, 'There is no entry in bursts table for the configuration ' + date_time[iset] + ' ' + conf[2]
        print,'Check the database integrity.'
        return
      endif

      Nbursts = nl-1
      brst_msg = string('Get metadata from the database for ', strtrim(string(Nbursts),2), ' bursts.') 

      state.nframes = 1         ; Single frame in CRISP files      

      for ifile = 0, Nbursts-1 do begin
        burst = strsplit(burst_ans[ifile+1],tab,/extract,/preserve_null)         
        scannum = burst[3] ; need this exact variable name to generate filename
        first_frame = burst[6]
        state.scannumber = scannum

        red_progressbar, ifile, Nbursts, /predict, brst_msg

        line = long(burst[7])      ; required for filename generation
        tuning = long(burst[2])   ; required for filename generation
        filter1 = fix(burst[9])   ; required for filename generation
        burst_id = burst[0]        
        
        state.prefilter = burst[9]

         ;; get information about prefilter if data are not darks
        if burst[9] ne '0' then begin
          query = 'SELECT * FROM filters WHERE filter1 = ' + burst[9] + ';'
          red_mysql_cmd,handle,query,filt_ans,nl,debug=debug
          skip_pref = 0B
          if nl eq 1 then begin
            if state.is_wb then begin
              print,inam, ', Warning: There is no entry in filters table for prefilter = ' +  burst[9]
              print, 'Perhaps WB flats were taken with more prefilters than NB flats and science data.'
              skip_pref = 1B
            endif else begin
              print,inam, ': There is no entry in filters table for prefilter = ' +  burst[9]
              print, 'Check the database integrity.'
              return
            endelse
          endif
        endif else skip_pref = 1B
        if ~skip_pref then begin
          filt = strsplit(filt_ans[1],tab,/extract,/preserve_null)
          waveunit = fix(filt[5])
          state.pf_wavelength = float(filt[4])*10.^waveunit ; nm
        endif else  state.pf_wavelength = 0.

        if state.is_wb then begin
          state.fpi_state =  burst[7] + '_' + strtrim(string(tuning, format='(I+)'),2)
          state.tuning = burst[7] + '_+0'
          state.tun_wavelength = state.pf_wavelength
        endif else begin
          state.tuning = burst[7] + '_' + strtrim(string(tuning, format='(I+)'),2)
          state.fpi_state = state.tuning
          state.tun_wavelength = float(line + tuning/1000.)*1e-10 ; Å
        endelse

        if datatype eq 'Darks' then begin
          state.tun_wavelength = 0
          state.tuning = '0'
          state.fpi_state = '0'
          state.fullstate = '0_lc0'
          state.prefilter=''
        endif

        if datatype eq 'Polcal' then table='crisp_polcal_frames' $
          else table = 'crisp_frames'
        query = 'SELECT * FROM ' + table + ' WHERE burst_id = ' + burst_id + ';'
        red_mysql_cmd,handle,query,frames_ans,nl,debug=debug
        if nl eq 1 then begin
          print, inam, ': There is no entry in ' + table + ' table for the burst_id ' + burst_id
          print, 'Check the database integrity.'
          return
        endif
        Nframes = nl-1

        if datatype eq 'Polcal' then begin ;keyword_set(polcal)
          st = replicate( {CRISP_POLCAL_STATE}, Nframes )
        endif else begin
          st = replicate( {CRISP_STATE}, Nframes )
        endelse

          ; generate dirname (required Y...s, datatype, camera)    
        v=execute(dir_gen)
        if ~v then begin
          print,inam,': Failed to generate dirname'
          print,'Check the database integrity.'
          return
        endif 

        for iframe=0, Nframes-1 do begin              
          frame = strsplit(frames_ans[iframe+1],tab,/extract,/preserve_null)
          framenum = long(frame[3])
          lc_state = fix(frame[8])
          state.framenumber = long(framenum)
          if state.is_wb then begin
            state.fullstate = burst[9] + '_' + burst[7] + '_+0'
          endif else begin
            state.fullstate = burst[9] + '_' + burst[7] + '_' + strtrim(string(tuning, format='(I+)'),2)
          endelse
          if lc_state ne -1 then state.lc = lc_state ; 
          if ~state.is_wb and datatype ne 'Darks' then state.fullstate += '_lc' + frame[8]
          if datatype eq 'Polcal' then begin
            qw_state = fix(frame[9])
            lp_state = fix(frame[10])
            state.fullstate = 'lp' + string(lp_state, format = '(i03)') + '_' $
                              + 'qw' + string(qw_state, format = '(i03)') + '_' $
                              + state.fullstate
            state.qw = qw_state
            state.lp = lp_state
          endif          
          
             ; generate filename              
          v=execute(fnm_gen)
          if ~v then begin
            print, inam, ': Failed to generate filename \r'
            print, 'Check the database integrity.'
          endif
          state.filename = dir_root + dir + fnm

          st[iframe] = state
        endfor     ; iframe
                               
        red_append, states, st        
      endfor                    ; ifile (burst)
    endfor  ; cams (configs)
  endfor  ;datasets

  if use_strings then begin
    if strmatch(strings[0],'*CRISP/data*') then begin
      ;;We have to generate filenames as link_data routine does.
      Nst = n_elements(states)
      files=strarr(Nst)
      if strmatch(strings[0],'*nostate*') then begin        
        for ifile=0,Nst-1 do $
          files[ifile] =  states[ifile].camera + '/' + states[ifile].detector $
                   + '_' + string(states[ifile].scannumber, format = '(i05)') $
                   + '_' + strtrim(states[ifile].prefilter, 2) $
                   + '_' + string(states[ifile].framenumber, format = '(i07)')      
      endif else begin
        for ifile=0,Nst-1 do $
          files[ifile] = states[ifile].camera + '/' + states[ifile].detector $
                 + '_' + string(states[ifile].scannumber, format = '(i05)') $
                 + '_' + strtrim(states[ifile].fullstate, 2) $
                 + '_' + string(states[ifile].framenumber, format = '(i07)') 
      endelse
      ;; Filter the found states with respect to the file names in
      ;; strings
      Ncams = n_elements(ucams)
      for icam=0,Ncams-1 do begin
        ss = where(strmatch(strings,'*'+ucams[icam]+'*'))
        if n_elements(ss) eq 1 then if ss eq -1 then continue
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
      states = states[suba]
    endelse
    ;; Store in cache
    Nstr = n_elements(strings)
    for ifile=0,Nstr-1 do $
      rdx_cache, strings[ifile], { state:states[ifile] }
  endif
  
  free_lun,handle  
end

