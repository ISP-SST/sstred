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
;   2019-07-10 : OA. Created. (Derived from chromis__extractstates and red_readhead_db)
;
;-
pro chromis::extractstates_db, datasets, states 
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  instrument = 'CHROMIS'

  red_mysql_check, handle

  Nsets = n_elements(datasets)
  set_msg = string('Get DB values for ', strtrim(string(Nsets),2), ' datasets.')

  split_dir = strsplit(self.root_dir,'/',/extract)
  dir_root = '/' + strjoin(split_dir[0:-2],'/') + '/'

  for iset = 0, Nsets-1 do begin
    state = {CHROMIS_STATE}
    red_progressbar, iset, Nsets, /predict, set_msg     
    split_set = strsplit(datasets[iset], ' ', /extract)
    date_time = datasets[iset];split_set[0] + ' ' + split_set[1]
    
    ; we need Y,...,sec to generate the directory name
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
      dark_dirs = *self_dark_dir
      data_dirs = *self.data_dirs
      dirs = [flat_dirs, dark_dirs, data_dirs]
      in = where(strmatch(dirs,'*'+split_set[1]+'*'))
      dir = dirs[in]
      red_rawdir2db,dir=dir,/all
      red_mysql_cmd,handle,query,ans,nl
    endif
    ;parse a result of the query
    tab = string(9B)
    set = strsplit(ans[1],tab,/extract,/preserve_null)
    set_id = set[0]
    datatype = set[2]

    query = 'SELECT * FROM configs WHERE sets_id = ' + set_id + ';'
    red_mysql_cmd,handle,query,conf_ans,nl,debug=debug
    if nl eq 1 then begin
      print, inam, ': There is no entry in configs table for ' + instrument + ' ' + date_time[iset] + ' ' + cameras[icam] + ' dataset.\r'
      print,'Check the database integrity.'
      return
    endif
    Ncams = nl-1 ; number of cameras (configs) in the dataset

    for icam=1,Ncams do begin
      conf = strsplit(conf_ans[icam],tab,/extract,/preserve_null)
      config_id = conf[0]

      state.gain = float(conf[3])
      state.exposure = float(conf[4])
      state.nframes = fix(conf[7]) ; NAXIS3 
      camera = conf[2] ; need this exact variable name to generate dirname
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

        if cal_id ne 0 then begin ;  otherwise it can be darks ...               
               ; get calibration data for CHROMIS to convert wheel*_hrz* into line+tuning
          query = 'SELECT prefilter, convfac, du_ref, lambda_ref FROM calibrations WHERE id = ' + cal_id + ';'
          red_mysql_cmd,handle,query,calib_ans,nl,debug=debug
          if nl eq 1 then begin
            print,inam, ': There is no entry in calibrations table for id ' + cal_id
            print,'Check the database integrity.'
            return
          endif

          calib = strsplit(calib_ans[1],tab,/extract,/preserve_null)
          nbpref = calib[0]
          convfac = float(calib[1])
          du_ref = float(calib[2])
          lambda_ref = float(calib[3])
          wheel = fix(burst[7])  ; required for filename generation
          hrz = long(burst[2])   ; required for filename generation

          ; get information about prefilter
            query = 'SELECT * FROM filters WHERE prefilter = ' + nbpref + ';'
          red_mysql_cmd,handle,query,filt_ans,nl,debug=debug
          if nl eq 1 then begin
            print,inam,': There is no entry in filters table for prefilter = ' + nbpref
            print,'Check the database integrity.'
            return
          endif
          filt = strsplit(filt_ans[1],tab,/extract,/preserve_null)
          state.pf_wavelength = float(filt[4])
          
          if set[2] eq 'flats' and strmatch(camera, '*-[WD]') then begin
                ; WB flats 
            state.fpi_state = nbpref + '_' + nbpref + '_+000' ;fake state
            state.tun_wavelength = states[ifile].pf_wavelength
            states.tuning = string(round(states[ifile].tun_wavelength*1d10) $
                   , format = '(i04)')  + '_+0'
          endif else begin
                 ;convert 'wheel_hrz' to 'line_+tuning'
            dlambda = convfac * (hrz-du_ref)
            lambda_ref_string = string(round(lambda_ref*1d10), format = '(i04)')
            tuning_string = strtrim(round(dlambda*1d13), 2)
            if strmid(tuning_string, 0, 1) ne '-' then tuning_string = '+'+tuning_string 
            state.fpi_state = nbpref + '_' + lambda_ref_string + '_' + tuning_string
            state.tuning = lambda_ref_string + '_' + tuning_string
            ;; Also return the tuning in decimal form [m]:
            state.tun_wavelength = lambda_ref + dlambda
          endelse
        endif

        if set[2] eq 'darks' then begin
          state.tun_wavelength = 0
          state.tuning = ''
          state.fpi_state = ''
        endif
        
        ; do we need this?
        ;;if states[ifile].tuning eq '0000_+0' then states[ifile].tuning = ''

          ; generate dirname (required Y...s, datatype, camera)
        if datatype eq 'science' then datatype='data'
        v=execute(dir_gen)
        if ~v then begin
          print,inam,': Failed to generate dirname'
          return
        endif        
        
        ; generate filename (required detector, scannum, first_frame, wheel, hrz variables)
        v=execute(fnm_gen)
        if ~v then begin
          print,inam,': Failed to generate filename'
          return
        endif
        zz=strpos(fnm,'hrz')
        if strmid(fnm,zz+3,3) eq '***' then stop

        state.filename = dir_root + dir + fnm        

        ;; The fullstate string
        undefine, fullstate_list
        if ~keyword_set(strip_settings) then red_append, fullstate_list, state.cam_settings
        if state.prefilter ne '' then red_append, fullstate_list, state.prefilter
        if state.tuning ne '' then begin     
          if keyword_set(strip_wb) then begin
            if state.is_wb eq 0 then $
               red_append, fullstate_list, state.tuning
          endif else begin
            red_append, fullstate_list, state.tuning
          endelse
        endif
        state.fullstate = strjoin(fullstate_list, '_')

        st[ifile] = state

      ;; Store in cache
      ;rdx_cache, strings[ifile], { state:states[ifile] }

      endfor                    ; bursts

      ; append states for different cams (configs)
      red_append, states, st 
      
    endfor ; cams (configs)
  endfor  ;states
  
end
