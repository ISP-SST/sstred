; docformat = 'rst'

;+
; This function constructs fits headers with help of information taken
; from the database. It might substitute red_readhead function.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Oleksii Andriienko, Institute for Solar Physics
; 
; 
; :Returns:
;
;   A list of FITS headers or a single FITS header (if a single
;   filename is provided as input).
;
;
; :Params:
; 
;   files: in, type=strarr
; 
;     An array of filenames.
; 
; 
; :Keywords:
; 
;    date_beg: out, optional, type=list
;
;      A list of string arrays (or a single array if a single
;      filename is provided as input) with timestamps for exposure start.
;
;    framenumbers: out, optional, type=list
;      
;      A list of integer arrays (or a single array if a single
;      filename is provided as input) with framenumbers.
;
;    status: out, optional, type=list
;
;      A list of structures with return status (0 for success, -1 for
;      failure) and error messages. If a single filename is provided
;      status is a signed integer.
;   
; 
; 
; :History:
; 
; 2019-05-28 : OA. First version.
;
; 2019-06-06 : OA. Rewritten. Loop for cameras added.
; 
;-

function red_readhead_db, files, date_beg=date_beg, framenumbers=framenumbers, status=status

  inam = red_subprogram(/low, calling = inam1)

  Nfls = n_elements(files)
  if Nfls eq 0 then begin
    print, inam, ': There is no filenames in the list.'
    status = -1
    return, 0B
  endif

  headers = list()
  status = list()
  date_beg = list()
  framenumbers = list()
  st = {code: 0, msg: ''}

  instr = ['CHROMIS', 'Crisp']
  for ins=0,1 do begin
    instr_indx = where(strpos(files, instr[ins]) ne -1)
    ; We expect that among submitted filenames there are ones with 'Crisp'
    ; or 'CHROMIS' substrings.
    if n_elements(instr_indx) eq 1 then if instr_indx eq -1 then continue ; no files with the instrument, switch to other one
    
    red_mysql_check, handle

    files_instr = files(instr_indx)
    sets = stregex(files_instr, '(20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]).*([0-2][0-9]:[0-5][0-9]:[0-6][0-9])',/extract)
    set_indx = uniq(sets,sort(sets))
    sets = sets(set_indx)   
    Nsets = n_elements(set_indx)
    
    red_mkhdr,hdr, '' ; create a dummy header
    red_fitsaddkeyword, hdr, 'ORIGIN', 'Institute for Solar Physics'
    red_fitsaddkeyword, hdr, 'TELESCOP', 'Swedish 1-meter Solar Telescope'
    red_fitsaddkeyword, hdr, 'OBJECT', 'Sun'
    red_fitsaddkeyword, hdr, 'DATAMIN', 0
    red_fitsaddkeyword, hdr, 'DATAMAX', 4095
    

    set_msg = string('Get DB values for ', strtrim(string(Nsets),2), ' datasets.')
    
    for iset = 0, Nsets-1 do begin
      red_progressbar, iset, Nsets, /predict, set_msg
      files_set = files_instr(where(stregex(files_instr,sets(iset)) ne -1))
      instrument = strupcase(instr[ins])
      date = stregex(sets(iset), '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]',/extract)
      time = stregex(sets(iset),'[0-2][0-9]:[0-5][0-9]:[0-6][0-9]',/extract)
      date_time = date + ' ' + time
      query='SELECT * FROM datasets WHERE date_obs = "' + date_time + '" AND instrument = "' + instrument + '";'
      red_mysql_cmd,handle,query,ans,nl,debug=debub
      if nl eq 1 then begin
        st.msg += 'There is no entries in the database for ' + date_time + ' ' + instrument + ' dataset. \r'
        st.code = -1
        status.add, st
        framenumbers.add, [0]
        date_beg.add, ['0']
        headers.add,hdr
        continue
      endif

      ;parse a result of the query
      tab = string(9B)
      set = strsplit(ans[1],tab,/extract,/preserve_null)
      set_id = set[0]
      dt = strsplit(set[1],' ',/extract)
      red_fitsaddkeyword, hdr, 'DATE-OBS', dt[0] + 'T' + dt[1]
;      red_fitsaddkeyword, hdr, 'OBSERVER',  set[2]
;      red_fitsaddkeyword, hdr, 'DESCRIPTION', reply[3]
      red_fitsaddkeyword, hdr,'INSTRUME', instrument, 'Name of instrument'
      red_fitsaddkeyword, hdr,'DATATYPE', set[2]
      red_fitsaddkeyword, hdr,'SOLARNET', float(set[4]), 'Fully SOLARNET-compliant=1.0, partially=0.5'     

      query = 'SELECT * FROM configs WHERE sets_id = ' + set_id + ';'  ;' AND camera = "' + camera + '" ;'
      red_mysql_cmd,handle,query,conf_ans,nl,debug=debug
      if nl eq 1 then begin
        print, inam, ': There is no entry for set_id', set_id, ' dataset in the configs table.'
        print,'Check the database integrity.'
        return, 0B
      endif
            
      Ncams = nl-1
      for icam = 1,Ncams do begin
        conf = strsplit(conf_ans[icam],tab,/extract,/preserve_null)
        config_id = conf[0]
        camera = conf[2]
        inn = where(strpos(files_set,camera) ne -1)
        if n_elements(inn) eq 1 then if inn eq -1 then continue
        files_cam = files_set[inn]
      
        red_fitsaddkeyword, hdr, 'CAMERA', conf[2]      
        red_fitsaddkeyword, hdr, 'XPOSURE', float(conf[4]), '  [s]'
        red_fitsaddkeyword, hdr, 'NAXIS1', fix(conf[5]), ' naxis1 of original image'
        red_fitsaddkeyword, hdr, 'NAXIS2', fix(conf[6]), ' naxis2 of original image'
        red_fitsaddkeyword, hdr, 'BITPIX', fix(conf[8]), ' bitpix of original image'
        red_fitsaddkeyword, hdr, 'CADAVG', float(conf[9]), ' [s]'
        red_fitsaddkeyword, hdr, 'DETOFFS', fix(conf[10]), ' [counts] or camera specific unit'
        naxs3 = fix(conf[7])
        if naxs3 eq 0 then $
           red_fitsaddkeyword, hdr, 'NAXIS', 2 $
        else begin 
          red_fitsaddkeyword, hdr, 'NAXIS', 3
          red_fitsaddkeyword, hdr, 'NAXIS3', naxs3, ' naxis3 of original image'
        endelse          
      
        det_id = conf[11]
          ; get detector information
        query = 'SELECT manufacturer, model, serial_number, name, detfirm FROM detectors WHERE id = ' + det_id + ';'
        red_mysql_cmd,handle,query,det_ans,nl,debug=debug
        if nl eq 1 then begin
          print, inam, ':There is no entry in detectors table for detector_id = ' + det_id
          print,'Check the database integrity.'
          return, 0B
        endif
        detectors = strsplit(det_ans[1],tab,/extract,/preserve_null)
        detector = detectors[3] ; need this exact variable name to generate filename
        red_fitsaddkeyword, hdr, 'DETECTOR', detector, 'Detector identifier'
        red_fitsaddkeyword, hdr, 'DETFIRM', detectors[4]
        red_fitsaddkeyword, hdr, 'DETMODEL', detectors[0] + ' ' + detectors[1] + ' serial ' + detectors[2]
        
           ; get filename generating string
        tmplt_id = conf[13]
        query = 'SELECT template FROM filename_templates WHERE id = ' + tmplt_id + ';'
        red_mysql_cmd,handle,query,tmplt_ans,nl,debug=debug
        if nl eq 1 then begin
          print, inam, ':There is no entry in filename_templates table for fnm_templt_id = ' + tmplt_id
          print,'Check the database integrity.'
          return, 0B
        endif
        fnm_gen = tmplt_ans[1]
        
            ; get infromation from burst table                       
        query = 'SELECT * FROM bursts WHERE config_id = ' + config_id + ';' 
        red_mysql_cmd,handle,query,burst_ans,nl,debug=debug
        if nl eq 1 then begin
          print, inam, ':There is no entry in bursts table for the configuration ' + date_time[iset] + ' ' + conf[2]
          print,'Check the database integrity.'
          return, 0B
        endif
        Nbursts = nl-1
        brst_msg = string('Get DB values for ', strtrim(string(Nbursts),2), ' bursts.')
      
        for ifile = 0, Nbursts-1 do begin
          burst = strsplit(burst_ans[ifile+1],tab,/extract,/preserve_null)         
          scannum = burst[3]
          first_frame = burst[6]                    

          red_progressbar, ifile, Nbursts, /predict, brst_msg
          red_fitsaddkeyword, hdr, 'SCANNUM',  long(scannum)          

          if ins eq 0 then begin ; CHROMIS
            egex = string('(',scannum,').*(',first_frame,')', format='(A,I05,A,I07,A)')
            inn = where(stregex(files_cam,egex) ne -1)
            if n_elements(inn) eq 1 then if inn eq -1 then continue
            if n_elements(inn) gt 1 then begin
              print, inam, ': Repeated filenames in the list?.'
              stop
            endif
            file_burst = files_cam[inn] ; must be one file
          
            red_fitsaddkeyword, hdr, 'EXTEND','T'
            red_fitsaddkeyword, hdr, 'OBS-HDU', 1
            red_fitsaddkeyword, hdr, 'TAB_HDUS', 'TABULATIONS;DATE-BEG'
            dt = strsplit(burst[1],' ',/extract)
            red_fitsaddkeyword, hdr, 'DATE', dt[0] + 'T' + dt[1]
            red_fitsaddkeyword, hdr, 'DETGAIN', float(conf[3]), '  [dB] or camera specific unit'
            red_fitsaddkeyword, hdr, 'DETTEMP', float(burst[4]), '/ [C] Current temperature of the detector'
            cal_id = burst[8]

            if cal_id ne 0 then begin ;  otherwise it can be darks ...               
               ; get calibration data for CHROMIS to convert wheel*_hrz* into line+tuning
              query = 'SELECT prefilter, convfac, du_ref, lambda_ref FROM calibrations WHERE id = ' + cal_id + ';'
              red_mysql_cmd,handle,query,calib_ans,nl,debug=debug
              if nl eq 1 then begin
                print, inam, ': There is no entry in calibrations table for id ' + cal_id 
                print,'Check the database integrity.'
                return,0B
              endif

              calib = strsplit(calib_ans[1],tab,/extract,/preserve_null)
              nbpref = calib[0]
              convfac = float(calib[1])
              du_ref = float(calib[2])
              lambda_ref = float(calib[3])
              wheel = fix(burst[7]) ; required for filename generation
              hrz = long(burst[2])  ; required for filename generation
              
              if set[2] eq 'flats' and strmatch(camera, '*-[WD]') then begin
                                ; WB flats 
                red_fitsaddkeyword, hdr, 'STATE', nbpref + '_' + nbpref + '_+000' ;fake state
              endif else begin
                   ;convert 'wheel_hrz' to 'line_+tuning'
                dlambda = convfac * (hrz-du_ref)
                lambda_ref_string = string(round(lambda_ref*1d10), format = '(i04)')
                tuning_string = strtrim(round(dlambda*1d13), 2)
                if strmid(tuning_string, 0, 1) ne '-' then tuning_string = '+'+tuning_string
                state = nbpref + '_' + lambda_ref_string + '_' + tuning_string              
                red_fitsaddkeyword, hdr, 'STATE', state
              endelse
              red_fitsaddkeyword, hdr, 'FILTER1', nbpref

                ; get information about prefilter
              query = 'SELECT * FROM filters WHERE prefilter = ' + nbpref + ';'
              red_mysql_cmd,handle,query,filt_ans,nl,debug=debug
              if nl eq 1 then begin
                print, inam, ': There is no entry in filters table for prefilter = ' + nbpref
                print,'Check the database integrity.'
                return,0B
              endif
              filt = strsplit(filt_ans[1],tab,/extract,/preserve_null)
              red_fitsaddkeyword, hdr, 'FILTER1', filt[1]
              red_fitsaddkeyword, hdr, 'WAVEMAX', float(filt[2]), ' [nm] Prefilter max wavelength (0.5 peak)'
              red_fitsaddkeyword, hdr, 'WAVEMIN', float(filt[3]), ' [nm] Prefilter min wavelength (0.5 peak)'
              red_fitsaddkeyword, hdr, 'WAVELNTH', float(filt[4]), ' [nm] Prefilter peak wavelength'
              red_fitsaddkeyword, hdr, 'WAVEUNIT', fix(filt[5]), ' WAVELNTH in units 10^WAVEUNIT m = nm'
              switch filt[1] of
                '3625' :
                '3934' :
                '3950' :
                '3969' :
                '3978' :
                '3999' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','Ca II H & K'
                  break
                end
                '4862' :
                '4846' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','H-beta'
                  break
                end
              endswitch

            endif               ; cal_id ne 0

               ; generate filename (required detector, scannum, first_frame, wheel, hrz variables)
            v=execute(fnm_gen)
            if ~v then begin
              print, inam, ': Failed to generate filename.'
              print,'Check the database integrity.'
              return, 0B
            endif
            red_fitsaddkeyword, hdr, 'FILENAME', fnm

            burst_id = burst[0]
            query = 'SELECT * FROM chromis_frames WHERE burst_id = ' + burst_id + ';'
            red_mysql_cmd,handle,query,frames_ans,nl,debug=debug
            if nl eq 1 then begin
              print, inam, ': There is no entry in chromis_frames table for the burst_id ' + burst_id 
              print,'Check the database integrity.'
              return, 0B
            endif 
            Nframes = nl-1

            framenums = intarr(Nframes)
            date_bgs = strarr(Nframes)
            for iframe=0, Nframes-1 do begin
              frame = strsplit(frames_ans[iframe+1],tab,/extract,/preserve_null)
              framenums[iframe] = frame[3]
              frac = strsplit(frame[2],'.',/extract)
              date_bgs[iframe] = frame[1]+'.'+frac[1]
            endfor          

            red_fitsaddkeyword, hdr, 'DATE_BEG', dt[0] + 'T' + date_bgs[0], ' First in table' ; date + time
            bg = red_Time_conv(date_bgs[0])
            nd = red_Time_conv(date_bgs[Nframes-1])            
            d_end = nd + float(conf[4]) ; last frame beginning + exposure time
            tm = strsplit(string(d_end, format='(f12.6)'),'.',/extract)
            frac = tm[1]
            time = red_Time_conv(long(d_end),/f2s)
            red_fitsaddkeyword, hdr, 'DATE_END', dt[0] + 'T' + time + '.' + frac, ' Last in table + XPOSURE'
            avrg = bg + (d_end - bg)/2.
            tm = strsplit(string(avrg, format='(f12.6)'),'.',/extract)
            frac = tm[1]
            time = red_Time_conv(long(avrg),/f2s)
            red_fitsaddkeyword, hdr, 'DATE_AVG', dt[0] + 'T' + time + '.' + frac, ' Average time from table'

            status.add,st
            framenumbers.add, framenums
            date_beg.add, date_bgs
            headers.add,hdr

          endif else begin      ; CRISP

            egex = string('(',scannum,').*(',first_frame,')', format='(A,I05,A,I07,A)')
            inn = where(stregex(files_cam,egex) ne -1)
            if n_elements(inn) eq 1 then if inn eq -1 then continue
            files_burst = files_cam[inn]
            
            line = fix(burst[7])  ; required for filename generation
            tuning = long(burst[2]) ; required for filename generation
            burst_id = burst[0]            

            if line ne 0 then begin ; otherwise darks or WB flats
              red_fitsaddkeyword, hdr, 'FILTER1', burst[9], 'CRISP prefilter' 
              red_fitsaddkeyword, hdr, 'STATE', burst[9] + '_' + burst[7] + '_' + string(tuning, format='(I+05)'), 'NB tuning state'

                                ; get information about prefilter
              query = 'SELECT * FROM filters WHERE prefilter = ' + burst[9] + ';'
              red_mysql_cmd,handle,query,filt_ans,nl,debug=debug
              if nl eq 1 then begin
                st.msg += 'There is no entry in filters table for prefilter = ' +  burst[9] + '.\r'
                st.code = -1
                status.add, st
                framenumbers.add, [0]
                date_beg.add, ['0']
                headers.add,hdr
                continue
              endif
              filt = strsplit(filt_ans[1],tab,/extract,/preserve_null)
              red_fitsaddkeyword, hdr, 'WAVEMAX', float(filt[2]), ' [nm] Prefilter max wavelength (0.5 peak)'
              red_fitsaddkeyword, hdr, 'WAVEMIN', float(filt[3]), ' [nm] Prefilter min wavelength (0.5 peak)'
              red_fitsaddkeyword, hdr, 'WAVELNTH', float(filt[4]), ' [nm] Prefilter peak wavelength'
              red_fitsaddkeyword, hdr, 'WAVEUNIT', fix(filt[5]), ' WAVELNTH in units 10^WAVEUNIT m = nm'

              switch filt[1] of
                '5173': begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','Mg b ' + burst[7]
                  break
                end
                '5578' :
                '6173' : 
                '5250' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','Fe I ' + burst[7]                     
                  break
                end
                '6302' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','Fe I 6301+6302'
                  break
                end
                '5382' : 
                '5578' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','C I ' + burst[7]
                  break
                end
                '5897' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','Na D ' + burst[7]
                  break
                end
                '6563' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','H-alpha'
                  break
                end
                '7772' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','O I ' + burst[7]
                  break
                end
                '8542' : begin
                  red_fitsaddkeyword, hdr, 'WAVEBAND','Ca II ' + burst[7]
                  break
                end
              endswitch
              
            endif

                                ; for Science WB we have to change state 
            if strmatch(conf[2], '*-[W]') then begin
              red_fitsaddkeyword, hdr, 'FP_STATE', burst[9] + '_' + burst[7] + '_' + string(tuning, format='(I+05)'), 'NB tuning state'
              red_fitsaddkeyword, hdr, 'STATE', burst[9] + '_' + burst[7] + '_+0000'
            endif

            if set[2] eq 'polcal' then table='crisp_polcal_frames' $
            else table = 'crisp_frames'
            query = 'SELECT * FROM ' + table + ' WHERE burst_id = ' + burst_id + ';'
            red_mysql_cmd,handle,query,frames_ans,nl,debug=debug
            if nl eq 1 then begin
              print, inam, ': There is no entry in ' + table + ' table for the burst_id ' + burst_id
              print, 'Check the database integrity.'
              return,0B
            endif
            Nframes = nl-1

            for iframe=0, Nframes-1 do begin              
              frame = strsplit(frames_ans[iframe+1],tab,/extract,/preserve_null)
              framenum = frame[3]
              frac = strsplit(frame[2],'.',/extract)
              date_bgs = frame[1]+'.'+frac[1]
              red_fitsaddkeyword, hdr, 'DATE', dt[0] + 'T' + date_bgs, ' Creation of fz file'
              red_fitsaddkeyword, hdr, 'DATE_BEG', dt[0] + 'T' + date_bgs, ' Begin exposure'            
              d_end = red_Time_conv(date_bgs) + float(conf[4])
              tm = strsplit(string(d_end, format='(f12.6)'),'.',/extract)
              frac = tm[1]
              time = red_Time_conv(long(d_end),/f2s)
              red_fitsaddkeyword, hdr, 'DATE_END', dt[0] + 'T' + time + '.' + frac, ' End exposure'
              lc_state = frame[8]
              red_fitsaddkeyword, hdr, 'FRAMENUM', long(framenum), 'Frame number'
              red_fitsaddkeyword, hdr, 'LCSTATE', fix(lc_state), 'Liquid crystal state'
              if set[2] eq 'polcal' then begin
                qw_state = frame[9]
                lp_state = frame[10]
              endif
                                ; generate filename
              filter1 = fix(burst[9])
              v=execute(fnm_gen)
              if ~v then begin
                print, inam, ': Failed to generate filename \r'
                 print, 'Check the database integrity.'
                 return,0B
              endif               
              if ~strmatch(files_burst, '*'+fnm) then continue
              red_fitsaddkeyword, hdr, 'FILENAME', fnm
              status.add,st               
              headers.add,hdr
              framenumbers.add, framenum
              date_beg.add, date_bgs
            endfor              ; iframe            
          endelse               ; CRISP / CHROMIS

        endfor                  ; ifile
      endfor ;icam
    endfor   ; iset
  endfor ;ins

  if Nfls eq 1 then begin
    status = st.code
    framenumbers = framenum
    date_beg = date_bgs
    return, hdr
  endif
  
  return, headers
end

