; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;   
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
; 
; 
; 
; 
; 
;-
pro red_quicklook_batch, dates, work_dir = work_dir, overwrite = overwrite, min_nscan = min_nscan

  if ~keyword_set(work_dir) then work_dir = '/scratch/olexa/quicklook/'
  if ~keyword_set(min_nscan) then min_nscan = 5
  if ~red_check_dates_range(dates, range) then return
  year = strtrim(range[0],2)
  month_beg = range[1]
  month_end = range[2]
  day_beg = range[3]
  day_end = range[4]  
  
  for mm = month_beg, month_end do begin
    month = string(mm, format = '(I02)')    

    if mm lt month_end then begin
      if mm eq 2 then begin
        if (fix(year) mod 4) eq 0 then d_end = 29 else d_end = 28 
      endif else begin
        if (mm mod 2) eq 0 then begin
          if mm eq 8 then d_end = 31 else d_end = 30
        endif else d_end = 31
      endelse
    endif else d_end = day_end
    if mm eq month_beg then d_beg = day_beg else d_beg = 1
    
    for dd = d_beg, d_end do begin
      day = string(dd,format = '(I02)')
      date = year + '-' + month + '-' + day
      dir = work_dir + date + '/'
      if ~file_test(dir) then begin
        print,"Working directory for ", date, " day has not been set up or there were no observations on the day. Run 'red_setup_quicklook_dirs' first."
        continue
      endif

      root_dir = '/data/' + year + '/' + year + '-' + month + '/' + date + '/'
      nthreads=20
      if file_test(dir + 'CHROMIS') then begin
        cd, dir + 'CHROMIS'
        a = chromisred("config.txt", /develop)
        obs_dirs = file_search(root_dir + 'CHROMIS-data/*')
        Ndirs = n_elements(obs_dirs)
        if  Ndirs eq 0 then begin
          print, 'There is no data for ', date, ' CHROMIS.'
          continue
        endif
        if  Ndirs eq 1 then begin
          qq = strsplit(obs_dirs,'/',/extract)
          time_stamps = [qq[-1]]
        endif
        if Ndirs gt 1 then begin
          qq = strsplit(obs_dirs,'/',/extract)
          zz = qq.ToArray(dimension=0)
          time_stamps = zz[*,-1]
        endif

        for idir = 0,Ndirs-1 do begin
          tdir = dir + time_stamps[idir]
          if file_test(tdir) and ~keyword_set(overwrite) then begin
            print, 'Quicklook directory ', time_stamps[idir], ' exists. Use /overwrite to redo it.'
            continue
          endif
          a->quicklook, datasets = time_stamps[idir], /core_and_wings, /destretch, /neuralnet, $
                        /cube_save, min_nscan = min_nscan, /no_plot_r0, overwrite = overwrite
        endfor
        cd, '..'
        undefine,a
      endif else print, print,"Working directory for ", date, "/CHROMIS has not been set up. Run 'red_setup_quicklook_dirs' first."

      if file_test(dir + 'CRISP') then begin
        cd, dir + 'CRISP'
        a = crispred("config.txt", /develop)
        obs_dirs = file_search(root_dir + 'Science/*')
        Ndirs = n_elements(obs_dirs)
        if  Ndirs eq 0 then begin
          print, 'There is no data for ', date, ' CHROMIS.'
          continue
        endif
        if  Ndirs eq 1 then begin
          qq = strsplit(obs_dirs,'/',/extract)
          time_stamps = [qq[-1]]
        endif
        if Ndirs gt 1 then begin
          qq = strsplit(obs_dirs,'/',/extract)
          zz = qq.ToArray(dimension=0)
          time_stamps = zz[*,-1]
        endif
        
        for idir = 0,Ndirs-1 do begin
          tdir = dir + time_stamps[idir]
          if file_test(tdir) and ~keyword_set(overwrite) then begin
            print, 'Quicklook directory ', time_stamps[idir], ' exists. Use /overwrite to redo it.'
            continue
          endif
          a->quicklook, datasets = time_stamps[idir], /core_and_wings, /destretch, /neuralnet, $
                        min_nscan = min_nscan, /no_plot_r0, overwrite = overwrite ;, /cube_save
        endfor
        cd, '..'
        undefine,a
      endif else print, print,"Working directory for ", date, "/CRISP has not been set up. Run 'red_setup_quicklook_dirs' first."

    endfor   ;; days
  endfor   ;;months
    
end
