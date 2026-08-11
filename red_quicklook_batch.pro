; docformat = 'rst'

;+
; Create quicklooks for a range of dates.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Oleksii Andriienko, ISP   
; 
; :Params:
; 
;    dates : in, type = str
; 
;       Dates in format YYYY-MM-DD:MM-DD or YYYY-MM-DD 
; 
; :Keywords:
; 
;    work_dir : in, optional, type=string
;   
;       Base directory where directories for quicklooks were created.
;
;    data_dir : in, optional, type=string
;
;       Root directory for data (default: '/data')
;
;    instruments : in, optional, type=string array
;
;       Instruments' names to be used for quicklook movies.
;
;    choose_states : in, optional, type=boolean
;
;        Set this keyword to choose spectral points to be used for
;        quicklook movies. If not set then 'core_and_wings'
;        keyword will be used in a call of
;        red__quicklook.pro. Preferably use 'choose_states'
;        if you are sure that same spectral sequences were
;        used for the whole range of observational dates.
;
;    overwrite : in, optional, type=boolean
;
;       Overwrite existing quicklooks.
;
;    min_nscan : in, optional, type=integer
;
;       Minimal number of scans in a data set.
;
;    no_db : in, optional, type=boolean
;
;       Set this keyword to prevent using database to get metadata.
;
;    nthreads : in, optional, type=integer
;
;       Number of threads.
;
;    do_wb : in, optional, type=boolean
;
;       Set this keyword to generate wb quicklook movies.
; 
; 
; :History:
; 
;     2020-06-20 : OA. First version.
;
;     2022-03-25 : OA. Added calls to red__setupworkdir. Added
;                  data_dir, choose_states, do_wb and no_db keywords.
; 
;-
pro red_quicklook_batch, dates $
                         , work_dir = work_dir $
                         , data_dir = data_dir $
                         , instruments = instruments $
                         , choose_states = choose_states $
                         , overwrite = overwrite $
                         , min_nscan = min_nscan $
                         , no_db = no_db $
                         , nthreads = nthreads $
                         , do_wb = do_wb

  if keyword_set(choose_states) then core_and_wings = 0B else core_and_wings = 1B
  if ~keyword_set(nthreads) then nthreads = 20  
  if ~keyword_set(data_dir) then data_dir = '/data'
  if ~file_test(data_dir) then begin
    print,"There are no data in ", data_dir,'.'
    print, 'Please provide a correct data directory.'
    return
  endif
  if ~keyword_set(instruments) then instruments = ['CHROMIS','CRISP'] else instruments = [instruments]
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
      root_dir = '/data/' + year + '/' + year + '-' + month + '/' + date + '/'
      if ~file_test(root_dir) then begin
        print,'There is no data for ', date, '.'
        continue
      endif

      dir = work_dir + date + '/'    
      file_mkdir, dir
      cd, dir
      for ii=0,n_elements(instruments)-1 do begin
        if file_test(dir + instruments[ii]) then continue $
        else $ 
          red_setupworkdir, instruments = instruments[ii] $
                          , search_dirs = root_dir $
                          , out_dir = dir $
                          , /no_observer_metadata 
      endfor

      if file_test(dir + 'CHROMIS') and strmatch(instruments,'CHROMIS')then begin
        cd, dir + 'CHROMIS'
        a = chromisred("config.txt", no_db = no_db)
        obs_dirs = file_search(root_dir + 'CHROMIS-data/*', count=Ndirs)
        if Ndirs eq 0 then begin
          print,'There are no data for CHROMIS on ', date,'.'
          cd,'..'
          spawn, 'rm -r CHROMIS'
          print,'The directory has been removed.'
          goto, skip_chromis
        endif else begin
          if Ndirs eq 1 then begin
            qq = strsplit(obs_dirs,'/',/extract)
            time_stamps = [qq[-1]]
          endif else begin
            qq = strsplit(obs_dirs,'/',/extract)
            zz = qq.ToArray(dimension=0)
            time_stamps = zz[*,-1]
          endelse
        endelse

        for idir = 0,Ndirs-1 do begin
          tdir = dir + time_stamps[idir]
           
          a->quicklook, datasets = time_stamps[idir] $
                        , core_and_wings = core_and_wings $
                        , use_states = chromis_states $
                        , /destretch $
                        , /derotate $
                        , /neuralnet $
                        , min_nscan = min_nscan $
                        , /no_plot_r0 $
                        , overwrite = overwrite $
                        , nthreads = nthreads $
                        , format = 'mov' $
                        , cam = 'Chromis-N'

          if keyword_set(do_wb) then $
            a->quicklook, datasets = time_stamps[idir] $
                        , core_and_wings = core_and_wings $
                        , use_states = chromis_wb_states $
                        , /destretch $
                        , /derotate $
                        , /neuralnet $
                        , /cube_save $
                        , min_nscan = min_nscan $
                        , /no_plot_r0 $
                        , overwrite = overwrite $
                        , nthreads = nthreads $
                        , format = 'mov' $
                        , cam = 'Chromis-W'

          if keyword_set(core_and_wings) then undefine, chromis_states, chromis_wb_states
        endfor        

        cd, '..'
        ff = file_search(dir + 'CHROMIS/quicklook/*', count=Nq)
        if Nq eq 0 then begin
          spawn, 'rm -r CHROMIS'
          print, 'There are no datasets longer than ', min_nscan, ' for CHROMIS on ', dd,'.'
          print, 'The directory is removed.'
        endif
        undefine,a

      endif else print,"Working directory for ", date, "/CHROMIS has not been set up."

skip_chromis:
      if file_test(dir + 'CRISP') and strmatch(instruments,'CRISP') then begin
        cd, dir + 'CRISP'
        a = crispred("config.txt", no_db = no_db)
        obs_dirs = file_search(root_dir + 'Science/*', count=Ndirs)
        if Ndirs eq 0 then begin
          print,'There are no data for CRISP on ', date,'.'
          cd,'..'
          spawn, 'rm -r CRISP'
          print,'The directory has been removed.'
          goto, skip_crisp
        endif else begin
          if Ndirs eq 1 then begin
            qq = strsplit(obs_dirs,'/',/extract)
            time_stamps = [qq[-1]]
          endif else begin
            qq = strsplit(obs_dirs,'/',/extract)
            zz = qq.ToArray(dimension=0)
            time_stamps = zz[*,-1]
          endelse
        endelse  

        for idir = 0,Ndirs-1 do begin
          tdir = dir + time_stamps[idir]
          
          a->quicklook, datasets = time_stamps[idir] $
                        , core_and_wings = core_and_wings $
                        , use_states = crisp_states $
                        , /destretch $
                        , /derotate $
                        , /neuralnet $
                        , min_nscan = min_nscan $
                        , /no_plot_r0 $
                        , overwrite = overwrite $                        
                        , nthreads = nthreads $
                        , format = 'mov' $
                        , cam = 'Crisp-R'

          if keyword_set(do_wb) then $
            a->quicklook, datasets = time_stamps[idir] $
                        , core_and_wings = core_and_wings $
                        , use_states = crisp_wb_states $
                        , /destretch $
                        , /derotate $
                        , /neuralnet $
                        , /cube_save $
                        , min_nscan = min_nscan $
                        , /no_plot_r0 $
                        , overwrite = overwrite $                        
                        , nthreads = nthreads $
                        , format = 'mov' $
                        , cam = 'Crisp-W'

          if keyword_set(core_and_wings) then undefine, crisp_states, crisp_wb_states
        endfor        

        cd, '..'
        ff = file_search(dir + 'CRISP/quicklook/*', count=Nq)
        if Nq eq 0 then begin
          spawn, 'rm -r CRISP'
          print, 'There are no datasets longer than ', min_nscan, ' for CRISP on ', dd,'.'
          print, 'The directory is removed.'
        endif
        undefine,a

      endif else print,"Working directory for ", date, "/CRISP has not been set up."

skip_crisp:
      cd, '..'
      ff = file_search(dir + '/*', count=Ni)
      if Ni eq 0 then $
        file_delete, dir

    endfor   ;; days
  endfor   ;;months
    
end
