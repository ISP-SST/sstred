; docformat = 'rst'

;+
; It checks if there are observations for a given range of dates and creates
; directories for quicklooks.
; 
; :Categories:
;
;    SST pipeline 
; 
; :Author:
; 
;    Oleksii Andriienko, ISP
; 
; :Params:
; 
;    dates : in, type = str
; 
;      Dates in format YYYY-MM-DD:MM-DD or YYYY-MM-DD 
; 
; :Keywords:
; 
;    instruments : in, optional, type=string
;   
;      Instruments for which quicklooks will be created
;
;    out_dir : in, optional, type=string
;
;      The base directory in which directories for quicklooks will be
;      created.
; 
; 
; :History:
; 
;     2020-06-20 : OA. First version.
; 
;-
pro red_setup_quicklook_dirs, dates, instruments=instruments, out_dir=out_dir

  ;; we are at Alba Nova all data should be located in /data

  if ~file_test('/data') then begin
    print,"The data should be located in '/data' directory."
    return
  endif

  if ~keyword_set(out_dir) then out_dir = '/scratch/olexa/quicklook/'
  if ~keyword_set(instruments) then instruments = ['CHROMIS','CRISP']

  if ~red_check_dates_range(dates, range) then return
  year = strtrim(range[0],2)
  month_beg = range[1]
  month_end = range[2]
  day_beg = range[3]
  day_end = range[4]

  if ~file_test('/data/' + year) then begin
    print,'There is no data for ', year, ' year.'
    return
  endif

  for mm = month_beg, month_end do begin
    month = string(mm, format = '(I02)')
    dir = '/data/' + year + '/' + year + '-' + month + '/'
    if ~file_test(dir) then begin
      print,'There is no data for ', year + '-' + month, ' month.'
      continue
    endif

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
      dir = '/data/' + year + '/' + year + '-' + month + '/' + date
      if ~file_test(dir) then begin
        print,'There is no data for ', year + '-' + month + '-' + day, ' day.'
        continue
      endif
      work_dir = out_dir + year + '-' + month + '-' + day
      file_mkdir, work_dir
      cd,work_dir
      red_setupworkdir, instruments = instruments, search_dirs = dir, out_dir = work_dir
      cd,'..'
    endfor
  endfor
    
end
