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
function red_check_dates_range, dates, range

  if ~strmatch(dates,'*:*') and strlen(dates) ne 10 then begin
    print,'Dates range should be in YYYY-MM-DD:MM-DD format.'
    return,0
  endif
  
  if strlen(dates) eq 10 then begin ;just one date
    dt = strsplit(dates,'-',/extract)
    if n_elements(dt) ne 3 then begin
      print,'Date should be in YYYY-MM-DT format.'
      return,0
    endif
    if strlen(dt[0]) ne 4 or strlen(dt[1]) ne 2 or strlen(dt[2]) ne 2 then begin
      print,'Date should be in YYYY-MM-DT format.'
      return,0
    endif
    year = dt[0]
    month = dt[1]
    day = dt[2]
    dir = '/data/' + year + '/' + year + '-' + month + '/' + year + '-' + month + '-' + day
    if ~file_test(dir) then begin
      print,'There is no data for ', dates, ' day.'
      return,0
    endif
    range = [fix(dt[0]), fix(dt[1]), fix(dt[1]), fix(dt[2]), fix(dt[2])]
    return,1
    
  endif else begin

    yy = strsplit(dates,':',/extract)
    if strlen(yy[0]) ne 10 or strlen(yy[1]) ne 5 then begin
      print,'Date should be in YYYY-MM-DT format.'
      return,0
    endif
    dt = strsplit(yy[0],'-',/extract)
    if n_elements(dt) ne 3 then begin
      print,'Dates range should be in YYYY-MM-DT:MM-DT format.'
      return,0
    endif
    if strlen(dt[0]) ne 4 or strlen(dt[1]) ne 2 or strlen(dt[2]) ne 2 then begin
      print,'Dates range should be in YYYY-MM-DT:MM-DT format.'
      return,0
    endif
    year = dt[0]
    month_beg = fix(dt[1])
    day_beg = fix(dt[2])

    dt = strsplit(yy[1],'-',/extract)
    if n_elements(dt) ne 2 then begin
      print,'Dates range should be in YYYY-MM-DT:MM-DT format.'
      return,0
    endif
    if strlen(dt[0]) ne 2 or strlen(dt[1]) ne 2 then begin
      print,'Dates range should be in YYYY-MM-DT:MM-DT format.'
      return,0
    endif
    month_end = fix(dt[0])
    day_end = fix(dt[1])

    switch 1 of
      month_end lt month_beg : 
      day_beg gt 31 or day_end gt 31 :
      month_beg or month_end gt 12 : begin
        print, 'Dates range is wrong: ', dates, ' Please input correct dates range.'
        return,0
      end
      month_end eq month_beg : begin
        if day_end lt day_beg then begin
          print, 'Dates range is wrong: ', dates, ' Please input correct dates range.'
          return,0
        endif
      end
    endswitch

  endelse

  range = [fix(year),month_beg,month_end,day_beg,day_end]
  return,1
    
end
