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
pro red_remove_from_db, date, sets=sets, all=all

  year = strmid(date,0,4)
  if year ge 2023 then db = 'sst_db_' + year else db = 'sst_db'  
  red_mysql_check, handle, database=db
  debug = 0B
  tab = string(9B)

  dt = strsplit(date,'-')
  if n_elements(dt) ne 3 then begin
    print,'Date should be in YYYY-MM-DT format.'
    return
  endif

  if keyword_set(sets) and keyword_set(all) then begin
    print, "Use only one of the keywords either 'sets' or 'all'."
    return
  endif
  
  if ~keyword_set(sets) then begin
    query='select date_obs from datasets;'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    for iset=1,nl-1 do begin
      if strmatch(ans[iset],'*'+date+'*') then begin
        set = strsplit(ans[iset],' ',/extract)
        red_append,sets, set[1]
      endif
    endfor
    if ~keyword_set(all) then begin
      tmp = red_select_subset(sets, indx = indx, qstring = "Select datasets' timestamps")
      sets = sets[indx]
    endif
  endif
  
  Nsets = n_elements(sets)
  if n_elements(sets) eq 1 then sets = [sets]
  for iset = 0, Nsets-1 do begin
    tt = strsplit(sets[iset],':')
    if n_elements(tt) ne 3 then begin
      print,'Dataset timestamp should be in hh:mm:ss format.'
      return
    endif
    query = 'select * from datasets where date_obs = "' + date + ' ' + sets[iset] +'";'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    if nl eq 1 then begin
      print,'There is no entry for ' + date + ' ' + sets[iset] +' dataset in the database.;'
      continue
    endif 
    ss = strsplit(ans[1],tab,/extract,/preserve_null)
    sets_id = ss[0]
    datatype = ss[2]
    instrument = ss[3]
    query = 'select id from configs where sets_id=' + sets_id + ';'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    if nl eq 1 then begin
      query = 'delete from datasets where id = ' + sets_id + ';'
      red_mysql_cmd, handle, query, ans, nl, debug=debug
      continue
    endif
    config_ids = ans[1:nl-1]    
    
    Nconf = nl-1
    for iconf = 0, Nconf-1 do begin
      config_id = config_ids[iconf]
      query = 'select id from bursts where config_id=' + config_id + ' order by id asc limit 1;'
      red_mysql_cmd, handle, query, ans, nl, debug=debug
      if nl eq 1 then begin
        query = 'delete from configs where id = ' + config_id + ';'
        red_mysql_cmd, handle, query, ans, nl, debug=debug
        continue
      endif else begin
        burst_id_first = ans[1]
        query = 'select id from bursts where config_id=' + config_id + ' order by id desc limit 1;'
        red_mysql_cmd, handle, query, ans, nl, debug=debug
        burst_id_last = ans[1]
        
        case instrument of
          'CHROMIS' : begin
            query = 'delete from chromis_frames where burst_id between ' + burst_id_first + ' and ' + burst_id_last + ';'
            red_mysql_cmd, handle, query, ans, nl, debug=debug
          end
          'CRISP' : begin
            if datatype eq 'polcal' then table = 'crisp_polcal_frames' else table = 'crisp_frames'
            query = 'delete from ' + table + ' where burst_id between ' + burst_id_first + ' and ' + burst_id_last + ';'
            red_mysql_cmd, handle, query, ans, nl, debug=debug
          end
        endcase
      endelse

      query = 'delete from bursts where config_id = ' + config_id + ';'
      red_mysql_cmd, handle, query, ans, nl, debug=debug
      query = 'delete from configs where id = ' + config_id + ';'
      red_mysql_cmd, handle, query, ans, nl, debug=debug
    endfor                      ;iconf    
    
    query = 'delete from datasets where id=' + sets_id + ';'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
  endfor                        ; iset
  
end  
