pro red_clean_db_tables

  red_mysql_check, handle
  debug = 0
  
  tables = ['chromis_frames','crisp_frames','crisp_polcal_frames','bursts','configs','datasets', $
    'filename_templates', 'dirname_templates', 'filters','calibrations']
  Ntables = n_elements(tables)
  
  for itable = 0, Ntables-1 do begin
    query = 'DELETE FROM ' + tables[itable] + ';'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
    query = 'alter table ' + tables[itable] + ' auto_increment=0;'
    red_mysql_cmd, handle, query, ans, nl, debug=debug
  endfor
  
end  