; docformat = 'rst'

;+
; Submit a mysql REPLACE command.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;  handle, in, type=integer
;
;      The logical unit of the pipe (opened by openmysql).
; 
;  table, in, type=string
; 
;      The table in which to do the update.
; 
;  fields, in, type=strarr
; 
;      The field (column) names.
;
;  values, in, type=strarr
; 
;      The values.
;
;
; 
; :Keywords:
; 
;   verbose, in, optional, type=boolean
;   
;     Be verbose!
; 
; 
; :History:
; 
;   2018-07-02 : MGL. First version.
; 
; 
;-
pro red_mysql_replace, handle, table, fields, values, verbose = verbose

  if n_elements(fields) ne n_elements(values) then stop 
  
  query = 'REPLACE INTO '+table+' (' + strjoin(fields, ',') + ')' $
          + ' VALUES ('+strjoin(values, ',')+') ;'

  if keyword_set(verbose) then print, query
  
  red_mysql_check, handle, /reopen
  red_mysql_cmd, handle, query, DEBUG=verbose

end
