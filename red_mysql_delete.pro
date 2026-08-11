; docformat = 'rst'

;+
; Delete a row from a table.
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
;   where_condition, in, type=string
;
;      The WHERE condition.
;
; 
; :Keywords:
; 
;   verbose, in, optional, type=boolean
;   
;     Be verbose!
; 
; :History:
; 
;   2018-07-03 : MGL. First version.
; 
; 
;-
pro red_mysql_delete, handle, table, where_condition, verbose = verbose

  query = 'DELETE FROM '+table+' WHERE + '+where_condition+' ;'

  if keyword_set(verbose) then print, query
  
  red_mysql_check, handle
  red_mysql_cmd, handle, query, DEBUG=verbose

end
