;+
; Check that a connection to mysql is open; open one if it isn't.
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
;    Based on mysqlcheck.pro by Marshall Perrin.
; 
; :Params:
;
;    handle : in, out, type=integer
;
;      The lun for the database pipe.
; 
; :Keywords:
; 
;    database : in, optional, type=string, default='sst_db'
;
;      Database name to open.
; 
; :History:
; 
;   2005-10-05 : MP. First version.
; 
;   2018-06-28 : MGL. Adapted to SST style and renamed
;                red_mysql_check. Set default database.
;
;   2019-10-08 : OA. Added site check to run mysql client on polar
;                remotely from La Palma (it's much faster to populate
;                the database this way).
;
;-
pro red_mysql_check, handle, reopen=reopen, database=database
  
  if n_elements(database) eq 0 then database="sst_db"
  
  if keyword_set(reopen) $
     and n_elements(handle) ne 0 then begin
    free_lun, handle
    ;;close,handle
    ;;delvarx, handle
  endif
  
  if n_elements(handle) gt 0 then return

  red_currentsite,site=site
  case site of
    ;'AlbaNova': red_mysql_open, handle, database
    'La Palma': red_mysql_open_remote, handle, database
    else: red_mysql_open, handle, database
  endcase
  
end
