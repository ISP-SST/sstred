; docformat = 'rst'

;+
; Write metadata for a raw-data file into the database.
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
;   dbinfo, in, type=structarr
; 
;     An array of structs, one per file, with metadata info.
; 
; 
; :Keywords:
; 
;   force, in, optional, type=boolean
;   
;     If metadata for a file is already in the database, force overwriting it.
; 
; 
; :History:
; 
;   2018-07-02 : MGL. First version.
; 
;-
pro red_rawfile2db, dbinfo, force = force

  inam = red_subprogram(/low, calling = inam1)
  
  if n_elements(dbinfo) eq 0 then return
  
  instrument = dbinfo[uniq(dbinfo.instrume, sort(dbinfo.instrume))].instrume
  if n_elements(instrument) gt 1 then stop

  ;; Make sure database is open
  red_mysql_check, SQLhandle

  ;; The DATE-OBS keyword defines the data set (together with the
  ;; instrument). 
  ;; Find unique combos of instrument, date and timestamp. Loop over
  ;; such combos.
  uset = dbinfo[uniq(dbinfo.date_obs, sort(dbinfo.date_obs))].date_obs
  Nsets = n_elements(uset)
  
  for iset = 0, Nsets-1 do begin

    indx = where(dbinfo.date_obs eq uset[iset], Nfiles)

    print
    print, inam+' : Write DB info for ' + uset[iset] + ', ' $
           + strtrim(Nfiles, 2) + ' files.'
    print

    ;; For each such combo, look in dataset table for a corresponding
    ;; line.
    where_condition = 'ts='+red_mysql_quote(strjoin(strsplit(uset[iset],'T',/extract),' '))
    dataset = red_mysql_select(SQLhandle $
                               , column_names = datasets_column_names $
                               , column_types = datasets_column_types $
                               , table = 'datasets' $
                               , select_expression = '*' $
                               , where_condition = where_condition $
                               , count = Nsetlines, /verbose)

    case Nsetlines of
      0 : doupdate = 1                  ; No match
      1 : doupdate = 1                  ; Only the heads.
      2 : doupdate = keyword_set(force) ; It's there, update the if /force.
      else : stop                       ; Could there even be more than one line?
    endcase
    
    if doupdate then begin
      ;; Update the line
      undefine, fields, values

      red_append, fields, 'ts'
      red_append, values, red_mysql_quote(strjoin(strsplit(uset[iset],'T',/extract),' '))

      red_append, fields, 'observer'
      red_append, values, red_mysql_quote(dbinfo[0].observer)

      red_append, fields, 'description'
      red_append, values, red_mysql_quote(dbinfo[0].object $
                                          + ' ' + dbinfo[0].waveband)
      
      ;; Add more fields and values that correspond to columns in the
      ;; table.

      ;; Write the lines into the database
      red_mysql_replace, SQLhandle, 'datasets', fields, values
      
    endif


;    ;; Look in the img_files table for lines that correspond to the
;    ;; files. 
;    where_condition = 'data_id='+strtrim(dataset.id, 2) $
;                      + ' && file_number in (' $
;                      + strjoin(strtrim(ddbinfo[indx[ifile]].framenum,2),',') $
;                      + ')'
;    img_files = red_mysql_select(SQLhandle $
;                                 , column_names = img_files_column_names $
;                                 , column_types = img_files_column_types $
;                                 , table = 'img_files' $
;                                 , select_expression = '*' $
;                                 , where_condition = where_condition $
;                                 , ngood = Nlines, /verbose)
;    print, Nfiles eq Nlines


;    stop                        ; OK so far
    
    for ifile = 0, Nfiles-1 do begin
      
      ;; For each file, look in img_files table. Select lines with
      ;; dataset_id from dataset table.
      where_condition = 'data_set='+red_mysql_quote(strtrim(dataset.ts, 2)) $
                        + ' && filename='+red_mysql_quote(dbinfo[indx[ifile]].filename)
      img_files = red_mysql_select(SQLhandle $
                                   , column_names = img_files_column_names $
                                   , column_types = img_files_column_types $
                                   , table = 'img_files' $
                                   , select_expression = '*' $
                                   , where_condition = where_condition $
                                   , count = Nlines, /verbose)

      case Nlines of
        0 : doupdate = 1                  ; No match
        1 : doupdate = 1                  ; Only the heads.
        2 : doupdate = keyword_set(force) ; It's there, update the if /force.
        else : stop                       ; Could there even be more than one line?
      endcase

      if doupdate then begin
        ;; Update the line
        undefine, fields, values
 
        red_append, fields, 'INSTRUMENT'
        red_append, values, red_mysql_quote(dbinfo[0].INSTRUME)
        
        red_append, fields, 'DATA_SET'
        red_append, values, red_mysql_quote(strtrim(dataset.ts, 2))

        red_append, fields, 'FILENAME'
        red_append, values, red_mysql_quote(dbinfo[0].FILENAME)

        red_append, fields, 'DETECTOR'
        red_append, values, strtrim(red_romannumber(dbinfo[0].detector), 2)

        ;; Add more fields and values that correspond to columns in the
        ;; table.

        ;; Write the fields into the database
        red_mysql_replace, SQLhandle, 'img_files', fields, values 
      endif
      stop
    
      for iframe = 0, Nframes-1 do begin
        
        ;; For each frame, look in img_frames table. Select lines with
        ;; file_id from img_file table.
        where_condition = 'file_id='+strtrim(img_files.id, 2) $
                          + ' && frame_number='+dbinfo[indx[ifile]].framenumber
        img_frames = red_mysql_select(SQLhandle $
                                      , column_names = img_frames_column_names $
                                      , column_types = img_frames_column_types $
                                      , table = 'img_frames' $
                                      , select_expression = '*' $
                                      , where_condition = where_condition $
                                      , ngood = Nlines, /verbose)

        case Nlines of
          0 : stop            
          ;; Something's wrong, should have at least the heads
          1 :  doupdate = 1
          ;; Only the heads. Add a line with info selected based on the
          ;; column heads.
          2 : doupdate = keyword_set(force)
          ;; It's there, update the if /force.
          else : stop
          ;; Could there even be more than one line?
        endcase
        
      endfor                    ; iframe
    endfor                      ; ifile
  endfor                        ; iset

  
end
