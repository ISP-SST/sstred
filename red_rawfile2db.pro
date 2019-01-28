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
  red_mysql_check, handle

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
    where_condition = 'ts='+red_mysql_quote(uset[iset]) ;$
                     ;+ ' && instrument='+red_mysql_quote(instrument)
    dataset = red_mysql_select(handle $
                               , column_names = datasets_column_names $
                               , column_types = datasets_column_types $
                               , table = 'datasets' $
                               , select_expression = '*' $
                               , where_condition = where_condition $
                               , count = Nsetlines)
    
    print,iset, '--> Nsetlines=', strtrim(Nsetlines,2)

    ;; case Nsetlines of
    ;;   0 : stop                          ; Something's wrong, should have at least the heads
    ;;   1 : doupdate = 1                  ; Only the heads.
    ;;   2 : doupdate = keyword_set(force) ; It's there, update the if /force.
    ;;   else : stop                       ; Could there even be more than one line?
    ;; endcase

    case Nsetlines of
      0 : doupdate = 1                  ; Only the heads.
      1 : doupdate = keyword_set(force) ; It's there, update the if /force.
      else : stop                      
    endcase
    
    
    if doupdate then begin
      ;; Update the line
      undefine, fields, values
      
      red_append, fields, 'TS'
      red_append, values, red_mysql_quote(uset[iset])
      
      red_append, fields, 'OBSERVER'
      red_append, values, red_mysql_quote(strtrim('Robustini Esteban', 2))
      
      red_append, fields, 'DESCRIPTION'
      red_append, values, red_mysql_quote(strtrim('Science Ca II H & K', 2))
      
      
      query = 'REPLACE INTO '+'datasets'+' (' + strjoin(fields, ',') + ')' $
              + ' VALUES ('+strjoin(values, ',')+') ;'
      
      print,query
      red_mysql_check, handle
      red_mysql_cmd, handle, query, DEBUG=verbose

      ;;      ;; Add more fields and values that correspond to columns in the
 ;;      ;; table.

 ;;      ;; Write the lines into the database
 ;;      ;red_mysql_replace, handle,'datasets', fields, values 
    endif


print, '--------------------------------------------------------------------------------------------'


    where_condition = 'ts='+red_mysql_quote(uset[iset]) ;$
                     ;+ ' && instrument='+red_mysql_quote(instrument)
    dataset = red_mysql_select(handle $
                               , column_names = datasets_column_names $
                               , column_types = datasets_column_types $
                               , table = 'datasets' $
                               , select_expression = '*' $
                               , where_condition = where_condition $
                               , count = Nsetlines)
    
    print,iset, '--> Nsetlines=', strtrim(Nsetlines,2)

    
      for ifile = 0, Nfiles-1 do begin
        
        
        ;; For each file, look in img_files table. Select lines with
        ;; dataset_id from dataset table.


        where_condition = 'filename='+red_mysql_quote(strmid(dbinfo[ifile].FILENAME,53,93))$
                          + 'AND data_set='+red_mysql_quote(uset[iset])
        
        img_files = red_mysql_select(handle $
                                     , column_names = img_files_column_names $
                                     , column_types = img_files_column_types $
                                     , table = 'img_files' $
                                     , select_expression = '*' $
                                     , where_condition = where_condition $
                                     , count = Nlines)
        
        print,ifile, '--> nfiles=', strtrim(nlines,2)

     
 
        ;; case Nlines of
        ;;   0 : stop                        ; Something's wrong, should have at least the heads
        ;;   1 : doupdate = 1                ; Only the heads.
        ;;   2 : doupdate = keyword_set(force) ; It's there, update the if /force.
        ;;   else : stop                       ; Could there even be more than one line?
        ;; endcase

        
        case Nlines of
          0 : doupdate = 1              ; Only the heads.
          1 : doupdate = keyword_set(force) ; It's there, update the if /force.
          else : stop                   
        endcase
        
        if doupdate then begin
          ;; Update the line
          undefine, fields, values

          time=StrJoin( StrSplit(uset[iset], 'T', /Regex, /Extract, /Preserve_Null), ' ')
          dum = where(time eq dataset.ts)
          
          
          red_append, fields, 'FILENAME'
          red_append, values, red_mysql_quote(strmid(dbinfo[ifile].FILENAME,53,93))
          
          red_append, fields, 'DATA_SET'
          red_append, values, red_mysql_quote(strtrim(dataset[dum].ts, 2))

          
          red_append, fields, 'DETECTOR'
          red_append, values, red_mysql_quote(dbinfo[ifile].detector) ;strtrim(red_romannumber(dbinfo[0].detector), 2)
        
          
          red_append, fields, 'INSTRUMENT'
          red_append, values, red_mysql_quote(dbinfo[ifile].INSTRUME)
          
          
          red_append, fields, 'PREFILTER'
          red_append, values, red_mysql_quote(dbinfo[ifile].PREFILTER)
          
          red_append, fields, 'TUNING'
          red_append, values, red_mysql_quote(dbinfo[ifile].TUNING)

          red_append, fields, 'FRAME_ZERO'
          red_append, values, strtrim((dbinfo[ifile].framenumbers)[0],2)

          red_append, fields, 'TIME'
          red_append, values,red_mysql_quote(strmid(dbinfo[ifile].date,11,11))

          
          
          query = 'REPLACE INTO '+'img_files'+' (' + strjoin(fields, ',') + ')' $
                  + ' VALUES ('+strjoin(values, ',')+') ;'
          
          
          
          red_append, queries, query
 
      
          red_mysql_check, handle;, /reopen 
          red_mysql_cmd, handle, query, DEBUG=verbose
        
        endif

        
      ;; for iframe = 0, Nframes-1 do begin
        
      ;;   ;; For each frame, look in img_frames table. Select lines with
      ;;   ;; file_id from img_file table.
      ;;   where_condition = 'file_id='+strtrim(img_files.id, 2) $
      ;;                     + ' && frame_number='+dbinfo[indx[ifile]].framenumber
      ;;   img_frames = red_mysql_select(SQLhandle $
      ;;                                 , column_names = img_frames_column_names $
      ;;                                 , column_types = img_frames_column_types $
      ;;                                 , table = 'img_frames' $
      ;;                                 , select_expression = '*' $
      ;;                                 , where_condition = where_condition $
      ;;                                 , ngood = Nlines, /verbose)
        
      ;;   case Nlines of
      ;;     0 : stop            
      ;;     ;; Something's wrong, should have at least the heads
      ;;     1 :  doupdate = 1
      ;;     ;; Only the heads. Add a line with info selected based on the
      ;;     ;; column heads.
      ;;     2 : doupdate = keyword_set(force)
      ;;     ;; It's there, update the if /force.
      ;;     else : stop
      ;;     ;; Could there even be more than one line?
      ;;   endcase
        
          ;; endfor                    ; iframe

          
      endfor                    ; ifile
    endfor                      ; iset
    
stop
end
