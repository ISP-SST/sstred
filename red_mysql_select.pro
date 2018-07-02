;+
; Submit MySQL SELECT query and get struct with results from a table.
;
; Requires an open connection to MySQL server (established by use of
; openmysql) as well as valid permissions for whatever query or
; command is to be executed.
;
; :Categories:
;
;    SST pipeline
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics.
; 
; :Params:
; 
;    lun : in, type=integer
;
;      The logical unit of the pipe (opened by openmysql).
;
;    table : in, type=string
;
;      The table in which the query is done.
;
;    variables : out
;
;      A list of variables to recieve columns of output. Default type
;      is ascii, but use the format keyword to specify other data
;      types.
;
;    query : in, type=strarr
;
;      Query line(s) to send to pipe. Do not include the "FROM table"
;      part. Specify the table in its parameter. 
; 
; 
; :Keywords:
; 
;    count : out, type=integer
;
;      Number of valid lines found that were read.
;
;    column_names : out, optional, type=strarr
;
;      Column names, all columns.
;
;    column_types : out, optional, type=strarr
;
;      Column types, all columns.
;
;    distinct : in, optional, type=boolean
;
;      Do SELECT DISTINCT, returning only distinct (unique) matches. 
;
;    verbose : in, type=boolean
;
;      Flag turns on debugging output to standard out.
; 
; :History:
; 
;   2018-06-28 : MGL. First version.
; 
;   2018-07-02 : MGL. New keyword distinct.
; 
;-
function red_mysql_select, lun $
                           , column_names = column_names $
                           , column_types = column_types $
                           , count = count $
                           , distinct = distinct $
                           , select_expression = select_expression $
                           , table = table $
                           , verbose = verbose $
                           , where_condition = where_condition 

  inam = red_subprogram(/low, calling = inam1)

  if n_elements(table) eq 0 then table = 'img_frames'

  tab = string(9B)

  ;; Add checks: lun has to be an open pipe, query has to be a string,
  ;; table has to be a string.

  ;; Find column definitions
  red_mysql_cmd, lun, 'SHOW COLUMNS FROM '+table+';', columns_info, Ncols, debug=verbose
  Ncols--                       ; First line is the heads
  if Ncols le 0 then stop
  column_names = strarr(Ncols)
  column_types = strarr(Ncols)
  for icol = 0, Ncols-1 do begin
    split_row = strsplit(columns_info[icol+1],tab,/extract,/PRESERVE_NULL)
    column_names[icol] = split_row[0]
    case 1 of
      strmatch(split_row[1], '*float*')       : column_types[icol] = 'FLOAT'
      strmatch(split_row[1], '*double*')      : column_types[icol] = 'DOUBLE'
      strmatch(split_row[1], 'timestamp*')    : column_types[icol] = 'STRING'
      strmatch(split_row[1], 'tinytext*')     : column_types[icol] = 'STRING'
      strmatch(split_row[1], 'date*')         : column_types[icol] = 'STRING' ; date, datetime ;
      strmatch(split_row[1], '*int*unsigned') : column_types[icol] = 'ULONG'
      strmatch(split_row[1], '*int*')         : column_types[icol] = 'LONG'
      else                                    : column_types[icol] = 'STRING' ; default ;
    endcase
  endfor                        ; icol ;

  ;; Construct and submit the query.
  query = 'SELECT '
  if keyword_set(distinct) then query += 'DISTINCT '
  query += select_expression
  if n_elements(table) gt 0 then $
     query += ' FROM '+table
  if n_elements(where_statement) gt 0 then $
     query += ' WHERE '+where_condition
  query += ';'
  red_mysql_cmd, lun, query, result, Nlines, debug=verbose
  Nlines--                      ; First line is the heads
  if Nlines le 0 then stop

  ;; First digest the column headings (split on tabs)
  heads = strsplit(result[0],tab,/extract,/preserve_null)
  Ncols = n_elements(heads)
  ;; Define the struct
  for icol = 0, Ncols-1 do begin
    ;; Match the searched for column names with the defined ones
    pos = where(column_names eq heads[icol], Nwhere)
    if Nwhere eq 0 then stop
    case column_types[pos] of
      'DOUBLE' : value = 0D
      'FLOAT'  : value = 0.0
      'LONG'   : value = 0L
      'STRING' : value = ''
      'UINT'   : value = uint(0)
      'ULONG'  : value = ulong(0)
    endcase
    if icol eq 0 then begin
      strct = create_struct(heads[icol], value)
    endif else begin
      strct = create_struct(strct, heads[icol], value)
    endelse
  endfor                        ; icol
  
  ;; Make array of structs for the output
  output = replicate(strct, Nlines)
  
  ;; Parse the query results
  for iline = 0, Nlines-1 do begin
    split_row = strsplit(result[iline+1],tab,/extract,/PRESERVE_NULL)
    for icol = 0, Ncols-1 do begin
      ;; Match the searched for column names with the defined ones
      pos = where(column_names eq heads[icol], Nwhere)
      if Nwhere eq 0 then stop
      case column_types[pos] of
        'DOUBLE' : value = double(split_row[icol])
        'FLOAT'  : value = float(split_row[icol])
        'LONG'   : value = long(split_row[icol])
        'STRING' : value = split_row[icol]
        'UINT'   : value = uint(split_row[icol])
        'ULONG'  : value = ulong(split_row[icol])
      endcase
      output[iline].(icol) = value
    endfor                      ; icol
  endfor                        ; iline

  count = Nlines

  return, output
  
end

red_mysql_check, SQLhandle

detectors = red_mysql_select(SQLhandle $
                             , select = '*' $
                             , table = 'detectors' $
                             , count = Ndetectors)

umanufacturers = red_mysql_select(SQLhandle $
                                  , select = 'DISTINCT manufacturer' $
                                  , table = 'detectors' $
                                  , count = Nmanufacturers)
print, umanufacturers.manufacturer


end
