;+
; Submit MySQL query and get response as vectors of data. 
;
; Submit a simple SQL query to MySQL server, using the connection
; previously opened with openmysql. Retrieve the result into a row of
; variables, much as readcol does.
; 
; Requires an open connection to MySQL server (established by use of
; openmysql) as well as valid permissions for whatever query or
; command is to be executed.
;
; Callin sequence: red_mysql_query,lun,query,[varables...],[format='(a,f,...)']
;
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics.
;
;    Based on mysqlquery.pro by Will Grundy.
; 
; :Params:
; 
;    lun : in, type=integer
;
;      The logical unit of the pipe (opened by openmysql).
;
;    variables : out
;
;      A list of variables to recieve columns of output. Default type
;      is ascii, but use the format keyword to specify other data
;      types.
;
;    query : in, type=strarr
;
;      Command(s) to send to pipe.
; 
; 
; :Keywords:
; 
;    heads : out, type=string
;
;      Array of column heads.
;
;    Ngood : out, type=integer
;
;      Number of valid lines found that were read.
;   
;    format : in, type=string
;
;      Specify format of output variables (default is ascii).
;
;    verbose : in, type=boolean
;
;      Flag turns on debugging output to standard out.
; 
; :History:
; 
;   2002-02-14 : WG. Adapted from earlier version called mysql.pro
;
;   2002-02-27 : WG. Changed behavior so 'NULL' becomes NaN instead of
;                making the line be ignored when it occurs in a
;                numerical field.
;
;   2002-03-25 : WG. Changed to split on tab instead of white space,
;                so that strings with internal spaces (but not tabs)
;                can be retrieved.
;
;   2003-01-10 : MWB. Fixed multi-line query error. Only one query per
;                call is allowed.
;
;   2006-11-09 : PLC. Changed strsplit call to use /preserve_null flag
;                this means a field can now return an empty string.
;                The behavior for non-'a' type fields is not defined
;                in this case.
;
;   2007-07-15 : MWB. Added NGOOD output keyword
;
;   2010-03-11 : MWB. Change behavior so that if there is only one
;                valid line the result is returned as a scalar instead
;                of a one-element vector.
;
;   2014-03-07 : MWB. Increased number of returned variables to 35.
;
;   2016-11-14 : MWB. Fix to LL format code handling, propagated from
;                readcol
; 
;   2018-06-28 : MGL. Renamed "red_mysql_query" for inclusion in SST
;                data pipelines. Adaptions to SST style. Removed
;                keyword cmd.
; 
;-
pro red_mysql_query, lun, query $
                     , v1, v2, v3, v4, v5, v6, v7, v8, v9, v10 $
                     , v11, v12, v13, v14, v15, v16, v17, v18, v19, v20 $
                     , v21, v22, v23, v24, v25, v26, v27, v28, v29, v30 $
                     , v31, v32, v33, v34, v35, $
                     format=fmt, heads=heads, verbose=verbose, ngood=Ngood

  ;; Zero the output variables, so there's no risk of accidentally
  ;; re-processing a previous result if something went wrong.
  heads = ''
  Ncol = n_params() - 2
  if Ncol gt 0 then begin
    vv = 'v' + strtrim( indgen(Ncol)+1, 2)
    for k=0,Ncol-1 do begin
      st = vv[k] + ' = ""'
      tst = execute(st)
    endfor
  endif

  ;; lun and query string are mandatory, as are at least one output
  ;; variable
  if n_params() lt 3 then begin
    print,'Usage: red_mysql_query,lun,query,v1,[v2,v3,v4,...]'
    return
  endif

;   query = gettok(query,';')  ; keep only up to 1st ";"
;   red_mysql_cmd,lun,query+';',result,Nlines,debug=verbose
  red_mysql_cmd, lun, query, result, Nlines, debug=verbose

  ;; First digest the column headings (split on tabs)
  heads = strsplit(result[0],'	',/extract,/PRESERVE_NULL)
  ;; Next process everything else into a series of variables, using
  ;; code shamelessly lifted from the astro library's readcol.pro.
  ;; (thank you kindly, Landsman et al.!)
  nskip = 0
  if N_elements(fmt) gt 0 then begin ;FORMAT string supplied?
    ;; Grind format string into usable form
    zparcheck, 'MYSQL', fmt, 2, 7, 0, 'FORMAT string'
    frmt = strupcase(strcompress(fmt,/REMOVE))
    remchar, frmt, '('
    remchar, frmt, ')'
    pos = strpos(frmt, 'X', 0)
    while pos ne -1 DO begin
      pos = strpos( frmt, 'X', pos+1)
      nskip = nskip + 1
    endwhile
  endif else begin
    ;; Default is ascii format (least likely to fail)
    frmt = 'A'
    if Ncol gt 1 then for i = 1,Ncol-1 do frmt = frmt + ',A'
  endelse
  Nfmt = Ncol + nskip
  idltype = intarr(Nfmt)
  ;; Create output arrays according to specified formats
  k = 0L
  hex = bytarr(Nfmt)
  for i = 0L, Nfmt-1 DO begin
    fmt1 = gettok( frmt, ',' )
    if fmt1 eq '' then fmt1 = 'A' ; Default is ascii format
    case strmid(fmt1,0,1) of
      'A': idltype[i] = 7
      'D': idltype[i] = 5
      'F': idltype[i] = 4
      'I': idltype[i] = 2
      'B': idltype[i] = 1
      'L': idltype[i] = strmid(fmt1,0,2) EQ 'LL' ? 14 : 3
      'U': if strmid(fmt1,1,1) NE 'L' then idltype[i] = 12 else $
         idltype[i] = strmid(fmt1,2,1) EQ 'L' ? 15 : 13
      'Z': begin 
        idltype[i] = 3          ;Hexadecimal
        hex[i] = 1b
      end
      'X':  idltype[i] = 0
      else: message,'Illegal format '+fmt1+' in field '+strtrim(i,2)
    endcase
    ;; Define output arrays
    if idltype[i] ne 0 then begin
      st = vv[k] + '= make_array(Nlines,TYPE = idltype[i] )'
      tst = execute(st)
      k = k+1
    endif
  endfor    
  Ngood = 0L
  temp = '	'

  for iline=1L,Nlines-1 DO begin ; Skip first line (headers)
    k = 0
    temp = strtrim(result[iline],1)
    var = strsplit(temp,'	',/extract,/PRESERVE_NULL)
    for ifmt = 0L,Nfmt-1 DO begin
      if ( idltype[ifmt] ne 0 ) then begin ;Expecting data?
        if ifmt+1 gt n_elements(var) then begin
          Ngood=Ngood-1
          goto, badline
        endif
        if ( idltype[ifmt] ne 7 ) then begin               ; Check for numeric data
          tst = strnumber(var[ifmt],val,hex=hex[ifmt])     ; Valid number?
          ;; Instead of failing on 'NULL', need to return 'NaN'
          if strmatch(var[ifmt],'NULL') then begin
            var[ifmt] = 'NaN'
            tst = 1
          endif
          if tst eq 0 then begin ;If not, skip this line
            Ngood = Ngood-1
            goto, BADLINE
          endif
          st = vv[k] + '[Ngood] = val'
        endif else st = vv[k] + '[Ngood] = strtrim(var[ifmt],2)'
        tst = execute(st)
        k = k + 1
      endif
    endfor                      ; ifmt

    BADLINE: Ngood = Ngood+1

  endfor                        ; iline
  
  if Ngood eq 0 then begin
    message,'ERROR - No valid lines found for specified format',/INFORM
    return
  endif else begin
    ;; Compress arrays to match actual number of valid lines
    for i = 0, Ncol-1 do tst = execute(vv[i] + '=trimrank('+ vv[i]+ '[0:Ngood-1])')
  endelse
  
end
