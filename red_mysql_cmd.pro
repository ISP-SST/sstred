; docformat = 'rst'

;+
; Send a command to open database and collect the answer.
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
;    Based on mysqlcmd.pro by Marc W. Buie, Lowell Observatory.
; 
; :Params:
; 
;    answer : out, optional, type=strarr
;
;       Returned information from query (possible string array).
;
;    lun : in, type=integer
;
;       The logical unit of the pipe (opened by red_mysql_open).
;
;    nlines : out, optional, type=integer
;
;       Number of lines of returned information (may be zero).
; 
;    query : in, type=strarr
;
;       Command line(s) to send to pipe.
;
; :History:
; 
;  2002-04-10 : MWB. First version.
; 
;  2018-06-28 : MGL. Renamed "red_mysql_cmd" for inclusion in SST data
;               pipelines. Adaptions to SST style.
; 
;-
pro red_mysql_cmd, lun, query, answer, nlines, debug=debug

  if red_badpar(lun,[2,3],0,caller='mysqlcmd (lun) ') then return
  if red_badpar(query,7,[0,1],caller='mysqlcmd (query) ',default='') then return
  if red_badpar(debug,[0,1,2,3],0,caller='mysqlcmd (DEBUG) ',default=0) then return
  
  ;; Send the query
  for i=0,n_elements(query)-1 do begin
    printf,lun,query[i]
    if debug then print,query[i]
  endfor

  ;; Send a special null query that will return two lines of EOT
  printf,lun,"select 'EOT';"
  ;;   if debug then print,"select 'EOT';"
  flush,lun

  ;; Read from pipe until two lines of EOT are seen.
  line=''
  nlines=0L
  done=0
  repeat begin
    readf,lun,line,format='(a)'
    if debug then print,line
    if line eq 'EOT' then begin
      readf,lun,line,format='(a)'
      if line eq 'EOT' then begin
        break
      endif else begin
        if nlines eq 0 then $
           answer = ['EOT',line] $
        else $
           answer = [answer,'EOT',line]
        nlines=nlines+2
      endelse
    endif else begin
      if nlines eq 0 then $
         answer = line $
      else $
         answer = [answer,line]
      nlines=nlines+1
    endelse
  endrep until done
  
  if debug and nlines gt 0 then begin
    print,'Up to the first 5 lines of the query response'
    for i=0,5<(nlines-1) do print,answer[i]
  endif
  
  if debug then print,strn(nlines),' lines of information returned from query.'
  
end
