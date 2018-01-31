; docformat = 'rst'

;+
; Logs parameter and version info for a pipeline step.
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Keywords:
; 
;   top_info_strings : in, optional, type=strarr
; 
;      Lines of text that will appear near the top of the log file. 
; 
; 
;   mid_info_strings : in, optional, type=strarr
; 
;      Lines of text that will appear in the middle the log file,
;      after the parameters and before the git info. 
; 
; 
;   bottom_info_strings : in, optional, type=strarr
; 
;      Lines of text that will appear at the bottom of the log file. 
; 
; 
;   selfinfo : in, optional, type=strarr
; 
;      Structure of the form { parameter_name1:parameter_value1,
;      parameter_name2:parameter_value2, ...}. Will be pretty printed
;      into the log file. 
; 
;   logfile : in, out, optional, type=string
; 
;      The name (including path) of the log file.
; 
;   add : in, optional, type=boolean
; 
;      Set this to add info to an already existing logfile. The
;      contents of the keywords *_info_strings will be used, in order.
;      (Everything is written after any pre-existing text, though.)
; 
; 
; 
; :History:
; 
;   2013-08-23 : MGL. First version.
; 
;   2013-09-02 : MGL. Use red_timestamp rather than the Coyote
;                Graphics timestamp program.
; 
;   2013-09-11 : MGL. Added keywords "add" and "logfile" to make it
;                possible to add info with a later call. Git info is
;                now written to a separate file. Rearranged calls to
;                scope_varfetch in order to reduce risk of it throwing
;                an error when the variable being fetched is
;                undefined. Use red_findpro rather than findpro.
; 
;   2013-09-12 : MGL. Documentation fix.
; 
;   2014-03-05 : THI. Descriptive info about current commit, diff against master
; 
;-
pro red_writelog $
   , selfinfo = selfinfo $
   , top_info_strings = top_info_strings $
   , mid_info_strings = mid_info_strings $
   , bottom_info_strings = bottom_info_strings $
   , logfile = logfile $
   , add = add

 
  ;; Time-stamp
  time = red_timestamp(/utc)

  ;; The name of the program that wants to be logged
  sctrace = scope_traceback(/structure)
  Nlevels = n_elements(sctrace)
  name = strlowcase(sctrace[(Nlevels-2) >0].routine)

  if keyword_set(add) then begin
     openu, lu, logfile, /get_lun
     printf, lu, ' '
     printf, lu, 'Info added by "'+name+'" at : '+time
     printf, lu, ' '
     if n_elements(top_info_strings) gt 0 then printf, lu, top_info_strings, format='(a0)'
     printf, lu, ' '
     if n_elements(mid_info_strings) gt 0 then printf, lu, mid_info_strings, format='(a0)'
     printf, lu, ' '
     if n_elements(bottom_info_strings) gt 0 then printf, lu, bottom_info_strings, format='(a0)'
     printf, lu, ' '
     free_lun, lu
     return
  endif


  ;; Location of the crispred source files:
  srcdir = file_dirname( routine_filepath("red_writelog"), /mark )

  ;; Open the log file
  logdir = './pipeline-log/'
  file_mkdir, logdir
  logfile = logdir+time+'_'+name+'.log'
  openw, lu, logfile, /get_lun

;  printf, lu, whoami()

  ;; Some initial info
  printf, lu, 'The SST/CRISP reduction pipeline step "'+name+'" was called.'
  printf, lu, 'Started at : '+time
  printf, lu, ' '
  
  ;; Top strings
  if n_elements(top_info_strings) gt 0 then printf, lu, top_info_strings, format='(a0)'

  ;; Find argument names
  anames = routine_info(name,/parameters) 
  ;; Find scope level
;  slevel = scope_level()-1 

  ;; Parameter info
  if anames.num_args gt 0 then begin
     printf, lu, ' '
     printf, lu, '------------------ ARGUMENT INFO -----------------'
     printf, lu, ' '
     for i = 0, anames.num_args-1 do begin
        defined = n_elements(SCOPE_VARFETCH(anames.args[i],level=-1))
        if defined then begin
           printf, lu, anames.args[i]
           printf, lu, 'size() : ', size(SCOPE_VARFETCH(anames.args[i],level=-1))
           tmp = SCOPE_VARFETCH(anames.args[i],level=-1)
           if n_elements(tmp) ne 0 then begin
              if n_elements(tmp) gt 10 then begin
                 printf, lu, tmp[0:9]+'...'
              endif else begin
                 printf, lu, tmp
              endelse
           endif
        endif else begin
           printf, lu, anames.args[i], '  (Not defined.)'
        endelse
        printf, lu, ' '
     endfor
     printf, lu, '---------------------------------------------------'
     printf, lu, ' '
  endif

  if anames.num_kw_args gt 0 then begin
     printf, lu, ' '
     printf, lu, '------------------ KEYWORD INFO -----------------'
     printf, lu, ' '
     for i = 0, anames.num_kw_args-1 do begin
        defined = n_elements(SCOPE_VARFETCH(anames.kw_args[i],level=-1))
        if defined then begin
           printf, lu, anames.kw_args[i]
           printf, lu, 'size() : ', size(SCOPE_VARFETCH(anames.kw_args[i],level=-1))
           tmp = SCOPE_VARFETCH(anames.kw_args[i],level=-1)
           if n_elements(tmp) gt 10 then begin
              printf, lu, tmp[0:9]+'...'
           endif else begin
              printf, lu, tmp
           endelse
        endif else begin
           printf, lu, anames.kw_args[i], '  (Not defined.)'
        endelse
        printf, lu, ' '
     endfor
     printf, lu, '---------------------------------------------------'
     printf, lu, ' '
  endif

  ;; SELF info
;  help, self, /obj, output=selfinfo, level = -1 ; This does not work,
;  therefore needs selfinfo as a keyword to this routine.
  if n_elements(selfinfo) gt 0 then begin
     printf, lu, ' '
     printf, lu, '------------------ SELF INFO -----------------'
     printf, lu, ' '
     printf, lu, selfinfo, format='(a0)'
     printf, lu, '---------------------------------------------------'
     printf, lu, ' '
  endif

  ;; Mid strings
  if n_elements(mid_info_strings) gt 0 then printf, lu, mid_info_strings, format='(a0)'



  ;; Bottom strings
  if n_elements(bottom_info_strings) gt 0 then printf, lu, bottom_info_strings, format='(a0)'

  free_lun, lu


  ;; Git info
  gitlogfile = logdir+time+'_'+name+'.git.log'
  openw, glu, gitlogfile, /get_lun
  printf, glu, ' '
  printf, glu, '------------------ GIT INFO -----------------'
  printf, glu, ' '
  gitcmd = 'cd '+srcdir+'; git describe --always --abbrev=12 --long --dirty=\ \(Modified\)'
  printf, glu, '$ '+gitcmd
  spawn, gitcmd, gitinfo
  printf, glu, gitinfo, format='(a0)'
  printf, glu, ' '
  gitcmd = 'cd '+srcdir+'; git rev-parse master'
  printf, glu, '$ '+gitcmd
  spawn, gitcmd, gitinfo
  printf, glu, gitinfo, format='(a0)'
  printf, glu, ' '
  gitcmd = 'cd '+srcdir+'; git diff master'
  printf, glu, '$ '+gitcmd
  spawn, gitcmd, gitinfo
  printf, glu, gitinfo, format='(a0)'
  printf, glu, ' '
  free_lun, glu

end
