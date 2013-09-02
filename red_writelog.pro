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
; 
; 
; :Params:
; 
;   name : in, type=string
;   
;      The name of the step. (Typically the name of the subroutine
;      performing the step.) Is used to create the log file name.
;   
;   
; 
; :Keywords:
; 
;  top_info_strings : in, optional, type=strarr
; 
;      Lines of text that will appear near the top of the log file. 
; 
; 
;  mid_info_strings : in, optional, type=strarr
; 
;      Lines of text that will appear in the middle the log file,
;      after the parameters and before the git info. 
; 
; 
;  bottom_info_strings : in, optional, type=strarr
; 
;      Lines of text that will appear at the bottom of the log file. 
; 
; 
;  self_info : in, optional, type=strarr
; 
;      Structure of the form { parameter_name1:parameter_value1,
;      parameter_name2:parameter_value2, ...}. Will be pretty printed
;      into the log file. 
; 
; 
; :history:
; 
;   2013-08-23 : MGL. First version.
; 
;   2013-09-02 : MGL. Use red_timestamp rather than the Coyote
;                Graphics timestamp program.
; 
;-
pro red_writelog $
   , selfinfo = selfinfo $
   , top_info_strings = top_info_strings $
   , mid_info_strings = mid_info_strings $
   , bottom_info_strings = bottom_info_strings

  ;; The name of the program that wants to be logged
  sctrace = scope_traceback(/structure)
  Nlevels = n_elements(sctrace)
  name = strlowcase(sctrace[(Nlevels-2) >0].routine)

  logdir = './pipeline-log/'
  file_mkdir, logdir
 
  ;; Time-stamp
  time = red_timestamp(/utc,/no,11)

  ;; 
  fname = logdir+time+'_'+name+'.log'

  ;; Location of the crispred source files:
  findpro, 'crispred', /NoPrint, dirlist = srcdir

  openw, lu, fname, /get_lun
  
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
        printf, lu, anames.args[i]
        printf, lu, 'size() : ', size(SCOPE_VARFETCH(anames.args[i],level=-1))
        defined = n_elements(SCOPE_VARFETCH(anames.args[i],level=-1))
        if defined then begin
           tmp = SCOPE_VARFETCH(anames.args[i],level=-1)
           if n_elements(tmp) ne 0 then begin
              if n_elements(tmp) gt 10 then begin
                 printf, lu, tmp[0:9]+'...'
              endif else begin
                 printf, lu, tmp
              endelse
           endif
        endif
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
        printf, lu, anames.kw_args[i]
        printf, lu, 'size() : ', size(SCOPE_VARFETCH(anames.kw_args[i],level=-1))
        defined = n_elements(SCOPE_VARFETCH(anames.kw_args[i],level=-1))
        if defined then begin
           tmp = SCOPE_VARFETCH(anames.kw_args[i],level=-1)
           if n_elements(tmp) gt 10 then begin
              printf, lu, tmp[0:9]+'...'
           endif else begin
              printf, lu, tmp
           endelse
        endif
        printf, lu, ' '
     endfor
     printf, lu, '---------------------------------------------------'
     printf, lu, ' '
  endif

  ;; SELF info
;  help, self, /obj, output=selfinfo, level = -1 ; This does not work
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


  ;; Git info
  printf, lu, ' '
  printf, lu, '------------------ GIT INFO -----------------'
  printf, lu, ' '
  gitcmd = 'cd '+srcdir+'; git rev-parse HEAD'
  printf, lu, '$ '+gitcmd
  spawn, gitcmd, gitinfo
  printf, lu, gitinfo, format='(a0)'
  printf, lu, ' '
  gitcmd = 'cd '+srcdir+'; git rev-parse origin'
  printf, lu, '$ '+gitcmd
  spawn, gitcmd, gitinfo
  printf, lu, gitinfo, format='(a0)'
  printf, lu, ' '
  gitcmd = 'cd '+srcdir+'; git diff origin'
  printf, lu, '$ '+gitcmd
  spawn, gitcmd, gitinfo
  printf, lu, gitinfo, format='(a0)'
  printf, lu, ' '
  printf, lu, '---------------------------------------------'
  printf, lu, ' '

  ;; Bottom strings
  if n_elements(bottom_info_strings) gt 0 then printf, lu, bottom_info_strings, format='(a0)'

  free_lun, lu


end
