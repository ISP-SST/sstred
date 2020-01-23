; docformat = 'rst'

;+
; Submit jobs to the (redux) momfbd queue.
;
; :Categories:
;
;   Image restoration, mombfd
;
; :Author:
;
;   Mats Löfdahl, Institute for Solar Physics, mats@astro.su.se.
;
; :Params:
;
;    cfgdir : in, type=string
;
;      The path to a directory. The momfbd cfg files in this directory
;      (or a subset thereof) will be submitted to the redux momfbd
;      queue (see the shell command rdx_sub).
;
; :Keywords:
;
;    force : in, optional, type=boolean
;
;       Overwrite existing output.
;
;    options : in, optional, type=string
;
;       Options to rdx_sub, added to the command line.
;
;    no_check : in, optional, type=boolean
;
;       Submit without checking.
;
;    port : in, optional, type=integer
;
;       Specify the port of the redux manager.
;
;    priority : in, optional, type=integer, default=10
;
;       Specify the priority of the job. (Use with good judgement!)
;
;    scannos : in, optional, type="intarr or string"
;
;       An array of scan numbers or a comma/dash-separated string of
;       scan numbers.
;
;    user : in, optional, type=string
;
;       User name to be shown by rdx_stat.
;
;    verbose  : in, optional, type=boolean
;
;       Make the rdx_sub command verbose.
;
; :History:
;
;   2010-09-21 : Written by Mats Löfdahl
;
;   2012-11-21 : Switched to docformat 'rst'.
;
;   2013-08-30 : MGL. Renamed for inclusion in crispred pipeline.
;
;   2013-09-04 : MGL. Use red_momfbd_check, not momfbd_check.
;
;   2014-01-22 : MGL. Adapt to string functions moved to the str_
;                namespace.
;
;   2014-08-12 : MGL. Don't use last().
;
;   2019-10-23 : MGL. Rewrite from scratch, now for running with the
;                redux momfbd code.
; 
;-
pro red_momfbd_submit, cfgdir $
                       , force = force $
                       , options = options $
                       , no_check = no_check $
                       , port = port $
                       , priority = priority $
                       , scannos = scannos $
                       , user = user $
                       , verbose = verbose 

  inam = red_subprogram(/low, calling = inam1)

  if n_elements(cfgdir) eq 0 then begin
    print, inam + ' : No cfg directory given.'
    retall
  endif
  pwd = getenv('PWD')+'/'
  cfgdir = red_strreplace(cfgdir, pwd, '') ; Want relative path
  dir_parts = strsplit(cfgdir, path_sep(), /extract, count = Nparts)

  if n_elements(scannos) eq 0 then begin
    cfgfiles = file_search(cfgdir+'/*.cfg', count = Nfiles)
    if Nfiles eq 0 then begin
      print, inam + ' : No cfg files in directory '+cfgdir
      retall
    endif
  endif else begin
    if size(scannos, /tname) eq 'STRING' then scannos = rdx_str2ints(scannos)
    cfgfiles = file_search(cfgdir+'/*_'+string(scannos, format = '(I05)')+'.cfg', count = Nfiles)
    if Nfiles eq 0 then begin
      print, inam + ' : No cfg files with the given scan numbers in directory '+cfgdir
      retall
    endif
  endelse

  
  cmd_start = 'rdx_sub'
  
  if keyword_set(force)        then cmd_start += ' -f'
  if keyword_set(verbose)      then cmd_start += ' -v'
  if keyword_set(no_check)     then cmd_start += ' --no-check'
  if n_elements(user) gt 0     then cmd_start += ' --user '+ strtrim(user, 2)
  if n_elements(port) gt 0     then cmd_start += ' --port '+ strtrim(port, 2)
  if n_elements(priority) gt 0 then cmd_start += ' --priority '+ strtrim(priority, 2)
  if n_elements(options) gt 0  then cmd_start += ' ' + strtrim(options, 2)

  tstamp = (stregex(cfgdir, '([0-9][0-9]:[0-9][0-9]:[0-9][0-9])', /sub, /extract))[1]
  date = (stregex(pwd, '([0-9][0-9][0-9][0-9][-.][0-9][0-9][-.][0-9][0-9])', /sub, /extract))[1]

  print
  if Nfiles gt 0 then print, inam + ' : Submitting cfg files in '+cfgdir
  
  for ifile = 0, Nfiles-1 do begin

    cmd = cmd_start

    cfgfile = file_basename(cfgfiles[ifile])
    cfgfile = cfgfiles[ifile]
    cmd += ' --config '+cfgfile

    logfile = red_strreplace(cfgfile, '.cfg', '.log')
    cmd += ' --log-file '+logfile
    
    undefine, name_parts 
    pref = (stregex(cfgfile, '_([0-9][0-9][0-9][0-9])_', /sub, /extract))[1]
    scanno = (stregex(cfgfile, '_([0-9][0-9][0-9][0-9][0-9]).cfg', /sub, /extract))[1]
    if date   ne '' then red_append, name_parts, date
    if tstamp ne '' then red_append, name_parts, tstamp
    if pref   ne '' then red_append, name_parts, pref
    if scanno ne '' then red_append, name_parts, scanno
    if n_elements(name_parts) gt 0 then begin
      name = strjoin(name_parts, '_')
      cmd += ' --name '+name
    endif
    
    print, cmd

    stop
    
    spawn, 'cd '+cfgdir+ ' & ' + cmd 
    
  endfor                        ; ifile
    

end
 
