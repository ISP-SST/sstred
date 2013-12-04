; docformat = 'rst'

;+
; Lists momfbd config files and offers to submit them.
; 
; :Categories:
;
;    SST observations
; 
; :Author:
; 
;    Mats Löfdahl, 2013-12-04
; 
; :Keywords:
; 
;    dir : in, optional, type=string, default='momfbd/'
; 
;      The subdirectory in which to look for config files.
; 
;    port : in, optional, type=integer, default=none
; 
;      The port used for the momfbd manager.
;
;    nthreads : in, optional, type=integer, default="A large fraction of the number of CPUs"
; 
;      The number of threads used by the jsub command.
;
;    dryrun :  in, optional, type=boolean
;
;      Set this to only view the jsub commands without actually
;      submitting anything.
;
; :History:
;
;
;
;
;-
pro red_momfbdjobs, dir = dir, port = port, nthreads = nthreads, dryrun = dryrun

  if n_elements(nthreads) eq 0 then begin
     Ncpu = !cpu.hw_ncpu
     If Ncpu le 2 then Nthreads = 2 else Nthreads = round(Ncpu*.75) 
  endif

  if n_elements(dir) eq 0 then dir = 'momfbd'
  if ~strmatch(dir,'*/') then dir += '/'

  if n_elements(port) eq 0 then begin
     portstring = '' 
  endif else begin
     portstring = ' -p '+string(port)
  endelse

  ;; Find all momfbd config files in dir
  direlements = ['momfbd', '??:??:??', '*', 'cfg'] ; default dir elements
  direlements[0] = strsplit(dir, '/', /extract)    ; elements in given dir
  searchdir = strjoin(direlements, '/')+'/'        
  cfiles = file_search(searchdir+'*.cfg', COUNT = Ncfg)

  ;; Find out the status of jobs corresponding to the config files.
  ;; Use that to build a list of cfg files to select from.
  selectionlist = strarr(Ncfg)
  defaultindx = bytarr(Ncfg)

  ;; Current momfbd job queue
  spawn, 'jstat '+portstring+' -j', joblist

  for i = 0, Ncfg-1 do begin
     
     ;; Index and config file
     sel = cfiles[i]+' '
     
     cdir = file_dirname(cfiles[i])+'/'
     cfile = file_basename(cfiles[i])
     scanno = (strsplit(cfile, '.', /extract))[3]
     rdir = cdir+'results/'

     ;; Currently in queue? Look for momfbd name tag given when submiting.
     name = strjoin((strsplit(cdir,'/',/extract))[1:2],'/')+'/'+strjoin((strsplit(cfile,'.',/extract))[3],'.')
     Nqueue = round(total(strmatch(joblist, name)))
     if Nqueue gt 0 then sel += 'Q:'+strtrim(Nqueue, 2)

     
     ;; Already processed, generated momfbd formatted output?
     mfiles = file_search(rdir+'cam*.'+scanno+'.*.*.lc?.momfbd', COUNT = Nmomfbd)
     if Nmomfbd gt 0 then sel += 'M:'+strtrim(Nmomfbd, 2)

     ;; Already processed, generated ana format output?
     afiles = file_search(rdir+'cam*.'+scanno+'.*.*.lc?.f?', COUNT = Nana)
     if Nana gt 0 then sel += 'A:'+strtrim(Nana, 2)

     ;; As default, suggest processing if not already in queue or
     ;; already processed.
     if Nmomfbd + Nana + Nqueue eq 0 then defaultindx[i] = 1

     selectionlist[i] = sel

  endfor
  
  default = red_collapserange(where(defaultindx),ld='',rd='')
  cfiles = red_select_subset(selectionlist, qstring = 'Submit some jobs?', default = default, count = Nselect)
  print, Nselect
  
  for i = 0, Nselect-1 do begin

     cdir = file_dirname(cfiles[i])+'/'
     cfile = file_basename(cfiles[i])
     lfile = file_basename(cfiles[i], '.cfg', /FOLD_CASE)+'.log'

     ;; The name tag shown in "jsub -j"
     name = strjoin((strsplit(cdir,'/',/extract))[1:2],'/')+'/'+strjoin((strsplit(cfile,'.',/extract))[3],'.')
     cmd = 'cd '+cdir+' ; jsub'+portstring+' -s -v -mt '+strtrim(Nthreads, 2)+' -cfg '+cfile+' -name '+name+' -lg '+lfile
     print, cmd
     if ~keyword_set(dryrun) then spawn, 'cmd'

  endfor

  stop
end
