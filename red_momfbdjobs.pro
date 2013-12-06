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
;    scriptfile : in, optional, type=string, default="momfbdjobs.sh"
;
;      The name of the file in which to write the jsub commands.
;
; :History:
;
;    2013-12-05 : MGL. Removed dryrun keyword. New keyword scriptfile.
;                 Now write jsub commands to a script file instead of
;                 submitting them right away. Now does pushd and popd
;                 for each config file instead of just once before the
;                 whole set of jsub commands.
;
;
;-
pro red_momfbdjobs, dir = dir $
                    , port = port $
                    , nthreads = nthreads $
                    , scriptfile = scriptfile

  if n_elements(nthreads) eq 0 then begin
     Ncpu = !cpu.hw_ncpu
     If Ncpu le 2 then Nthreads = 2 else Nthreads = round(Ncpu*.75) 
  endif

  if n_elements(scriptfile) eq 0 then scriptfile = 'momfbdjobs.sh'

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

     spawn, 'grep DATE_OBS '+cfiles[i]+' | cut -d= -f2', date_obs 

     ;; Currently in queue? Look for momfbd name tag given when submiting.
     name = date_obs + '/' + strjoin((strsplit(cdir,'/',/extract))[1:2],'/') + '/' $
            + strjoin((strsplit(cfile,'.',/extract))[3],'.')
     Nqueue = round(total(strmatch(joblist, name)))
     if Nqueue gt 0 then sel += 'Q:' + strtrim(Nqueue, 2)

     
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
  tmp = red_select_subset(selectionlist, qstring = 'Submit some jobs?' $
                          , default = default, count = Nselect, indx = sindx)

  if Nselect gt 0 then begin

     cfiles = cfiles[sindx]

     openw, slun, scriptfile, /get_lun

     for i = 0, Nselect-1 do begin

        cdir = file_dirname(cfiles[i])+'/'
        cfile = file_basename(cfiles[i])
        lfile = file_basename(cfiles[i], '.cfg', /FOLD_CASE)+'.log'

        spawn, 'grep DATE_OBS '+cfiles[i]+' | cut -d= -f2', date_obs 

        ;; The name tag shown in "jsub -j"
        name = date_obs + '/' + strjoin((strsplit(cdir,'/',/extract))[1:2],'/') + '/' $
               + strjoin((strsplit(cfile,'.',/extract))[3],'.')
        cmd = 'jsub' + portstring + ' -s -v -mt ' + strtrim(Nthreads, 2) $
              + ' -cfg ' + cfile + ' -name ' + name + ' -lg ' + lfile

        print, cmd

        printf, slun, 'sleep 0.5'
        printf, slun, 'pushd '+cdir
        printf, slun, cmd
        printf, slun, 'popd'

     endfor

     free_lun, slun
     spawn, 'chmod a+x '+scriptfile
  endif else begin
     print,  'None selected.'
  endelse

end
