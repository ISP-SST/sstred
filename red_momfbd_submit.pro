; docformat = 'rst'

;+
; Submit a momfbd job
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
;   IMNO : in, type=string
;
;     A string with an image number or image number range that momfbd
;     understands. The job will be submitted with
;     momfbd_submit.IMNO.cfg as the config file.
;
;     It is the responsibility of the calling program to make sure
;     IMNO does not clash with running jobs or jobs whose output
;     should not be overwritten.
;
; :Keywords:
;
;    port : in, type=integer
;
;      Submit job to manager running with -p PORT
;
;    cfgfile : in, optional, type=string
;
;      A file name, this file will be copied to
;      momfbd_submit.IMNO.cfg. One of CFGFILE and CFGLINES must be
;      given. 
;
;    cfglines : in, optional, type=strarr
;
;      Lines to append to momfbd_submit.IMNO.cfg. (Or create a new
;      file if needed.) One of CFGFILE and CFGLINES must be given.
;
;    dir : in, optional, type=string
;
;      Working dir is DIR.
;
;    force : in, optional, type=boolean
;
;      If force is set, then submit with -f flag to overwrite existing
;      output.
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
;
;-
pro red_momfbd_submit, imno $
                       , port=port $
                       , cfgfile=cfgfile $
                       , cfglines=cfglines $
                       , dir=dir $
                       , force=force

  thiscfgfile='momfbd_submit.'+strtrim(string(imno),2)+'.cfg'
  jobname='momfbd_submit.'+strtrim(string(imno),2)

  if n_elements(port) ne 0 then begin
     red_momfbd_check, PORT = port, NMANAGERS=nmanagersrunning, NSLAVES = nslavesrunning
     portstring=' -p '+strtrim(string(port),2)+' '
  end else begin
     red_momfbd_check, NMANAGERS=nmanagersrunning, NSLAVES = nslavesrunning
     portstring=' '
  end

  if nmanagersrunning eq 0 then begin
     print, 'ERROR in momfbd_submit: No "manager '+portstring+'" running.'
     retall
  end

  if nslavesrunning eq 0 then begin
     print, 'WARNING in momfbd_submit: No slaves running for "manager '+portstring+'".'
  end

  if n_elements(dir) eq 0 then begin
     dir = './'
  end else begin
     if strmid(dir,strlen(dir)-1,1) ne '/' then dir += '/'
  end

  if n_elements(cfgfile) ne 0 then begin
     ;; Copy indicated config file
     spawn,'cp '+cfgfile+' '+dir+thiscfgfile
     append = 1
  endif else begin
     append = 0
  endelse

  if n_elements(cfglines) ne 0 then begin
     ;; Append to config file
     openw, flun, dir+thiscfgfile, append = append, /get_lun
     for i=0,(size(cfglines,/dim))[0]-1 do begin
        ;; Kludge to get output where momfbd_results.pro expects it.
        ;; Must be fixed to deal properly with multi-object jobs.
        printf,flun,cfglines[i]
        if strmatch(cfglines[i], 'object*') then printf,flun,'  OUTPUT_FILE=momfbd_submit.'+string(imno, format = '(i0)')
     endfor
     free_lun,flun
  end

  if ~file_test(dir+thiscfgfile) then begin
     print, 'ERROR in momfbd_submit: No config info given.'
     retall
  end
  
  logfile = red_strreplace(thiscfgfile, '.cfg', '.log')

  spawnstring='cd '+dir+'; jsub '+portstring+' -cfg '+thiscfgfile+' -name '+jobname+' -n '+strtrim(string(imno), 2)+' -lg '+logfile
  if keyword_set(force) then spawnstring += ' -f '

  print, spawnstring
  spawn, spawnstring

end
 
