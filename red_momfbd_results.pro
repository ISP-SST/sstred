; docformat = 'rst'

;+
; Returns results from momfbd job with image numbers IMNO.
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
;     A string with an image number or image number range that momfbd understands.
;
; :Keywords:
;
;   port : in, type=integer
;
;      Wait for manager running with -p PORT to finish before reading
;      output. 
;
;   dir : in, optional, type=string
;
;      Working dir.
;
;   deletecfg : in, optional, type=boolean, default=0
;
;     If set, then delete config file.
;
;   deleteraw : in, optional, type=boolean, default=0
;
;     If set, then delete raw data.
;
;   deleteout : in, optional, type=boolean, default=0
;
;     If set, then delete momfbd output.
;
;   deleteall : in, optional, type=boolean, default=0
;
;     If set, then delete config file, raw data, momfbd output.
;
;   nowait : in, optional, type=boolean, default=0
;
;     If set, don't wait before reading.
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
;
;-
function red_momfbd_results, imno $
                             , port=port $
                             , deletecfg=deletecfg $
                             , deleteraw=deleteraw $
                             , deleteout=deleteout $
                             , deleteall=deleteall $
                             , dir=dir $
                             , nowait = nowait

  if keyword_set(deleteall) then begin
     deletecfg=1
     deleteraw=1
     deleteout=1     
  end
  
  if n_elements(dir) eq 0 then begin
     dir = './'
  end else begin
     if last(dir) ne '/' then dir += '/'
  end

  ;; Wait for momfbd to finish
  repeat begin
     red_momfbd_check, PORT = port, NMANAGERS=nmanagersrunning, NSLAVES = nslavesrunning, NJOBS = njobs
     print, 'njobs=', njobs
     if ~keyword_set(nowait) then wait, 3
  endrep until njobs eq 0


  f0files = file_search(dir+'*.'+string(imno, format = '(i0)')+'.f0', count = Nf0)
  momfbdfiles = file_search(dir+'*.'+string(imno, format = '(i0)')+'.momfbd', count = Nmomfbd)

  if Nf0 eq 0 and Nmomfbd eq 0 then begin
     print,'WARNING in momfbd_submit: No output (yet).'
     return,0
  end

  if Nf0 gt 0 and Nmomfbd gt 0 then begin
     print, 'ERROR in momfbd_submit: Both .f0 and .momfbd files exist for imno=',imno
     stop
  end

  if Nmomfbd gt 1 then begin
     print, 'ERROR in momfbd_submit: Multiple .momfbd files exist for imno=',imno
     stop
  end

  print, 'Reading '+momfbdfiles[0]


  if Nf0 gt 0 then begin        ; ANA format output
     print,'Read all relevant ANA formatted files and make a structure of their contents.'
     if keyword_set(deleteout) then begin ; Delete if so requested
        spawn,'cd '+dir+'; rm -f '+strjoin(f0files,' ')
     end
  end else begin                ; MOMFBD format output
     mr=momfbd_read(momfbdfiles[0])
     if keyword_set(deleteout) then begin
        spawn,'rm '+momfbdfiles[0]
     end   
  end

  if keyword_set(deleteraw) then begin
     print,'Find raw data through cfg file and delete them.'
  end   

  if keyword_set(deletecfg) then begin
     cfgfiledefault='momfbd_submit.'+strtrim(string(imno),2)+'.cfg'
     spawn,'rm -f '+dir+cfgfiledefault
  end   

  return, mr

end
