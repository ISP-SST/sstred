; docformat = 'rst'

;+
; Move raw files from selected scans to a temporary directory for
; deletion.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Params:
; 
;   tdir : in, type=string
; 
;      The timestamp of the directory.
; 
; 
; :Keywords:
; 
;   doit : in, optional, type=boolean
;   
;      Without setting this, nothing will change on disk.
; 
; 
; :History:
; 
; 
; 
;-
pro red::selective_delete, tdir, doit = doit

  cams = *self.cameras
  instrument = (strsplit(cams[0], '-', /extract))[0]

  Ndir = n_elements(tdir)
  
  if ~file_test(tdir) then begin
    datadirs = file_search(self.root_dir + '/' + instrument.toupper() + '-*/' $
                           + tdir + '/' + cams, count = Ndir)
  endif 

  if Ndir eq 0 then stop

  deldirs = red_strreplace(datadirs, self.root_dir, self.root_dir+'DELETE/')
  files = file_search(datadirs[0]+'/*', count = Nfiles)

  ;; Quickly get the scan numbers
  red_extractstates, files, /basename, scan = scannos
  uscan = scannos[uniq(scannos, sort(scannos))]
  Nscans = n_elements(uscan)
  
  ;; We need the times for the first and last files of each scan
  beg_indx = lonarr(Nscans)
  end_indx = lonarr(Nscans)
  beg_time = dblarr(Nscans)
  end_time = dblarr(Nscans)
  for iscan = 0, Nscans-1 do begin
    indx = where(scannos eq uscan[iscan])
    beg_indx[iscan] = indx[0]
    end_indx[iscan] = indx[-1]
  endfor                        ; iscan

  beg_time = red_time2double(strmid(red_fitsgetkeyword_multifile(files[beg_indx], 'DATE-BEG'),11))
  end_time = red_time2double(strmid(red_fitsgetkeyword_multifile(files[end_indx], 'DATE-END'),11))

  ;; Download the r0 log file 
  red_logdata, self.isodate, r0time, r0 = r0data, /use_r0_time

  dt_log = median(red_differential(r0time))
  dt_data = median(end_time-beg_time)

  if dt_data lt dt_log*1.1 then stop
  
  r0_median = fltarr(Nscans)
  r0_mean   = fltarr(Nscans)
  r0_min    = fltarr(Nscans)
  r0_max    = fltarr(Nscans)
  for iscan=0, Nscans-1 do begin
    indx = where(r0time ge beg_time[iscan] and r0time le end_time[iscan])
    r0_median[iscan] = median(r0data[1, indx])
    r0_mean[iscan]   = mean(r0data[1, indx])
  endfor

  selstring = strtrim(long(uscan), 2) + ' : ' + string(r0_mean*100, format = '(f4.1)') + ' , ' + string(r0_median*100, format = '(f4.1)')

  print
  print, 'Scan# : r0 , ground r0 [cm]'
  print,selstring,format='(a0)'
  print
  print, 'See also r0 plots from quicklook.'
  print

  red_strflow, 'Type selection to keep (like K1-5,10-20,23,27) or delete (D1-5,10-20,23,27):'
  s = ''
  read, s

  ;; Do some checking on s first. If "D*", then delete all scans. If
  ;; empty, "K*", or beginning with anything but "D" or "K", then
  ;; return without deleting anything. Write appropriate messages
  ;; before returning.

  ;; If "D" or "K" followed by anything else but "*" or a string that
  ;; can be expanded by rdx_str2ints(), then complain but don't delete.

  case strupcase(strmid(s, 0, 1)) of
    'D' :  begin
      dscan = rdx_str2ints(strmid(s, 1))
      Ndel = n_elements(dscan)
      match2, long(uscan), long(dscan), suba, subb
      kindx = where(suba eq -1, Nkeep)
      if Nkeep gt 0 then kscan = long(uscan[kindx])
    end
    'K' : begin
      kscan = rdx_str2ints(strmid(s, 1))
      Nkeep = n_elements(kscan)
      match2, long(uscan), long(kscan), suba, subb
      dindx = where(suba eq -1, Ndel)
      if Ndel gt 0 then dscan = long(uscan[dindx])
    end
    else : stop
  endcase

  print, 'Delete these scans: ', dscan
  print, 'Keep these scans: ', kscan

  ;; Actually, don't delete. Move to temporary directories on the same
  ;; disk.
  if keyword_set(doit) then file_mkdir, deldirs

  for idir = 0, Ndir-1 do begin

    if idir gt 0 then begin
      files = file_search(datadirs[idir]+'/*', count = Nfiles)
      red_extractstates, files, /basename, scan = scannos
    endif
    
    ;; filenames to delete
    match2, long(scannos), dscan, suba, subb
    indx = where(suba ne -1)
    dfiles = files[indx]

    if keyword_set(doit) then file_move, dfiles, deldirs[idir]
    
  endfor                        ; idir

  print, 'Moved files from these directories:'
  hprint, datadirs
  print, 'to these temporary directories:'
  hprint, deldirs
  print
  red_strflow, 'Please inspect the result and delete the temporary directories when satisfied.'

  if ~keyword_set(doit) then begin
    print
    red_strflow, 'This was a test run. Set keyword /doit to actually create temporary directories and move the files.'
  endif
  
endif

a = crisp2red(/dev, /no)

a -> selective_delete, '08:47:50'

end
