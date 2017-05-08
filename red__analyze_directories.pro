; docformat = 'rst'

;+
; Analyze the data directories for a certain day. Both science data
; and various calibration data.
;
; Note, this method will overwrite the results of previous runs, in
; case data have been added.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;     2016-09-06 : MGL. First version, based on parts of red_plot_r0.
; 
;     2016-09-07 : MGL. Speed up. Remove unused code.
; 
;     2016-09-13 : MGL. Check for empty data directories.
;
; 
;-
pro red::analyze_directories, force = force
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 

  red_writelog, selfinfo = selfinfo


  root_dir = self.root_dir
  isodate = self.isodate
  analysis_dir = self.out_dir+'dir-analysis/'+isodate+'/'
  data_dirs = *self.data_dirs

  timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
  dateregex = '20[0-2][0-9][.-][01][0-9][.-][0-3][0-9]'

  print, inam + ' : Finding the timestamp directories.' 
  timedirs = red_find_matching_dirs(timeregex, rootdir = root_dir, count = Ntd)
  if Ntd gt 0 then timedirs = red_strreplace(timedirs,root_dir,'')

  file_mkdir, analysis_dir

  openw, flun, /get_lun, analysis_dir+'timedirs.txt'
  printf, flun, root_dir, ' '+strtrim(Ntd, 2)
  printf, flun, red_strreplace(timedirs,root_dir,''), format='(a0)'
  free_lun, flun

  
  print, inam + ' : Find start and stop times for each timestamp directory'
  for idir = 0, Ntd-1 do begin

     intfile = analysis_dir + timedirs[idir] + '/interval.txt'
     if keyword_set(force) or ~file_test(intfile) then begin

        ;; The current timestamp directory
        tdir = root_dir+'/'+timedirs[idir]
        tdir = red_strreplace(tdir,'/./','/')
        tdir = red_strreplace(tdir,'//','/')

        ;; Camera directories, pick the wideband camera
        cdirs = file_basename(file_search(tdir+'*', count = Ncams))
        If Ncams gt 0 then cindx = where(strmatch(cdirs,'*-W'), Nw) else Nw = 0
        

        if Nw eq 1 then begin

           cdir = cdirs[cindx[0]]
           
           fnames = file_search(tdir[0]+'/'+cdir+'/*cam*', count = Nf)
           ;; Use spawn rather than file_search in case there are a
           ;; huge number of files in the directory. May be necessary
           ;; for single-frame-per-file CRISP data? Maybe ask the OS
           ;; how many files there are before listing them? Use find
           ;; rather than ls?
;           spawn, 'cd '+tdir[0]+'/'+cdir+' ; ls ', fnames
;           Nf = n_elements(fnames)

           if Nf gt 0 then begin
              
              print, inam + ' : Analyzing ' + timedirs[idir] + ' with '+strtrim(Nf, 2) + ' files.'

              ;; If we got the fnames from file_search, they should
              ;; already be lexically sorted. This is what we need for
              ;; CHROMIS data. But not for CRISP?
;              fnames = red_sortfiles(tdir[0]+'/'+cdir+'/'+fnames)

              head_start = red_readhead(tdir[0]+'/'+cdir+'/'+file_basename(fnames[ 0]))
              head_stop  = red_readhead(tdir[0]+'/'+cdir+'/'+file_basename(fnames[-1]))

              ;; Start and stop time defines the interval
              tstart = red_time2double((strsplit(fxpar(head_start, 'DATE-BEG') $
                                                 ,'T',/extract))[1])
              tstop  = red_time2double((strsplit(fxpar(head_stop,  'DATE-END') $
                                                 ,'T',/extract))[1])
              tinterval = [tstart, tstop]/3600.

              file_mkdir, analysis_dir+timedirs[idir]

              openw, flun, /get_lun, intfile
              printf, flun, tinterval, format = '(f7.4)'
              free_lun, flun
              
           endif                ; Nf
        endif                   ; cdirs
     endif                      ; Actually do it?
  endfor                        ; idir

  ;; Now generate r0 statistics for the scans.
  print, inam + ' : Analyzing the science data.'
  red_logdata, isodate, r0time, r0 = r0data, /use_r0_time; Need the r0 sample times

  for idir = 0, Ntd-1 do begin

     print,strtrim(idir,2)+'/'+strtrim(Ntd,2)
     
     if ~strmatch(timedirs[idir],'*dark*',/fold) $
        and ~strmatch(timedirs[idir],'*flat*', /fold) $
        and ~strmatch(timedirs[idir],'*pinh*', /fold) $
        and ~strmatch(timedirs[idir],'*pfscan*', /fold) $
        and ~strmatch(timedirs[idir],'*polcal*', /fold) then begin

        ;; This is a data directory!

        print, timedirs[idir]

        ;; The current timestamp directory
        tdir = root_dir+'/'+timedirs[idir]
        tdir = red_strreplace(tdir,'/./','/')
        tdir = red_strreplace(tdir,'//','/')
        
        ;; Camera directories, pick the wideband camera
        cdirs = file_basename(file_search(tdir+'*', count = Ncams))
        if Ncams gt 0 then cindx = where(strmatch(cdirs,'*-W'), Nw) else Nw = 0
           
        if Nw eq 1 then begin

           cdir = cdirs[cindx[0]]

           instrument = (strsplit(cdir,'-',/extract))[0] ; CRISP or CHROMIS?

           ;; File with info about the scans in this time-stamped data directory
           scanfile = analysis_dir + timedirs[idir] + '/'+instrument+'-scans.txt'
           
           if keyword_set(force) or ~file_test(scanfile) then begin

              print, inam + ' : Analyse the scans in '+timedirs[idir]
              
              fnames = file_search(root_dir+timedirs[idir]+cdir + '/*cam*', count = Nfile)
              
              if Nfile gt 0 then begin

                 ;; Write the scanfile
                 openw, flun, /get_lun, scanfile
                 
                 self -> extractstates, fnames, states
;nums = filenos, scan = scannos, pref = prefilts
                 
                 indx = sort(states.framenumber)
                 fnames = fnames[indx]
                 states = states[indx]
                 
                 ;; How many different prefilters were used?
                 upref = (states[uniq(states.prefilter, sort(states.prefilter))]).prefilter
                 print, inam + 'WB prefilters: ' + strjoin(upref, ', ') + '.'
                 Npref = n_elements(upref)
                 
                 ;; How many scans (of each kind)?
                 Nscans = max(states.scannumber)+1
                 
                 for iscan = 0L, Nscans-1 do begin
                    for ipref = 0L, Npref-1 do begin
                       
                       print, iscan, Nscans, ipref, Npref, format = '(i0, "/",i0, " ",i0, "/",i0)'
                       
                       scanindx = where(states.scannumber eq iscan and states.prefilter eq upref[ipref], Nscan)

                       ;; Scan start and stop times
                       head_start = red_readhead(fnames[scanindx[ 0]])   
                       head_stop  = red_readhead(fnames[scanindx[-1]])   
                       tstart = red_time2double((strsplit(fxpar(head_start, 'DATE-BEG') $
                                                          ,'T',/extract))[1])
                       tstop  = red_time2double((strsplit(fxpar(head_stop,  'DATE-END') $
                                                          ,'T',/extract))[1])

                       print, iscan, upref[ipref]
                       print, fnames[scanindx[0]]
                       print, fnames[scanindx[-1]]
                       print, tstart, tstop
                       
                       r0indx = where(r0time gt tstart and r0time le tstop, Nr0scan)
                       
                       if Nr0scan gt 0 then begin
                          r0string = string(min(r0data[0, r0indx]), mean(r0data[0, r0indx]) $
                                            , median(r0data[0, r0indx]), max(r0data[0, r0indx]) $
                                            , format = '(f8.5, " ", f8.5, " ", f8.5, " ", f8.5)')
                          if (size(r0data, /dim) )[0] eq 2 then begin
                             r0string += "   " + string(min(r0data[1, r0indx]), mean(r0data[1, r0indx]) $
                                                        , median(r0data[1, r0indx]), max(r0data[1, r0indx]) $
                                                        , format = '(f8.5, " ", f8.5, " ", f8.5, " ", f8.5)')
                          endif

                          
                          r0scanfile = analysis_dir + red_strreplace(timedirs[idir], root_dir, '') $
                                       + '/r0data_scan' + strtrim(iscan, 2) + '_pre' + upref[ipref] $
                                       + '.txt'
                          openw, rlun, /get_lun, r0scanfile
                          for iii = 0, Nr0scan-1 do printf, rlun, r0time[r0indx[iii]], r0data[*, r0indx[iii]]
                          free_lun, rlun

                       endif else r0string = ''
                       
                       printf, flun, string(iscan, upref[ipref], tstart/3600., tstop/3600. $
                                            , format = '(i0, " ", a0, " ", f8.5," ", f8.5)') + r0string
                       
                       
                    endfor      ; ipref
                 endfor         ; iscan

                 free_lun, flun

              endif             ; Any files?                 
              
           endif                ; Actually do it?
        endif                   ; Nw
     endif                      ; Science dir
     
  endfor                        ; idir



end
