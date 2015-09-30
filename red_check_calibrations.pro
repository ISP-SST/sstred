; docformat = 'rst'

;+
; Warn if any calibration data is missing or incomplete.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2015-09-01
; 
; 
; :Params:
;
;    root_dir : in, type=string
;
;       The date-directory of the observations.
; 
; :Keywords:
;
;       all : in, optional, type=boolean, default=FALSE
;
;          Set this to make all of the darks, flats, polcal, and
;          pinholes keywords TRUE.
;
;       darks : in, optional, type=boolean, default=FALSE
;
;          Set this to check for darks data.
;
;       flats : in, optional, type=boolean, default=FALSE
;
;          Set this to check for flats data.
;
;       polcal : in, optional, type=boolean, default=FALSE
;
;          Set this to check for polcal data.
;
;       pinholes : in, optional, type=boolean, default=FALSE
;
;          Set this to check for pinholes data.
;
;       logfile : in, optional, type=string, default='check_calibrations_TIMESTAMP.txt'
;
;          The name of the file where the output will be logged. The
;          string TIMESTAMP in the default will be replaced with an
;          actual timestamp
;
; :History:
;
;    2015-09-01 : MGL. Looking for different kinds of data ready.
;                 Checking science data for corresponding dark frame
;                 data ready.
;
;    2015-09-03 : MGL. Checking science data for corresponding flat
;                 field data ready.
;
;    2015-09-04 : MGL. Checking science data for corresponding polcal
;                 data ready. Checking science data for corresponding
;                 pinhole data ready. 
;
;    2015-09-07 : MGL. Made the CRISPRED method into an ordinary
;                 subroutine. A skeleton method will get info from the
;                 config file and call this subroutine.
;
;    2015-09-18 : MGL. Include a timestamp in the default logfile
;                 name. Write error status at end of execution. Drop
;                 science frames that match '*.lcd.*'.
;
;    2015-09-21 : MGL. Drop frames of all kinds that match '*.lcd.*'.
;                 Catch case when there are no calibration data of a
;                 certain kind.
;
;    2015-09-28 : MGL. Variables that are 'uniq'-ued now all start
;                 with the letter 'u'. Now checks polcal in each
;                 CRISP-[TR] science directory for one prefilter at a
;                 time.
;
;    2015-09-30 : MGL. Fix a typo and work around strmatch() not
;                 thinking 'pinholecamstates = []' makes
;                 pinholecamstates defined. 
;                 
;
; 
;-
pro red_check_calibrations, root_dir $
                            , all = all $
                            , darks = darks $
                            , flats = flats $
                            , polcal = polcal $
                            , pinholes = pinholes $
                            , logfile = logfile

; The logic of this program is as follows:
;
; * First look for calibrations data of the specified types. Report if
;   there are none.
;
; * Then loop through the science data directories (essentially all
;   directories that are not calibrations data). For the specified
;   calibrations data types, see if they are required by the data and,
;   if so, they were found in the previous step. Report instances when
;   required data were not found.

  if n_elements(root_dir) eq 0 then begin
     print, 'red_check_calibrations : Provide a root directory.' 
     return
  endif

  if n_elements(logfile) eq 0 then begin
     tstamp = red_timestamp(/utc)
     logfile = 'check_calibrations_'+tstamp+'.txt'
  endif
  openw, llun, logfile, /get_lun

  ;; Increment one of these when finding a problem.
  Ndarkproblems = 0             
  Nflatproblems = 0            
  Npolcalproblems = 0          
  Npinholeproblems = 0         
  Nscienceproblems = 0   

  printf, llun
  printf, llun, 'Examine '+root_dir

  if keyword_set(all) then begin
     
     darks = 1
     flats = 1
     polcal = 1
     pinholes = 1

  endif
  
  ;; Make lists of calibrations directories
  darkdirs    = file_search(root_dir+'/dark*/*',   count = Ndarkdirs,    /fold)
  flatdirs    = file_search(root_dir+'/flat*/*',   count = Nflatdirs,    /fold)
  pinholedirs = file_search(root_dir+'/pinh*/*',   count = Npinholedirs, /fold)
  polcaldirs  = file_search(root_dir+'/polc*/*',   count = Npolcaldirs,  /fold)
  pfscandirs  = file_search(root_dir+'/pfscan*/*', count = Npfscandirs,  /fold)

  ;; We do not check for prefilter scan data as they are not required
  ;; calibrations. We also do not include them in the science data
  ;; directories, so there is no check if any calibration data needed
  ;; for the pfscan data are there.
  
  ;; The rest must be science directories
  nonsciencedirs = [darkdirs, flatdirs, pinholedirs, polcaldirs, pfscandirs]
  dirs = file_search(root_dir+'/*/*', count = Ndirs)
  for idir = 0, Ndirs-1 do begin
     if total(dirs[idir] eq nonsciencedirs) eq 0 then begin
        if n_elements(sciencedirs) eq 0 then begin
           sciencedirs = dirs[idir]
        endif else begin
           sciencedirs = [sciencedirs, dirs[idir]]
        endelse
     endif
  endfor                        ; idir
  Nsciencedirs = n_elements(sciencedirs)

  
  if keyword_set(darks) then begin
     ;; See what darks data there are.

     printf, llun
     printf, llun, 'Dark frame data.'
     printf, llun

     print, 'Look for darks.'

     ;; List all directories with actual dark frames in them. Assumes
     ;; directory structure similar to this:
     ;; root_dir+'darks/<timestamp>/<camera name>/'.
     for idir = 0, Ndarkdirs-1 do begin
        subdirs = file_search(darkdirs[idir]+'/*', count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
           if n_elements(darksubdirs) eq 0 then begin
              darksubdirs = subdirs
           endif else begin
              darksubdirs = [darksubdirs, subdirs]
           endelse
        endif
     endfor                     ; idir

     ;; Go through the list, find out what cameras and how many frames
     Ndarksubdirs = n_elements(darksubdirs)
     if Ndarksubdirs gt 0 then begin
        Ndarkframes = lonarr(Ndarksubdirs) 
        darkcams = strarr(Ndarksubdirs)
        for idir = 0, Ndarksubdirs-1 do begin
           fnames = file_search(darksubdirs[idir]+'/cam*')
           red_extractstates, file_basename(fnames), cam = cam
           ucam = cam[uniq(cam,sort(cam))]
           if n_elements(ucam) gt 1 then begin
              printf, llun, 'Warning: More than one camera in '+darksubdirs[idir]
           endif else begin
              darkcams[idir] = ucam
              Ndarkframes[idir] = n_elements(fnames)
           endelse
           printf, llun, string(darkcams[idir], '(a12)') $
                   + ' :'+string(Ndarkframes[idir], format = '(i6)') $
                   + ' ' + darksubdirs[idir]
        endfor                  ; idir
     endif else begin           ; Ndarksubdirs
        printf, llun, 'No darks data'
        darkcams = []
     endelse                    ; Ndarksubdirs
  endif                         ; darks
  
  
  if keyword_set(flats) then begin
     ;; See what flats data there is.

     printf, llun
     printf, llun, 'Flat field data.'
     printf, llun

     print, 'Look for flats.'

     ;; List all directories with actual flat frames in them. Assumes
     ;; directory structure similar to this:
     ;; root_dir+'flats/<timestamp>/<camera name>/'. 
     for idir = 0, Nflatdirs-1 do begin
        subdirs = file_search(flatdirs[idir]+'/*', count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
           if n_elements(flatsubdirs) eq 0 then begin
              flatsubdirs = subdirs
           endif else begin
              flatsubdirs = [flatsubdirs, subdirs]
           endelse
        endif
     endfor                     ; idir

     ;; Go through the list, find out what cameras, what states, and
     ;; how many frames.
     Nflatsubdirs = n_elements(flatsubdirs)
     if Nflatsubdirs gt 0 then begin
        for idir = 0, Nflatsubdirs-1 do begin
           fnames = file_search(flatsubdirs[idir]+'/cam*', count = Nflats)

           if strmatch(file_basename(flatsubdirs[idir]),'Crisp-?') then begin
              ;; CRISP data

              ;; Are there lcd files that should be dropped?
              if Nflats gt 0 then red_match_strings, fnames, '*.lcd.*' $
                                                     , nonmatching_strings = fnames $
                                                     , Nnonmatching = Nflats
              
              red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate
              ucam = cam[uniq(cam,sort(cam))]
              if strmatch(file_basename(flatsubdirs[idir]),'Crisp-W') then begin
                 ;; For CRISP WB data we only care about the prefilter
                 for ii = 0, n_elements(fullstate)-1 do $
                    fullstate[ii] = (strsplit(fullstate[ii],'.',/extract))[0] 
              endif             ; Crisp-W
              ufullstate = fullstate[uniq(fullstate,sort(fullstate))]
              camstates = ucam+'.'+ufullstate
           endif else begin
              ;; Blue data.
              red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate, /blue
              ucam = cam[uniq(cam,sort(cam))]
              ufullstate = fullstate[uniq(fullstate,sort(fullstate))]
              if ufullstate ne '' then camstates = ucam+'.'+ufullstate else camstates = ucam
              ;; We need the exposure time as part of the state
              h = fzhead(fnames[0])
              dT = strtrim(round((double(strmid(h, strpos(h, 'Te=')+20, 9)) $
                                  - double(strmid(h, strpos(h, 'Ts=')+20, 9)))*1000), 2) + 'ms'
              camstates += '.' + dT
           endelse

           if n_elements(ucam) gt 1 then begin
              printf, llun, 'Warning: More than one camera in '+flatsubdirs[idir]
           endif else begin
              if n_elements(flatcamstates) eq 0 then begin
                 flatcamstates = camstates
                 Nflatframes = replicate(0L, n_elements(camstates))
              endif else begin
                 flatcamstates = [camstates, flatcamstates]
                 Nflatframes = [replicate(0L, n_elements(camstates)), Nflatframes]
              endelse
              ;; Find out how many frames of each kind:
              for i = 0, n_elements(camstates)-1 do begin
                 Nflatframes[i] = n_elements(where(strmatch(fnames,'*'+fullstate[i]+'*')))
                 printf, llun, string(flatcamstates[i], '(a30)') $
                         + ' :'+string(Nflatframes[i], format = '(i6)') $
                         + ' ' + flatsubdirs[idir]
              endfor
           endelse
           
        endfor                  ; idir
        
     endif else begin           ; Nflatsubdirs
        flatcamstates = []
        printf, llun, 'No flats data'
     endelse                    ; Nflatsubdirs
  endif                         ; flats
  
  
  if keyword_set(polcal) then begin
     ;; See what polcal data there is.

     printf, llun
     printf, llun, 'Polcal data.'
     printf, llun

     print, 'Look for polcal data.'

     ;; List all directories with actual polcal data in them. Assumes
     ;; directory structure similar to this:
     ;; root_dir+'polcal/<prefilter>/<timestamp>/Crisp-?/'. 
     for idir = 0, Npolcaldirs-1 do begin
        subdirs = file_search(polcaldirs[idir]+'/*/*', count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
           if n_elements(polcalsubdirs) eq 0 then begin
              polcalsubdirs = subdirs
           endif else begin
              polcalsubdirs = [polcalsubdirs, subdirs]
           endelse
        endif
     endfor                     ; idir

     ;; Go through the list, find out what cameras and what
     ;; prefilters. Also check that the data appear complete. 
     Npolcalsubdirs = n_elements(polcalsubdirs)
     if Npolcalsubdirs gt 0 then begin
        polcalcamstates = strarr(Npolcalsubdirs)
        polcalok = replicate(1, Npolcalsubdirs)
        for idir = 0, Npolcalsubdirs-1 do begin

           fnames = file_search(polcalsubdirs[idir]+'/cam*', count = Npolcal)
           
           ;; Are there lcd files that should be dropped?
           if Npolcal gt 0 then red_match_strings, fnames, '*.lcd.*' $
                                                   , nonmatching_strings = fnames $
                                                   , Nnonmatching = Npolcal
           
           red_extractstates, file_basename(fnames), cam = cam, pref = pref, lc = lc, qw = qw, lp = lp
           ucam = cam[uniq(cam,sort(cam))]
           upref = pref[uniq(pref,sort(pref))]
           if n_elements(ucam) gt 1 then begin
              printf, llun, 'Warning: More than one camera in '+polcalsubdirs[idir]
           endif else begin

              polcalcamstates[idir] = ucam+'.'+upref

              ;; Find out if the polcal data are ok.

              ;; There should be two lp (linear polarizer) states and
              ;; the number of frames of each should be equal.
              ulp = lp[uniq(lp,sort(lp))]
              Nlp = lonarr(n_elements(ulp))
              for i = 0, n_elements(ulp)-1 do Nlp[i] = n_elements(where(strmatch(fnames,'*'+ulp[i]+'*')))
              polcalok[idir] = polcalok[idir] and Nlp[0] eq Nlp[1]

              ;; There should be four lc (liqud crystal) states and
              ;; the number of frames of each should be equal.
              ulc = lc[uniq(lc,sort(lc))]
              Nlc = lonarr(n_elements(ulc))
              for i = 0, n_elements(ulc)-1 do Nlc[i] = n_elements(where(strmatch(fnames,'*'+ulc[i]+'*')))
              polcalok[idir] = polcalok[idir] and Nlc[0] eq Nlc[1] and Nlc[1] eq Nlc[2] and Nlc[2] eq Nlc[3]

              ;; The number of qw (quarter wave plate) angles is not
              ;; fixed but the angles should be evenly distributed
              ;; over 360 degrees, repeating the 0=360 angle. And the
              ;; number of frames at each angle should be equal.
              uqw = qw[uniq(qw,sort(qw))]     
              Nqw = lonarr(n_elements(uqw))
              for i = 0, n_elements(uqw)-1 do Nqw[i] = n_elements(where(strmatch(fnames,'*'+uqw[i]+'*')))
              polcalok[idir] = polcalok[idir] and stddev(Nqw) lt 1e-1 ; Equal numbers of frames for all angles?
              polcalok[idir] = polcalok[idir] and min(uqw) eq 'qw000' ; Zero angle present?
              polcalok[idir] = polcalok[idir] and max(uqw) eq 'qw360' ; Repeated enpoint?
              dqw = deriv(long(strmid(uqw,2,3)))                      ; Differential angle
              polcalok[idir] = polcalok[idir] and stddev(dqw) lt 1e-1 ; Evenly distributed?
           endelse
           
           outline = string(polcalcamstates[idir], '(a15)') + ' : '
           if polcalok[idir] then begin
              outline += 'complete polcal data in ' 
           endif else begin
              outline += 'incomplete polcal data in ' 
           endelse
           outline += red_strreplace(polcalsubdirs[idir], root_dir, '')
           printf, llun, outline
        endfor                  ; idir
     endif else begin           ; Npocalsubdirs 
        polcalcamstates = []
        printf, llun, 'No polcal data'
     endelse                    ; Npocalsubdirs     
  endif                         ; polcal
  
  
  ;; See what pinholes data there is.
  if keyword_set(pinholes) then begin

     printf, llun
     printf, llun, 'Pinhole data.'
     printf, llun

     print, 'Look for pinholes.'

     ;; List all directories with actual pinhole data in them. Assumes
     ;; directory structure similar to this:
     ;; root_dir+'pinholes/<timestamp>/<camera name>/'. 
     for idir = 0, Npinholedirs-1 do begin
        subdirs = file_search(pinholedirs[idir]+'/*', count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
           if n_elements(pinholesubdirs) eq 0 then begin
              pinholesubdirs = subdirs
           endif else begin
              pinholesubdirs = [pinholesubdirs, subdirs]
           endelse
        endif
     endfor                     ; idir

     ;; Go through the list, find out what cameras, what states, and
     ;; how many frames.
     Npinholesubdirs = n_elements(pinholesubdirs)
     if Npinholesubdirs gt 0 then begin
        for idir = 0, Npinholesubdirs-1 do begin
           fnames = file_search(pinholesubdirs[idir]+'/cam*', count = Npinhole)

           if strmatch(file_basename(pinholesubdirs[idir]),'Crisp-?') then begin
              ;; CRISP data

              ;; Are there lcd files that should be dropped?
              if Npinhole gt 0 then red_match_strings, fnames, '*.lcd.*' $
                                                       , nonmatching_strings = fnames $
                                                       , Nnonmatching = Npinhole
              
              red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate
              ucam = cam[uniq(cam,sort(cam))]
              if strmatch(file_basename(pinholesubdirs[idir]),'Crisp-W') then begin
                 ;; For CRISP WB data we only care about the prefilter
                 for ii = 0, n_elements(fullstate)-1 do $
                    fullstate[ii] = (strsplit(fullstate[ii],'.',/extract))[0] 
              endif else begin
                 ;; For CRISP NB data we don't care about the LC state
                 for i=0,n_elements(fullstate)-1 do $
                    fullstate[i] = strjoin((strsplit(fullstate[i],'.',/extract))[0:1],'.')
              endelse
              ufullstate = fullstate[uniq(fullstate,sort(fullstate))]
              camstates = ucam+'.'+ufullstate
           endif else begin
              ;; Blue data.
              red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate, /blue
              ucam = cam[uniq(cam,sort(cam))]
              ufullstate = fullstate[uniq(fullstate,sort(fullstate))]
              if ufullstate ne '' then camstates = ucam+'.'+ufullstate else camstates = ucam
           endelse

           if n_elements(ucam) gt 1 then begin
              printf, llun, 'Warning: More than one camera in '+pinholesubdirs[idir]
           endif else begin
              if n_elements(pinholecamstates) eq 0 then begin
                 pinholecamstates = camstates
                 Npinholeframes = replicate(0L, n_elements(camstates))
              endif else begin
                 pinholecamstates = [camstates, pinholecamstates]
                 Npinholeframes = [replicate(0L, n_elements(camstates)), Npinholeframes]
              endelse
              ;; Find out how many frames of each kind:
              for i = 0, n_elements(camstates)-1 do begin
                 Npinholeframes[i] = n_elements(where(strmatch(fnames,'*'+ufullstate[i]+'*')))
                 printf, llun, string(pinholecamstates[i], '(a30)') $
                         + ' :'+string(Npinholeframes[i], format = '(i5)') $
                         + ' ' + pinholesubdirs[idir]
              endfor
           endelse
           
        endfor                  ; idir
     endif else begin           ; Npinholesubdirs
        pinholecamstates = []
        printf, llun, 'No pinhole data.'
     endelse                    ; Npinholesubdirs
  endif                         ; pinholes
  

  
  ;; Now examine the science data directories, find out cameras and
  ;; states and exposure times and see if we have matching
  ;; calibrations.
  printf, llun
  printf, llun, 'Check the following science directories: '
  printf, llun, '* '+red_strreplace(sciencedirs, root_dir, ''), format = '(a0)'
  printf, llun

  print
  print, 'Check the following science directories: '
  print, '* '+red_strreplace(sciencedirs, root_dir, ''), format = '(a0)'
  print

  for idir = 0, Nsciencedirs-1 do begin

     sciencesubdirs = file_search(sciencedirs[idir]+'/*', count = Nsubdirs)
     for isubdir = 0, Nsubdirs-1 do begin

        print, red_strreplace(sciencesubdirs[isubdir], root_dir, '')
        printf, llun, red_strreplace(sciencesubdirs[isubdir], root_dir, '')

        snames = file_search(sciencesubdirs[isubdir]+'/cam*', count = Nscienceframes)

        ;; Are there lcd files that should be dropped?
        if Nscienceframes gt 0 then red_match_strings, snames, '*.lcd.*' $
           , nonmatching_strings = snames, Nnonmatching = Nscienceframes

        if Nscienceframes gt 0 then begin
           
           if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-?') then begin
              ;; CRISP data
              red_extractstates, file_basename(snames), cam = cam, fullstate = fullstate, pref = pref, lc = lc
              upref = pref[uniq(pref,sort(pref))]
              ucam = cam[uniq(cam,sort(cam))]
              if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-W') then begin
                 ;; For CRISP WB data we only care about the prefilter
                 for ii = 0, n_elements(fullstate)-1 do $
                    fullstate[ii] = (strsplit(fullstate[ii],'.',/extract))[0]
              endif             ; Crisp-W
              ufullstate = fullstate[uniq(fullstate,sort(fullstate))]
              sciencecamstates = ucam+'.'+ufullstate
           endif else begin
              ;; Blue data.
              red_extractstates, file_basename(snames), cam = cam, fullstate = fullstate, /blue
              ucam = cam[uniq(cam,sort(cam))]
              ufullstate = fullstate[uniq(fullstate,sort(fullstate))]
              if ufullstate ne '' then sciencecamstates = ucam+'.'+ufullstate else sciencecamstates = ucam
              ;; We need the exposure time as part of the state
              h = fzhead(snames[0])
              dT = strtrim(round((double(strmid(h, strpos(h, 'Te=')+20, 9)) $
                                  - double(strmid(h, strpos(h, 'Ts=')+20, 9)))*1000), 2) + 'ms'
              sciencecamstates += '.' + dT
           endelse

           if n_elements(ucam) gt 1 then begin
              printf, llun, '   Warning: More than one camera in '+sciencesubdirs[isubdir]
              print, '      More than one camera in '+sciencesubdirs[isubdir]
              Nscienceproblems += 1
              break
           endif


           ;; Do we have darks for these data?
           if keyword_set(darks) then begin

              print, '   Check darks.'

              printf, llun, '   There are ' + strtrim(Nscienceframes, 2) + ' ' +ucam + ' science frames.'
              idark = (where(ucam eq darkcams, Nhits))[0]
              if Nhits gt 0 then begin
                 printf, llun, '   There are '+strtrim(Ndarkframes[idark], 2)+' dark frames for '+ucam+'.'
              endif else begin
                 printf, llun, '   No darks for '+ucam+'!   ---------- Warning!'
                 print, '      No darks for '+ucam+'!'
                 Ndarkproblems += 1
              endelse
           endif                ; darks
           


           ;; Do we have flats for these data?
           if keyword_set(flats) then begin

              print, '   Check flats.'
              
              ;; Find all camstates, loop through them and compare with the flats camstates.
;              sciencecamstates = cam+'.'+fullstate

              for i = 0, n_elements(ufullstate)-1 do begin
                 Nstateframes = n_elements(where(strmatch(snames,'*'+ufullstate[i]+'*')))
                 iflat = (where(sciencecamstates[i] eq flatcamstates, Nhits))[0]
                 outline = '   ' + sciencecamstates[i] + ' with ' + strtrim(Nstateframes, 2) + ' science frames : '
                 if Nhits gt 0 then begin
                    outline += strtrim(Nflatframes[iflat], 2) + ' flat frames.'
                 endif else begin
                    outline += ' No flat frames!   ---------- Warning!'
                    print, '      No flats for '+sciencecamstates[i]
                    Nflatproblems += 1
                 endelse       ; Nhits
                 printf, llun, outline
              endfor            ; i

           endif                ; flats

           
           ;; Do we have pinholes data?
           if keyword_set(pinholes) then begin

              print, '   Check pinholes.'

              for i = 0, n_elements(ufullstate)-1 do begin
                 Nstateframes = n_elements(where(strmatch(snames,'*'+ufullstate[i]+'*')))
                 if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-[TR]') then begin ; CRISP NB
                    ;; We don't care about the LC state for NB pinhole data:
                    sciencecamstate = strjoin((strsplit(sciencecamstates[i],'.',/extract))[0:2],'.')
                    ipinhole = (where(sciencecamstate eq pinholecamstates, Nhits))[0]
                    outline = '   ' + sciencecamstates[i] + ' with ' + strtrim(Nstateframes, 2) + ' science frames : '
                    if Nhits gt 0 then begin
                       ;; Pinhole data exists for this state
                       outline += strtrim(Npinholeframes[ipinhole], 2) + ' pinhole frames.'
                    endif else begin
                       ;; We don't always have pinholes for all
                       ;; wavelength points. Report nearest wavelength
                       ;; in that case.
                       
                       ;; First make sure we have data with the same
                       ;; prefilter:
                       if n_elements(pinholecamstates) gt 0 then begin
                          indx = where(strmatch(pinholecamstates, $
                                                strjoin((strsplit(sciencecamstate,'.',/extract))[0:1],'.')+'*'), Nmatch)
                       endif else Nmatch = 0
                       if Nmatch eq 0 then begin
                          outline += ' No pinhole frames!   ---------- Warning!'
                          print, '      No pinhole data for '+sciencecamstate
                          Npinholeproblems += 1
                       endif else begin
                          ;; These are the existing pinhole cam states
                          ;; from the current prefilter:
                          ppinholecamstates = pinholecamstates[indx]
                          sciencetuning = (strsplit(sciencecamstate,'.',/extract))[2]
                          sciencelambda = total(double(strsplit(sciencetuning,'_',/extract)) * [1d,1e-3])
                          pinholelambda = dblarr(Nmatch)
                          for ii = 0, Nmatch-1 do begin
                             pinholetuning = (strsplit(ppinholecamstates[ii],'.',/extract))[2]
                             pinholelambda[ii] = total(double(strsplit(pinholetuning,'_',/extract)) * [1d,1e-3])
                          endfor ; ii
                          dlambda = round(min(abs(sciencelambda - pinholelambda),mloc) * 1e3)
                          ipinhole = indx[mloc]
                          outline += strtrim(Npinholeframes[ipinhole], 2) + ' pinhole frames ' $
                                     + strtrim(dlambda, 2)+' mÅ away.'
                       endelse
                    endelse        ; Nhits
                 endif else begin  ; CRISP WB and blue
                    if ~strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-W') then begin
                       ;; For blue data we need to remove the exposure
                       ;; time (but keep the tilt filter tuning).
                       sciencecamstate = strjoin((strsplit(sciencecamstates[i], '.', /extract, count = Nsplit) $
                                                 )[0:Nsplit-2], '.')
                    endif else begin
                       sciencecamstate = sciencecamstates[i]
                    endelse
                    ipinhole = (where(sciencecamstate eq pinholecamstates, Nhits))[0]
                    outline = '   ' + sciencecamstate + ' with ' + strtrim(Nstateframes, 2) $
                              + ' science frames : '
                    if Nhits gt 0 then begin
                       outline += strtrim(Npinholeframes[ipinhole], 2) + ' pinhole frames.'
                    endif else begin
                       outline += ' No pinhole frames!   ---------- Warning!'
                       print, '      No pinhole data for '+sciencecamstate
                       Npinholeproblems += 1
                    endelse     ; Nhits
                 endelse
                 printf, llun, outline
              endfor            ; i

           endif                ; pinholes


           ;; Do we need polcal for these data?
           if keyword_set(polcal) then begin
              if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-[TR]') then begin
                 ;; Only for CRISP NB

                 Npref = n_elements(upref)
                 for ipref = 0, Npref-1 do begin

                    this_pref = upref[ipref]
                    indx = where(strmatch(fullstate, this_pref+'\.*\.*'), Nmatch)
                    
                    if Nmatch le 0 then begin
                       print, 'This should not happen!'
                       stop
                    endif

                    these_lc = lc[indx]
                    ulc = these_lc[uniq(these_lc,sort(these_lc))]


                    if n_elements(ulc) gt 1 then begin
                       print, '   Check polcal for '+this_pref+'.'
                       if n_elements(ulc) gt 4 then begin
                          help, ulc
                          stop
                       endif
                       ipolcal = where(ucam+'.'+this_pref eq polcalcamstates, Nhits)
                       if Nhits eq 0 then begin
                          printf, llun, '   no polcal data.   ---------- Warning!'
                          print, '      No polcal data'
                          Npolcalproblems += 1
                       endif else begin
                          for i = 0, Nhits-1 do begin
                             ;; Should print out also the polcal
                             ;; directory here.
                             if polcalok[ipolcal[i]] then begin
                                printf, llun, '   complete polcal data in ' $
                                        + red_strreplace(polcalsubdirs[ipolcal[i]], root_dir, '')
                             endif else begin
                                printf, llun, '   incomplete polcal data in ' $
                                        + red_strreplace(polcalsubdirs[ipolcal[i]], root_dir, '') $
                                        + '   ---------- Warning!'
                                print, '      Incomplete polcal data'
                                Npolcalproblems += 1
                             endelse
                          endfor ; i
                       endelse 
                    endif       ; lc
                 endfor         ; ipref
              endif             ; CRISP data
           endif                ; polcal

        endif                   ; Nscienceframes
     endfor                     ; isubdir   
  endfor                        ; idir
   
  print
  print, '**************************************************'
  if Ndarkproblems + Nflatproblems + Npolcalproblems + Npinholeproblems + Nscienceproblems eq 0 then begin
     print, 'red_check_calibrations : No problems detected.'
  endif else begin
     if Ndarkproblems gt 0 then $
        print, 'red_check_calibrations : '+strtrim(Ndarkproblems, 2)+' dark problems found.'
     if Nflatproblems gt 0 then $
        print, 'red_check_calibrations : '+strtrim(Nflatproblems, 2)+' flat problems found.'
     if Npinholeproblems gt 0 then $
        print, 'red_check_calibrations : '+strtrim(Npinholeproblems, 2)+' pinhole problems found.'
     if Npolcalproblems gt 0 then $
        print, 'red_check_calibrations : '+strtrim(Npolcalproblems, 2)+' polcal problems found.'
     if Nscienceproblems gt 0 then $
        print, 'red_check_calibrations : '+strtrim(Nscienceproblems, 2)+' science problems found.'
  endelse

  print, 'Detailed output in '+logfile+'.'
  print, '**************************************************'
 
  printf, llun
  printf, llun, '**************************************************'
  if Ndarkproblems + Nflatproblems + Npolcalproblems + Npinholeproblems + Nscienceproblems eq 0 then begin
     printf, llun, 'red_check_calibrations : No problems detected.'
  endif else begin
     if Ndarkproblems gt 0 then $
        printf, llun, 'red_check_calibrations : '+strtrim(Ndarkproblems, 2)+' dark problems found.'
     if Nflatproblems gt 0 then $
        printf, llun, 'red_check_calibrations : '+strtrim(Nflatproblems, 2)+' flat problems found.'
     if Npinholeproblems gt 0 then $
        printf, llun, 'red_check_calibrations : '+strtrim(Npinholeproblems, 2)+' pinhole problems found.'
     if Npolcalproblems gt 0 then $
        printf, llun, 'red_check_calibrations : '+strtrim(Npolcalproblems, 2)+' polcal problems found.'
     if Nscienceproblems gt 0 then $
        printf, llun, 'red_check_calibrations : '+strtrim(Nscienceproblems, 2)+' science problems found.'
  endelse

  printf, llun, 'Detailed output in '+logfile+'.'
  printf, llun, '**************************************************'
  
 

  free_lun, llun

end
