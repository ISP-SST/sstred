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
; 
;-
pro red_check_calibrations, root_dir $
                            , all = all $
                            , darks = darks $
                            , flats = flats $
                            , polcal = polcal $
                            , pinholes = pinholes $
                            , logfile = logfile

  if n_elements(root_dir) eq 0 then begin
     print, 'red_check_calibrations : Provide a root directory.' 
     return
  endif

  if n_elements(logfile) eq 0 then logfile = 'check_calibrations_output.txt'
  openw, llun, logfile, /get_lun

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

  ;; The rest must be science data
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

  
  ;; See what darks data there is.
  if keyword_set(darks) then begin

     printf, llun
     printf, llun, 'Dark frame data.'
     printf, llun

     print, 'Look for darks.'

     ;; List all directories with actual dark frames in them
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
     Ndarkframes = lonarr(Ndarksubdirs) 
     darkcams = strarr(Ndarksubdirs)
     for idir = 0, Ndarksubdirs-1 do begin
        fnames = file_search(darksubdirs[idir]+'/cam*')
        red_extractstates, file_basename(fnames), cam = cam
        cam = cam[uniq(cam,sort(cam))]
        if n_elements(cam) gt 1 then begin
           printf, llun, 'Warning: More than one camera in '+darksubdirs[idir]
        endif else begin
           darkcams[idir] = cam
           Ndarkframes[idir] = n_elements(fnames)
        endelse
        printf, llun, string(darkcams[idir], '(a12)') $
               + ' :'+string(Ndarkframes[idir], format = '(i6)') $
               + ' ' + darksubdirs[idir]
     endfor                     ; idir
     
  endif                         ; darks
  
  
  ;; See what flats data there is.
  if keyword_set(flats) then begin

     printf, llun
     printf, llun, 'Flat field data.'
     printf, llun

     print, 'Look for flats.'

     ;; List all directories with actual flat frames in them.
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
     for idir = 0, Nflatsubdirs-1 do begin
        fnames = file_search(flatsubdirs[idir]+'/cam*')

        if strmatch(file_basename(flatsubdirs[idir]),'Crisp-?') then begin
           ;; CRISP data
           red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate
           cam = cam[uniq(cam,sort(cam))]
           if strmatch(file_basename(flatsubdirs[idir]),'Crisp-W') then begin
              ;; For CRISP WB data we only care about the prefilter
              for ii = 0, n_elements(fullstate)-1 do $
                 fullstate[ii] = (strsplit(fullstate[ii],'.',/extract))[0] 
           endif                ; Crisp-W
           fullstate = fullstate[uniq(fullstate,sort(fullstate))]
           camstates = cam+'.'+fullstate
        endif else begin
           ;; Blue data.
           red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate, /blue
           cam = cam[uniq(cam,sort(cam))]
           fullstate = fullstate[uniq(fullstate,sort(fullstate))]
           if fullstate ne '' then camstates = cam+'.'+fullstate else camstates = cam
           ;; We need the exposure time as part of the state
           h = fzhead(fnames[0])
           dT = strtrim(round((double(strmid(h, strpos(h, 'Te=')+20, 9)) $
                               - double(strmid(h, strpos(h, 'Ts=')+20, 9)))*1000), 2) + 'ms'
           camstates += '.' + dT
        endelse

        if n_elements(cam) gt 1 then begin
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
        
     endfor                     ; idir
    

  endif                         ; flats
  
  
  ;; See what polcal data there is.
  if keyword_set(polcal) then begin

     printf, llun
     printf, llun, 'Polcal data.'
     printf, llun

     print, 'Look for polcal data.'

     ;; List all directories with actual polcal data in them.
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
     polcalcamstates = strarr(Npolcalsubdirs)
     polcalok = replicate(1, Npolcalsubdirs)
     for idir = 0, Npolcalsubdirs-1 do begin

        fnames = file_search(polcalsubdirs[idir]+'/cam*')
        
        red_extractstates, file_basename(fnames), cam = cam, pref = pref, lc = lc, qw = qw, lp = lp
        cam = cam[uniq(cam,sort(cam))]
        pref = pref[uniq(pref,sort(pref))]
        if n_elements(cam) gt 1 then begin
           printf, llun, 'Warning: More than one camera in '+polcalsubdirs[idir]
        endif else begin
           polcalcamstates[idir] = cam+'.'+pref
           ;; Find out if the polcal data are ok.
           ;; There should be two lp (linear polarizer) states and the
           ;; number of frames of each should be equal.
           lp = lp[uniq(lp,sort(lp))]
           Nlp = lonarr(n_elements(lp))
           for i = 0, n_elements(lp)-1 do Nlp[i] = n_elements(where(strmatch(fnames,'*'+lp[i]+'*')))
           polcalok[idir] = polcalok[idir] and Nlp[0] eq Nlp[1]
           ;; There should be four lc (liqud crystal) states and the
           ;; number of frames of each should be equal.
           lc = lc[uniq(lc,sort(lc))]
           Nlc = lonarr(n_elements(lc))
           for i = 0, n_elements(lc)-1 do Nlc[i] = n_elements(where(strmatch(fnames,'*'+lc[i]+'*')))
           polcalok[idir] = polcalok[idir] and Nlc[0] eq Nlc[1] and Nlc[1] eq Nlc[2] and Nlc[2] eq Nlc[3]
           ;; The number of qw (quarter wave plate) angles is not
           ;; fixed but the angles should be evenly distributed over
           ;; 360 degrees, repeating the 0=360 angle. And the number
           ;; of frames at each angle should be equal.
           qw = qw[uniq(qw,sort(qw))]     
           Nqw = lonarr(n_elements(qw))
           for i = 0, n_elements(qw)-1 do Nqw[i] = n_elements(where(strmatch(fnames,'*'+qw[i]+'*')))
           polcalok[idir] = polcalok[idir] and stddev(Nqw) lt 1e-1 ; Equal numbers of frames for all angles?
           polcalok[idir] = polcalok[idir] and min(qw) eq 'qw000'  ; Zero angle present?
           polcalok[idir] = polcalok[idir] and max(qw) eq 'qw360'  ; Repeated enpoint?
           dqw = deriv(long(strmid(qw,2,3)))                 ; Differential angle
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
     endfor                     ; idir
     
  endif                         ; polcal
  
  
  ;; See what pinholes data there is.
  if keyword_set(pinholes) then begin

     printf, llun
     printf, llun, 'Pinhole data.'
     printf, llun

     print, 'Look for pinholes.'

     ;; List all directories with actual pinhole data in them.
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
     for idir = 0, Npinholesubdirs-1 do begin
        fnames = file_search(pinholesubdirs[idir]+'/cam*')

        if strmatch(file_basename(pinholesubdirs[idir]),'Crisp-?') then begin
           ;; CRISP data
           red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate
           cam = cam[uniq(cam,sort(cam))]
           if strmatch(file_basename(pinholesubdirs[idir]),'Crisp-W') then begin
              ;; For CRISP WB data we only care about the prefilter
              for ii = 0, n_elements(fullstate)-1 do $
                 fullstate[ii] = (strsplit(fullstate[ii],'.',/extract))[0] 
           endif else begin
              ;; For CRISP NB data we don't care about the LC state
              for i=0,n_elements(fullstate)-1 do $
                 fullstate[i] = strjoin((strsplit(fullstate[i],'.',/extract))[0:1],'.')
           endelse
           fullstate = fullstate[uniq(fullstate,sort(fullstate))]
           camstates = cam+'.'+fullstate
        endif else begin
           ;; Blue data.
           red_extractstates, file_basename(fnames), cam = cam, fullstate = fullstate, /blue
           cam = cam[uniq(cam,sort(cam))]
           fullstate = fullstate[uniq(fullstate,sort(fullstate))]
           if fullstate ne '' then camstates = cam+'.'+fullstate else camstates = cam
        endelse

        if n_elements(cam) gt 1 then begin
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
              Npinholeframes[i] = n_elements(where(strmatch(fnames,'*'+fullstate[i]+'*')))
              printf, llun, string(pinholecamstates[i], '(a30)') $
                      + ' :'+string(Npinholeframes[i], format = '(i5)') $
                      + ' ' + pinholesubdirs[idir]
           endfor
        endelse
        
     endfor                     ; idir

  endif                         ; pinholes
  

  
  ;; Now examine the data directories, find out cameras and states
  ;; and exposure times and see if we have matching calibrations.
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
 
        if Nscienceframes gt 0 then begin
           
;           red_extractstates, file_basename(snames), cam = cam, fullstate = fullstate
;           cam = cam[uniq(cam,sort(cam))]
;           fullstate = fullstate[uniq(fullstate,sort(fullstate))]
           
           if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-?') then begin
              ;; CRISP data
              red_extractstates, file_basename(snames), cam = cam, fullstate = fullstate, pref = pref, lc = lc
              cam = cam[uniq(cam,sort(cam))]
              lc = lc[uniq(lc,sort(lc))]
              if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-W') then begin
                 ;; For CRISP WB data we only care about the prefilter
                 for ii = 0, n_elements(fullstate)-1 do $
                    fullstate[ii] = (strsplit(fullstate[ii],'.',/extract))[0]
              endif             ; Crisp-W
              fullstate = fullstate[uniq(fullstate,sort(fullstate))]
              sciencecamstates = cam+'.'+fullstate
           endif else begin
              ;; Blue data.
              red_extractstates, file_basename(snames), cam = cam, fullstate = fullstate, /blue
              cam = cam[uniq(cam,sort(cam))]
              fullstate = fullstate[uniq(fullstate,sort(fullstate))]
              if fullstate ne '' then sciencecamstates = cam+'.'+fullstate else sciencecamstates = cam
              ;; We need the exposure time as part of the state
              h = fzhead(snames[0])
              dT = strtrim(round((double(strmid(h, strpos(h, 'Te=')+20, 9)) $
                                  - double(strmid(h, strpos(h, 'Ts=')+20, 9)))*1000), 2) + 'ms'
              sciencecamstates += '.' + dT
           endelse

           if n_elements(cam) gt 1 then begin
              printf, llun, '   Warning: More than one camera in '+sciencesubdirs[isubdir]
              break
           endif


           ;; Do we have darks for these data?
           if keyword_set(darks) then begin

              print, '   Check darks.'

              printf, llun, '   There are ' + strtrim(Nscienceframes, 2) + ' ' +cam + ' science frames.'
              idark = (where(cam eq darkcams, Nhits))[0]
              if Nhits gt 0 then begin
                 printf, llun, '   There are '+strtrim(Ndarkframes[idark], 2)+' dark frames for '+cam+'.'
              endif else begin
                 printf, llun, '   Warning: No darks for '+cam+'!'
              endelse
           endif                ; darks
           


           ;; Do we have flats for these data?
           if keyword_set(flats) then begin

              print, '   Check flats.'
              
              ;; Find all camstates, loop through them and compare with the flats camstates.
;              sciencecamstates = cam+'.'+fullstate

              for i = 0, n_elements(fullstate)-1 do begin
                 Nstateframes = n_elements(where(strmatch(snames,'*'+fullstate[i]+'*')))
                 iflat = (where(sciencecamstates[i] eq flatcamstates, Nhits))[0]
                 outline = '   ' + sciencecamstates[i] + ' with ' + strtrim(Nstateframes, 2) + ' science frames : '
                 if Nhits gt 0 then begin
                    outline += strtrim(Nflatframes[iflat], 2) + ' flat frames.'
                 endif else begin
                    outline += 'No flat frames!'
                 endelse       ; Nhits
                 printf, llun, outline
              endfor            ; i

           endif                ; flats

           ;; Do we need polcal for these data?
           if keyword_set(polcal) then begin
              if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-[TR]') then begin
                 if n_elements(lc) gt 1 then begin
                    print, '   Check polcal.'
                    ipolcal = where(cam+'.'+pref eq polcalcamstates, Nhits)
                    if Nhits eq 0 then begin
                       printf, llun, '  no polcal data.'
                    endif else begin
                       for i = 0, Nhits-1 do begin
                          ;; Should print out also the polcal
                          ;; directory here.
                          if polcalok[ipolcal[i]] then begin
                             printf, llun, '   complete polcal data in ' $
                                     + red_strreplace(polcalsubdirs[ipolcal[i]], root_dir, '')
                          endif else begin
                             printf, llun, '   incomplete polcal data in ' $
                                     + red_strreplace(polcalsubdirs[ipolcal[i]], root_dir, '')
                          endelse
                       endfor   ; i
                    endelse 
                 endif          ; lc
              endif             ; CRISP data
           endif                ; polcal
           
           
           ;; Do we have pinholes data?
           if keyword_set(pinholes) then begin

              print, '   Check pinholes.'

              for i = 0, n_elements(fullstate)-1 do begin
                 Nstateframes = n_elements(where(strmatch(snames,'*'+fullstate[i]+'*')))
                 if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-[TR]') then begin ; CRISP NB
                    ;; We don't care about the LC state for NB pinhole data:
                    sciencecamstate = strjoin((strsplit(sciencecamstates[i],'.',/extract))[0:2],'.')
                    ipinhole = (where(sciencecamstate eq pinholecamstates, Nhits))[0]
                    outline = '   ' + sciencecamstates[i] + ' with ' + strtrim(Nstateframes, 2) + ' science frames : '
                    if Nhits gt 0 then begin
                       ;; Pinhole data exists for this state
                       outline += strtrim(Npinholeframes[ipinhole], 2) + ' pinhole frames.'
                    endif else begin
                       ;; We don't always have pinholes for all wavelength
                       ;; points. Report nearest wavelength in that case.
                       
                       ;; First make sure we have data with the same prefilter:
                       indx = where(strmatch(pinholecamstates, $
                                             strjoin((strsplit(sciencecamstate,'.',/extract))[0:1],'.')+'*'), Nmatch)
                       if Nmatch eq 0 then begin
                          outline += 'No pinhole frames!'
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
                    outline = '   ' + sciencecamstate + ' with ' + strtrim(Nstateframes, 2) + ' science frames : '
                    if Nhits gt 0 then begin
                       outline += strtrim(Npinholeframes[ipinhole], 2) + ' pinhole frames.'
                    endif else begin
                       stop
                       outline += 'No pinhole frames!'
                    endelse     ; Nhits
                 endelse
                 printf, llun, outline
              endfor            ; i

           endif                ; pinholes

        endif                   ; Nscienceframes
     endfor                     ; isubdir   
  endfor                        ; idir
  
  free_lun, llun
  
  print
  print, 'check_calibrations : Output in '+logfile+'.'

end
