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
; :Keywords:
; 
;
; :History:
;
;    2015-09-01 : MGL. Looking for different kinds of data ready.
;                 Checking science data for corresponding dark frame
;                 data ready.
;
;    2015-09-01 : MGL. Checking science data for corresponding flat
;                 field data ready.
; 
;-
pro red::check_calibrations, all = all $
                             , darks = darks $
                             , flats = flats $
                             , polcal = polcal $
                             , pinholes = pinholes $
                             , logfile = logfile

  if n_elements(logfile) eq 0 then logfile = 'check_calibrations_output.txt'
  openw, llun, logfile, /get_lun

  if keyword_set(all) then begin
     
     darks = 1
     flats = 1
     polcal = 1
     pinholes = 1

  endif
  
  ;; Make lists of calibrations directories
  darkdirs   = file_search(self.root_dir+'/dark*/*',   count = Ndarkdirs,    /fold)
  flatdirs   = file_search(self.root_dir+'/flat*/*',   count = Nflatdirs,    /fold)
  pinhdirs   = file_search(self.root_dir+'/pinh*/*',   count = Npinholedirs, /fold)
  polcaldirs = file_search(self.root_dir+'/polc*/*',   count = Npolcaldirs,  /fold)
  pfscandirs = file_search(self.root_dir+'/pfscan*/*', count = Npfscandirs,  /fold)

  ;; The rest must be science data
  nonsciencedirs = [darkdirs, flatdirs, pinhdirs, polcaldirs, pfscandirs]
  dirs = file_search(self.root_dir+'/*/*', count = Ndirs)
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
        dnames = file_search(darksubdirs[idir]+'/cam*')
        red_extractstates, file_basename(dnames), cam = cam
        cam = cam[uniq(cam,sort(cam))]
        if n_elements(cam) gt 1 then begin
           printf, llun, 'Warning: More than one camera in '+darksubdirs[idir]
        endif else begin
           darkcams[idir] = cam
           Ndarkframes[idir] = n_elements(dnames)
        endelse
        printf, llun, string(darkcams[idir], '(a12)') $
               + ' :'+string(Ndarkframes[idir], format = '(i6)') $
               + ' ' + darksubdirs[idir]
     endfor                     ; idir
     
  endif                         ; darks
  
  
  ;; See what flats data there is.
  if keyword_set(flats) then begin
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
     print, 'Look for polcal (not implemented yet).'
  endif                         ; polcal
  
  
  ;; See what pinholes data there is.
  if keyword_set(pinholes) then begin
     print, 'Look for pinholes (not implemented yet).'
  endif                         ; pinholes
  

  
  ;; Now examine the data directories, find out cameras and states
  ;; and exposure times and see if we have matching calibrations.
  printf, llun
  printf, llun, 'Examine the following science directories: ', sciencedirs
  printf, llun

  print
  print, 'Examine the following science directories: ', sciencedirs
  print

  for idir = 0, Nsciencedirs-1 do begin

     sciencesubdirs = file_search(sciencedirs[idir]+'/*', count = Nsubdirs)
     for isubdir = 0, Nsubdirs-1 do begin

        print, sciencesubdirs[isubdir]
        printf, llun, sciencesubdirs[isubdir]

        snames = file_search(sciencesubdirs[isubdir]+'/cam*', count = Nscienceframes)
 
        if Nscienceframes gt 0 then begin
           
;           red_extractstates, file_basename(snames), cam = cam, fullstate = fullstate
;           cam = cam[uniq(cam,sort(cam))]
;           fullstate = fullstate[uniq(fullstate,sort(fullstate))]
           
           if strmatch(file_basename(sciencesubdirs[isubdir]),'Crisp-?') then begin
              ;; CRISP data
              red_extractstates, file_basename(snames), cam = cam, fullstate = fullstate
              cam = cam[uniq(cam,sort(cam))]
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
              printf, llun, 'There are ' + strtrim(Nscienceframes, 2) + ' ' +cam + ' science frames.'
              idark = (where(cam eq darkcams, Nhits))[0]
              if Nhits gt 0 then begin
                 printf, llun, '   There are '+strtrim(Ndarkframes[idark], 2)+' dark frames for '+cam+'.'
              endif else begin
                 printf, llun, '   Warning: No darks for '+cam+'!'
              endelse
           endif                ; darks
           


           ;; Do we have flats for these data?
           if keyword_set(flats) then begin

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

              ;; Check first if polcal is needed (more than one LC
              ;; state).

           endif                ; polcal
           
           
           ;; Do we have pinholes data?
           if keyword_set(pinholes) then begin

              
              ;; Note that we don't always have pinholes for all
              ;; wavelength points. Report nearest wavelength in that
              ;; case? 

           endif                ; pinholes

        endif                   ; Nscienceframes
     endfor                     ; isubdir   
  endfor                        ; idir
  
  free_lun, llun
  
  print, 'check_calibrations : Output in '+logfile+'.'

end
