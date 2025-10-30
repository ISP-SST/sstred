; docformat = 'rst'

;+
; Summarize collected data.
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
;    dirs : in, optional, type=strarr, default='All data and flats dirs'
;
;      Set this to the time-stamp directories that you want summarized.
; 
; :History:
; 
;   2019-08-07 : MGL. First version.
; 
;   2023-08-14 : MGL. Recognize also mosaic data.
; 
;-
pro red::summarize_datadir, dirs

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if ptr_valid(self.flat_dir)  then begin
    red_append, available_dirs, *self.flat_dir
    top_flatdir = file_dirname((*self.flat_dir)[0])
  endif
  
  if ptr_valid(self.data_dirs) then begin
    red_append, available_dirs, *self.data_dirs
    top_datadir = file_dirname((*self.data_dirs)[0])
  endif

  if n_elements(available_dirs) eq 0 then return
  
  if n_elements(dirs) eq 0 then begin
    
    use_dirs = available_dirs
    Ndirs = n_elements(use_dirs)
    if Ndirs eq 0 then return

  endif else begin

    ;; Did we give absolute paths?
    match2, dirs, available_dirs, suba, subb
    indx = where(suba ne -1, Ndirs)
    if Ndirs gt 0 then red_append, use_dirs,  available_dirs[suba[indx]]

    if n_elements(top_flatdir) gt 0 then begin
      ;; Any flats timestamps ?
      match2, top_flatdir+'/'+dirs, available_dirs, suba, subb
      indx = where(suba ne -1, Ndirs)
      if Ndirs gt 0 then red_append, use_dirs,  available_dirs[suba[indx]]
    endif
    
    if n_elements(top_datadir) gt 0 then begin
      ;; Any science-data timestamps?
      match2, top_datadir+'/'+dirs, available_dirs, suba, subb
      indx = where(suba ne -1, Ndirs)
      if Ndirs gt 0 then red_append, use_dirs,  available_dirs[suba[indx]]
    endif
    
    Ndirs = n_elements(use_dirs)
    if Ndirs eq 0 then return

  endelse

  file_mkdir, 'summaries'

  for idir = 0, Ndirs-1 do begin
    
    print, use_dirs[idir]

    undefine, output
    
    red_append, output, use_dirs[idir]
    red_append, output, ' '

    ;; Camera directories
    cams = file_basename(file_search(use_dirs[idir]+'/*'))
    Ncams = n_elements(cams)
    Nfiles_cam = lonarr(Ncams)

    camWB = (cams[where(strmatch(cams,'*-W'),     Nwbcam)])[0]
    camNB = (cams[where(strmatch(cams,'*-[NTR]'), Nnbcam)])[0]

    instrument = strupcase((strsplit(camWB, '-', /extract))[0])
    
    if Nnbcam gt 0 then begin
      ;; All NB files
      files = self -> raw_search(use_dirs[idir] + '/' + camNB, count = Nfiles)
    endif else begin
      ;; All WB files
      files = self -> raw_search(use_dirs[idir] + '/' + camWB, count = Nfiles, fpi_states = '*')
    endelse
  
    ;; Get header info 
    hdr0 = red_readhead(files[0])
    hdrN = red_readhead(files[-1])
    
    time_beg = (strsplit(fxpar(hdr0, 'DATE-BEG'), 'T', /extract))[1]
    time_end = (strsplit(fxpar(hdrN, 'DATE-END'), 'T', /extract))[1] ; Note: from last file!
    red_append, output, 'Started : ' + time_beg
    red_append, output, 'Ended   : ' + time_end
    red_append, output, ' '

    red_append, output, strtrim(Ncams, 2)+' Cameras:'
    for icam = 0, Ncams-1 do begin
      spawn, 'ls '+use_dirs[idir] + '/' + cams[icam] + ' | wc -l', Nf
      Nfiles_cam[icam] = Nf
      red_append, output, cams[icam] + ' : ' + strtrim(Nf, 2) + ' files'
    endfor
    red_append, output, ' '

    xposure = fxpar(hdr0, 'XPOSURE', count = cnt)
    if cnt gt 0 then red_append, output, 'Exposure time: '+strtrim(string(xposure*1000, format = '(f8.1)'), 2) + ' ms'

    cadence = fxpar(hdr0, 'CADENCE', count = cnt)
    if cnt eq 0 then begin
      ;; No CADENCE keyword in header. Calculate it from DATE-BEG for
      ;; the first few files
      Ncad = 30
      times = dblarr(Ncad)
      for icad = 0, Ncad-1 do $
         times[icad] = red_time2double(fxpar(red_readhead(files[icad]), 'DATE-BEG'))
                                ;times[icad] = red_time2double(fxpar(red_readhead(files0[icad]), 'DATE-BEG'))
      cadence = median(red_differential(times[sort(times)]))
    endif
    red_append, output, 'Cadence: '+strtrim(string(cadence*1000, format = '(f8.1)'), 2)+ ' ms'

    detgain = fxpar(hdr0, 'DETGAIN', count = cnt)
    if cnt gt 0 then red_append, output, 'Detector gain: '+string(detgain, format = '(f4.1)')
    red_append, output, ' '

    Nscans = fxpar(hdrN, 'SCANNUM')+1 ; Note: from last file!
    red_append, output, 'Nscans: ' + strtrim(Nscans, 2) 

    if Nnbcam ne 0 then begin
      
      ;; NB files, first scan. Don't want to extract the states
      ;; from all files!
      files0 = self -> raw_search(use_dirs[idir] + '/' + camNB $
                                  , scanno = 0, count = Nfiles0)
      self -> extractstates, files0, states0

      ;; Find all unique states, including LC for CRISP
      indx = uniq(states0.fullstate, sort(states0.fullstate))
      ;; Sort in tuning order
      indx = indx[sort(states0[indx].tun_wavelength)]
      Nstates = n_elements(indx)
      
      if max(strmatch(TAG_NAMES(states0[0]), 'LC')) eq 1 then begin
        ulc = states0[uniq(states0.lc, sort(states0.lc))].lc
        Nlc = n_elements(ulc)
      endif else Nlc = 0
      red_append, output, 'Nstates: ' + strtrim(Nstates, 2)

      ;; Mosaic data?
      pos = strpos(file_basename(files[-1]), '_mos')
      if pos ne -1 then begin
        Nmos = long(strmid(file_basename(files[-1]), pos+4, 2))+1
        red_append, output, 'Ntiles: ' + strtrim(Nmos, 2)
      endif else begin
        Nmos = 1
      endelse
      
      if Nlc gt 0 then begin
        red_append, output, 'Nlc : '+strtrim(Nlc, 2)
        red_append, output, 'Polarization states: '+strjoin(strtrim(long(ulc), 2), ', ')
      endif
      red_append, output, 'Nstates x Nscans = ' + strtrim(Nstates*Nscans, 2)
      if Nmos gt 1 then begin
        red_append, output, 'Nstates x Nscans x Ntiles = ' + strtrim(Nstates*Nscans*Nmos, 2)
      endif
      red_append, output, ' '
      
      ufullstat = states0[indx].fullstate 
      if Nlc le 1 then begin
        ustat = states0[indx].prefilter + '_' + states0[indx].tuning
      endif else begin
        ustat = states0[indx].prefilter + '_' + states0[indx].tuning + '_LC' + strtrim(long(states0[indx].lc),2)
      endelse

      upref = states0(uniq(states0.prefilter, sort(states0.prefilter))).prefilter

      Npref = n_elements(upref)
      red_append, output, strtrim(Npref, 2)+' NB prefilters: '+strjoin(strtrim(upref, 2), ', ' )
      red_append, output, ' '

      red_append, output, 'The states: '
      Nframes = 0L
      for istate = 0, Nstates-1 do begin

        ;; Some info per state
        sindx=where(ufullstat[istate] eq states0.fullstate, Nmatch)
        Nexp = 0
        if Nmatch ne 0 then begin
          sfiles = files0[sindx]
          for isfile = 0, Nmatch-1 do begin
            hdr = red_readhead(files0[indx[istate]])
            Nexp += fxpar(hdr, 'NAXIS3') >1
          endfor                ; isfile
        endif
        Nexp /= Nmos            ; Divide by number of mosaic tiles

        case instrument of
          'CHROMIS' : begin
            stat = strmid(ustat[istate]+'   ', 0, 15)
          end
          'CRISP' : begin
            ;;stat = strmid(ufullstat[istate]+'   ', 0, 19)
            stat = strmid(ustat[istate]+'   ', 0, 19)
          end
          'CRISP2' : begin
            stat = strmid(ustat[istate]+'   ', 0, 19)
          end
        endcase

        Nframes += Nexp

        red_append, output, stat + ' : ' $
                    + strtrim(Nexp, 2) + ' fr/scan, exposure time: ' $
                    + strtrim(string(xposure*Nexp,format='(f8.3)'), 2)+' s/scan'
        
      endfor                    ; istate
      red_append, output, ' '

      red_append, output, 'Nframes/scan: '+strtrim(Nframes, 2)

    endif else begin

      ;; If there are no NB files, this is probably a WB or PD flats
      ;; dir. Let's get some info from the WB files.

      files0 = self -> raw_search(use_dirs[idir] + '/' + camWB $
                                  , scanno = 0, count = Nfiles0, fpi_states = '*')
      self -> extractstates, files0, states0

      ;; Find all unique states, including LC for CRISP
      indx = uniq(states0.fullstate, sort(states0.fullstate))
      ;; Sort in tuning order
      indx = indx[sort(states0[indx].tun_wavelength)]
      Nstates = n_elements(indx)
      
      ufullstat = states0[indx].fullstate 
      ustat = states0[indx].prefilter + '_' + states0[indx].tuning

      red_append, output, 'Nstates: ' + strtrim(Nstates, 2)

      upref = states0(uniq(states0.prefilter, sort(states0.prefilter))).prefilter
      Npref = n_elements(upref)
      red_append, output, strtrim(Npref, 2)+' WB prefilters: '+strjoin(strtrim(upref, 2), ', ' )
      red_append, output, ' '

      red_append, output, 'The states: '
      Nframes = 0L
      for istate = 0, Nstates-1 do begin

        ;; Some info per state
        sindx=where(ufullstat[istate] eq states0.fullstate, Nmatch)
        Nexp = 0
        if Nmatch ne 0 then begin
          sfiles = files0[sindx]
          for isfile = 0, Nmatch-1 do begin
            hdr = red_readhead(files0[indx[istate]])
            Nexp += fxpar(hdr, 'NAXIS3') >1
          endfor                ; isfile
        endif
        

        case instrument of
          'CHROMIS' : begin
            stat = strmid(ustat[istate]+'   ', 0, 15)
          end
          'CRISP' : begin
            stat = strmid(ufullstat[istate]+'   ', 0, 19)
          end
          'CRISP2' : begin
            stat = strmid(ufullstat[istate]+'   ', 0, 19)
          end
        endcase

        Nframes += Nexp

        red_append, output, stat + ' : ' $
                    + strtrim(Nexp, 2) + ' fr/scan, exposure time: ' $
                    + strtrim(string(xposure*Nexp,format='(f8.3)'), 2)+' s/scan'
        
      endfor                    ; istate
      red_append, output, ' '

      red_append, output, 'Nframes/scan: '+strtrim(Nframes, 2)

    endelse
    
    ;; Write the file
    openw, lun, /get_lun, 'summaries/'+file_basename(use_dirs[idir])+'.txt'
    printf, lun, output, format = '(a0)'
    free_lun, lun
    
  endfor                        ; idir
  
end

;; Testing:

pwd = getenv('PWD')
case 0 of
  strmatch(pwd,'*CRISP*') : begin
    a=crispred(/dev)
    a -> summarize_datadir, ['09:30:20', '11:18:58', '11:21:22']
  end
  strmatch(pwd,'*CHROMIS*') : begin
    a=chromisred(/dev)
    ;; One data dir and one flats dir:
    a -> summarize_datadir, ['09:28:36', '11:18:58']
  end
  else : stop
endcase

end
