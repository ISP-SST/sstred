; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
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
;   all : in, optional, type=boolean
;   
;      Process all data sets.
; 
;   no_descatter : in, optional, type=boolean 
;   
;      Don't do back-scatter compensation.
;
;   tdirs : in, optional, type=strarr
;
;      Time-stamp directories, should exist in data/ subdirectory (or
;      wherever keyword link_data points to).
;
; :History:
;
;   2013-07-01 : JdlCR : Created!
; 
;   2013-07-24 : Use red_show rather than show.
;
;   2014-01-02 : PS take out nremove (selection only done in
;                prepmomfbd) use linked data, to skip incomplete
;                scans.
;
;   2015-05-05 : MGL. With a *, run all directories.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
; 
;   2016-02-16 : MGL. New keyword "all".
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding
;                SOLARNET keywords.
; 
;   2018-04-18 : MGL. Adapt to new code-base. New keyword tdir.
;
;
;-
pro red::sum_data_intdif, all = all $
                          , cam = cam $
                          , tdir = tdir $
                          , link_dir = link_dir $
                          , no_descatter = no_descatter $
                          , nthreads = nthreads $
                          , overwrite=overwrite $
                          , pref = pref $
                          , show = show $
                          , t1 = t1 $
                          , verbose = verbose 


  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Use links created by link_data as default.
  if n_elements(link_dir) eq 0 then link_dir = self.out_dir + '/' + 'data'

  if n_elements(tdir) ne 0 then begin
    eindx = where(file_test(link_dir + '/' + tdir, /directory), Nexists)
    if Nexists eq 0 then begin
      print, inam + ' : ERROR: None of the given timestamp directories exist in ' + link_dir
      print, inam + ' :        Did you run link_data?'
      if keyword_set(debug) then stop else return
    endif
    ;; At least one tdir exists so proceed with that
    dirs =  link_dir + '/' + tdir[eindx]
    Ndirs = n_elements(dirs)
    print, inam + ' : Processing given time-stamp directories in ' $
           + link_dir + ' : '
    print, file_basename(dirs)
  endif else begin
    dirs = file_search(link_dir + '/*', count = Ndirs, /test_dir)
    if Ndirs eq 0 then begin
      print, inam + ' : ERROR: No data found in ' + link_dir
      print, inam + ' :      - did you run link_data?'
      if keyword_set(debug) then stop else return
    endif
    ;; There are subdirectories in link_dir.
    if Ndirs eq 1 or keyword_set(all) then begin
      ;; We want all subdirectories (or there is only one), so proceed
      ;; with them (it).
      print, inam + ' : Processing all time-stamp directories in ' $
             + link_dir + ' : '
      print, file_basename(dirs)
    endif else begin
      ;; Make a selection from the available time-stamp directories
      tmp = red_select_subset(file_basename(dirs)$
                              , qstring = inam + ' : Select timestamp directory ID:' $
                              , count = Ndirs, indx = sindx)
      if Ndirs eq 0 then begin
        print, inam + ' : No timestamp sub-folders selected.'
        return                  ; Nothing more to do
      endif 
      dirs = dirs[sindx]
      print, inam + ' : Processing selected time-stamp directories in ' $
             + link_dir + ' : '
      print, file_basename(dirs)
    endelse
    
  endelse
  


  ;; Find the narrowband cameras.
  self -> getdetectors
  cams = *self.cameras
  detectors = *self.detectors
  is_wb = strmatch(cams,'*-[DW]')
  is_nb = ~is_wb
  cams = cams[where(~is_wb, Nnb)]
  detectors = detectors[where(~is_wb, Nnb)]
  if Nnb eq 0 then stop
  Ncams = n_elements(cams)

  for idir = 0, Ndirs-1 do begin

    dir = dirs[idir]

;    ;; Get camera tags
;    self -> getdetectors, dir = dir
;    ctag = [self.camrtag, self.camttag]

    outdir = self.out_dir + '/cmap_intdif/' + file_basename(dir) + '/'
    file_mkdir, outdir
    
;    files = file_search(dir+ '/' + cams[0] + '/cam*', count = Nfiles)
    files = file_search(dir+ '/' + cams + '/cam*', count = Nfiles)

    if Nfiles eq 0 then begin
      print, inam + ' : ERROR, data folder is empty'
      if keyword_set(debug) then stop else return
    endif
    
 ;   files = red_sortfiles(temporary(files))

    ;; Extract tags from file names
;    state = red_getstates(files, /links)
    self -> extractstates, files, states
    
    ;; Remove frames
    ;;red_flagtuning, state, remove

    ;; Build my states
    mpref = states.prefilter
    upref = reform([mpref[uniq(mpref, sort(mpref))]])
    Npref = n_elements(upref)

    sel_pref = 1B
    if n_elements(pref) gt 0 then begin
      pos = where(upref eq pref, count)
      if count eq 0 then begin
        print, inam + ' : User supplied prefilter '+pref+' not found in '+dir
        continue
      endif else sel_pref = 0B
    endif 

    if(sel_pref) then begin
      if(Npref eq 1) then begin
        pref = upref[0]
      endif else begin
        while Npref ne 1 do begin
          tmp = red_select_subset(upref $
                                  , qstring = inam + ' : Select a single prefilter:' $
                                  , count = Npref, indx = sindx)
        endwhile
        pref = pref[sindx[0]]
;        print, inam + ' : Found prefilters:'
;        for ii = 0, Npref - 1 do print, string(ii, format='(I3)') + ' -> ' + upref[ii]
;        ip = 0L
;        read, ip, prompt = 'Select state id: '
      endelse
    endif
    print, inam + ' : Selected prefilter -> ' + pref
    
    ;; Isolate files from prefilter
    self->selectfiles, files = files, states = states $
                       , prefilter = pref, sel = sel, count = Nselect
    if Nselect eq 0 then stop
    
    selstates = states[sel]
    selfiles = files[sel]

;    idx = where(mpref eq pref, count)
;    mwav = state.wav[idx]
;    mdwav = state.dwav[idx]
;    mlc = state.lc[idx]
;    mscan = state.scan[idx]
;    mnum = state.nums[idx]
;    mfiles = state.files[idx]
;    mstar = state.star[idx]

;    uscan = mscan[uniq(mscan, sort(mscan))]
;    ulc = mlc[uniq(mlc, sort(mlc))]
;    uwav = mwav[uniq(mwav, sort(mdwav))]
;    udwav = float(mdwav[uniq(mdwav, sort(mdwav))] - double(pref))

    uscan = selstates[uniq(selstates.scannumber, sort(selstates.scannumber))].scannumber
    Nscans = n_elements(uscan)

    ulc = selstates[uniq(selstates.lc, sort(selstates.lc))].lc
    Nlc = n_elements(ulc)

    ;; Sorted in wavelength order
    indx = uniq(selstates.tun_wavelength, sort(selstates.tun_wavelength))
    Ntunings = n_elements(indx)
    uwav = selstates[indx].tuning
    udwav = selstates[indx].tun_wavelength

    ;; Start summing
    if n_elements(t1) ne 0 then begin
      lscan = long(uscan)
      tt1 = where(lscan eq t1, count)
      if count eq 0 then begin
        print, inam + ' : Scan', t1, ' not found -> ending at t1 = ', Nscans-1
        tt1 = Nscans-1L
      endif
    endif else tt1 = Nscans-1

    
    for icam = 0, Ncams-1 do begin
      ;; Did we select a camera?
      if n_elements(cam) ne 0 then begin
        if cams[icam] ne cam then begin
          print, inam + ' : Skipping cam -> '+cams[icam]+' != '+cam
          continue
        endif
      endif

      self -> selectfiles, files = selstates.filename, states = selstates $
                           , cam = cams[icam] $
                           , sel = sell, count = Nselect
      if Nselect eq 0 then stop

      
      ;; Descatter?
      if ~keyword_set(no_descatter) $
         and self.dodescatter $
         and (pref eq '8542' or pref eq '7772') then begin
        self -> loadbackscatter, detectors[icam], pref, bff, pff
      endif

      
;      tempdir = dir + '/' + cams[icam]+'/'
;      longi = strlen(file_dirname(mfiles[0])+'/'+self.camrtag)
;      mmfiles = tempdir + detectors[icam]+strmid(mfiles,longi,200)
      mstates = selstates[sell]
;      mmfiles = selstates[sell]

      mwav   = mstates.tuning
      mlc    = mstates.lc
      mscan  = mstates.scannumber
      mfiles = mstates.filename
      mstar  = lonarr(Nselect)  ; all zeros

      cfile = outdir + '/' + detectors[icam] + '.'+ pref + '.intdif.icube'
      dfile = outdir + '/' + detectors[icam] + '.'+ pref + '.intdif.save'

      ;; We could change the icube + save files to fitscubes with the save
      ;; info in the header/extensions.
      
      if file_test(dfile) and ~keyword_set(overwrite) then begin
        print, inam + ' : WARNING! Data files already exist:'
        print, dfile
        print, cfile
        print, inam + ' : You can overwrite them with /overwrite, otherwise we will resume the summing!'
        restore, dfile
        idx = where(done eq 0, count)
        if count eq 0 then begin
          print, inam + ' : Nothing to do for '+detectors[icam]
          continue
        endif else tt0 = idx[0]
        openu, lun, cfile, /get_lun
      endif else begin
        done = bytarr(Nscans)
        openw, lun, cfile, /get_lun
        tt0 = 0L
      endelse
;      dim = size(f0(mmfiles[0]),/dim)
      hdr = red_readhead(mstates[0].filename)
      dim = fxpar(hdr, 'NAXIS*')
      Nx = dim[0]
      Ny = dim[1]

      head = red_pol_lpheader(Nx, Ny, Nlc*(tt1+1L)*Ntunings)
      writeu,lun, head
      dat = assoc(lun, intarr(Nx, Ny, /nozer), 512)

      ;; Load dark file
;      self -> get_calib, selstates[0], darkdata = dd
      self -> get_calib, mstates[0], darkdata = dd

      ;; Load fitgains results
      cmf = self.out_dir + '/flats/spectral_flats/' + detectors[icam] + '.' + $
            pref + '.' + 'fit_results.sav'
      
      ;; Loop wavelengths within the scan
      for iscan = tt0[0], tt1[0] do begin
        
        for ilc = 0L, Nlc - 1 do begin

;          scstate = strarr(Ntunings)

          for ituning = 0L, Ntunings - 1 do begin

;            istate = strjoin([uscan[iscan], pref,uwav[ituning], ulc[ilc]],'.')
;            scstate[ituning] = istate
            
            idx = where((mwav EQ uwav[ituning]) AND (mlc EQ ulc[ilc]) AND $
                        (mscan EQ uscan[iscan]) AND mstar eq 0, count)

            if(count eq 0) then begin
              print, inam + ' : No files found for state ' + mstates[idx[0]].fullstate
            endif else begin
              print, inam + ' : Found ' + red_stri(count) + ' images for state ' + mstates[idx[0]].fullstate
            endelse


            ;; Read cavity free gains
            self -> get_calib, mstates[idx[0]], cgaindata = gg
            
            if keyword_set(verbose) then print, file_basename(mstates[idx].filename),format='(a0)'
            imsum = rdx_sumfiles(mstates[idx].filename, /check, nthreads = 2) - dd

            if ~keyword_set(no_descatter) $
               && self.dodescatter $
               && (pref eq '8542' || pref eq '7772') then begin
              imsum = rdx_descatter(temporary(imsum), bff, pff, /verbose, nthreads = nthreads)
            endif
            
            imsum = red_fillpix(temporary(imsum)*gg, nthreads=nthreads)

            ele = iscan*Nlc*Ntunings + ilc*Ntunings + ituning
            dat[ele] = fix(round(7.0 * imsum))

            if keyword_set(show) then begin
              if n_elements(mydum) eq 0 then begin
                mydum = 1
                red_show, red_histo_opt(imsum)
              endif
              red_show, red_histo_opt(imsum), /nowin
            endif

          endfor                ; ituning
        endfor                  ; ilc

        done[iscan] = 1B
        ;; Save incomplete results
        save, file=dfile, done, uwav, ulc, uscan, Ntunings, Nlc, Nscans, pref, Nx, Ny, udwav

      endfor                    ; iscan

      free_lun, lun

    endfor                      ; icam

    ;; Save final results
    save, file=dfile, done, uwav, ulc, uscan, Ntunings, Nlc, Nscans, pref, Nx, Ny, udwav

  endfor                        ; idir

end

