; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
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
;    aligncontkludge : in, optional, type=boolean
;
;       Set this to activate a kludge where the misaligned Ca
;       continuum narrowband image is stretched to its wideband image
;       before applying the stretching that aligns it to the anchor
;       wideband image.
;
;    ny_limit : in, optional, type=integer
;
;       Used in red_flipthecube_unpol to help determine whether the
;       cube should be flipped in memory.
;
;    rot_dir  : 
;   
;   
;   
;    square  : 
;   
;   
;   
;    tiles : 
;   
;   
;   
;    clips : 
;   
;   
;   
;    scans_only : in, type=boolean
;   
;       Set this to save data on a per-scan basis, not per time series.
;   
;    overwrite  : in, optional, type=boolean
;   
;       Set this to overwrite existing files.   
;
;    wbwrite : in, type=boolean
;
;       Set this to write also the global wideband image to the
;       crispex directory. (So far only implemented for /scans_only.) 
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2013-07-05 : Allow to bypass flat-ratio operations
;
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
;   2013-07-11 : MGL. Added keyword wbwrite. Set this to write also
;                the global wideband image to disk. So far only
;                implemented for /scans_only. 
; 
;   2013-07-12 : MGL. Bugfixes. Calculate cscl also when we skip the
;                first scan because it's already been processed.
;
;   2013-08-19 : JdlCR. Spectfile is created along with the crispex
;                cubes.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
;
;   2014-11-29 : JdlCR, added support for fullframe cubes (aka,
;                despite rotation and shifts, the entire FOV is inside
;                the image
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-09-26 : MGL. Adapt red::make_unpol_crispex for CHROMIS data.
;
;   2016-09-28 : MGL. New keyword aligncontkludge. New output file
;                names including date and scan selection. 
;
;   2016-11-02 : MGL. Changed default clips and tiles. Now works with
;                selection of scans done in polish_tseries.
;
;   2016-11-17 : MGL. Changed default clips and tiles again.
;   
;   2016-11-17 : JdlCR. Added prefilter calibration and absolute units
;                calibration. Removed the /noflat keyword for the time
;                being.
;
;   2016-11-28 : MGL. Renamed keyword aligncontkludge to aligncont,
;                make it align not only the continuum image but all nb
;                images based on output from earlier pipeline step
;                align_continuum.
;
;   2016-12-04 : JdlCR. Fixed contalign bug when only one scan is present.
;
;   2017-04-18 : MGL. New keyword momfbddir.
;
;   2017-05-24 : MGL. Add prpara info.
;
;   2017-12-13 : MGL. New keyword Ny_limit.
;
;-
pro chromis::make_crispex, aligncont = aligncont $
                           , clips=clips $
                           , float = float $
                           , momfbddir = momfbddir $
                           , no_timecor=no_timecor $
                           , nostretch=nostretch $ $
                           , Ny_limit = Ny_limit $
                           , overwrite = overwrite $
                           , rot_dir = rot_dir $
                           , scans_only = scans_only $
                           , selscan=selscan $
                           , tiles=tiles $
                           , verbose=verbose $
                           , wbwrite = wbwrite $
                           , square = square 
;                           , noflats=noflats $
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Make prpara
  if n_elements(aligncont   ) ne 0 then red_make_prpara, prpara, 'aligncont'    , aligncont 
  if n_elements(clips       ) ne 0 then red_make_prpara, prpara, 'clips'        , clips         
  if n_elements(float       ) ne 0 then red_make_prpara, prpara, 'float'        , float  
  if n_elements(momfbddir   ) ne 0 then red_make_prpara, prpara, 'momfbddir'    , momfbddir    
  if n_elements(no_timecor  ) ne 0 then red_make_prpara, prpara, 'no_timecor'   , no_timecor 
  if n_elements(nostretch   ) ne 0 then red_make_prpara, prpara, 'nostretch'    , nostretch 
  if n_elements(np          ) ne 0 then red_make_prpara, prpara, 'np'           , np           
  if n_elements(overwrite   ) ne 0 then red_make_prpara, prpara, 'overwrite'    , overwrite
  if n_elements(rot_dir     ) ne 0 then red_make_prpara, prpara, 'rot_dir'      , rot_dir         
  if n_elements(scans_only  ) ne 0 then red_make_prpara, prpara, 'scans_only'   , scans_only          
  if n_elements(selscan     ) ne 0 then red_make_prpara, prpara, 'selscan'      , selscan 
  if n_elements(square      ) ne 0 then red_make_prpara, prpara, 'square'       , square       
  if n_elements(tiles       ) ne 0 then red_make_prpara, prpara, 'tiles'        , tiles        
  if n_elements(wbwrite     ) ne 0 then red_make_prpara, prpara, 'wbwrite'      , wbwrite

  ;; Default keywords
  if n_elements(momfbddir) eq 0 then momfbddir = 'momfbd' 
  if n_elements(rot_dir) eq 0 then rot_dir = 0B

  ;; Camera/detector identification
  self->getdetectors
  wbindx = where(strmatch(*self.cameras,'Chromis-W'))
  wbcamera = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx = where(strmatch(*self.cameras,'Chromis-N')) 
  nbcamera = (*self.cameras)[nbindx[0]]
  nbdetector = (*self.detectors)[nbindx[0]]
  ;; Should be generalized to multiple NB cameras if CHROMIS gets
  ;; polarimetry. We don't need to identify any PD cameras for
  ;; restored data.

  ;; Find timestamp subdirs
  search_dir = self.out_dir +'/'+momfbddir+'/'
  timestamps = file_basename(file_search(search_dir + '*' $
                                         , count = Ntimestamps, /test_dir))
  if Ntimestamps eq 0 then begin
    print, inam + ' : No timestamp sub-directories found in: ' + search_dir
    return
  endif

  ;; Select timestamp folders
  selectionlist = strtrim(indgen(Ntimestamps), 2)+ '  -> ' + timestamps
  tmp = red_select_subset(selectionlist $
                          , qstring = inam + ' : Select timestamp directory ID:' $
                          , count = Ntimestamps, indx = sindx)
  if Ntimestamps eq 0 then begin
    print, inam + ' : No timestamp sub-folders selected.'
    return                      ; Nothing more to do
  endif
  timestamps = timestamps[sindx]
  print, inam + ' : Selected -> '+ strjoin(timestamps, ', ')
  
  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0
  red_logdata, self.isodate, time_pig, pig = metadata_pig, rsun = rsun

  ;; Loop over timestamp directories
  for itimestamp = 0L, Ntimestamps-1 do begin

    timestamp = timestamps[itimestamp]
    datestamp = self.isodate+'T'+timestamp

    ;; Find prefilter subdirs
    search_dir = self.out_dir +'/'+momfbddir+'/'+timestamp+'/'
    prefilters = file_basename(file_search(search_dir + '*' $
                                           , count = Nprefs, /test_dir))
    if Nprefs eq 0 then begin
      print, inam + ' : No prefilter sub-directories found in: ' + search_dir
      continue                  ; Next timestamp
    endif
    
    ;; Select prefilter folders
    selectionlist = strtrim(indgen(Nprefs), 2)+ '  -> ' + prefilters
    tmp = red_select_subset(selectionlist $
                            , qstring = inam + ' : Select prefilter directory ID:' $
                            , count = Nprefs, indx = sindx)
    if Nprefs eq 0 then begin
      print, inam + ' : No prefilter sub-folders selected.'
      continue                  ; Go to next timestamp
    endif
    prefilters = prefilters[sindx]
    print, inam + ' : Selected -> '+ strjoin(prefilters, ', ')

    ;; Loop over WB prefilters
    for ipref = 0L, Nprefs-1 do begin


      search_dir = self.out_dir + '/'+momfbddir+'/' + timestamp $
                   + '/' + prefilters[ipref] + '/cfg/results/'
      files = file_search(search_dir + '*.'+['f0', 'momfbd', 'fits'], count = Nfiles)      
      
      ;; Find all the global WB images 
;      self -> selectfiles, files = files, states = states $
;                           , cam = wbcamera, ustat = '' $
;                           , sel = wbgindx, count = Nscans $
;                           , complement = complement, Ncomplement = Ncomplement
      ;; We have no special state (or absence of state) to identify
      ;; the global WB images but we do know that their exposure times
      ;; are much larger than the ones corresponding to the individual
      ;; NB states.
      self -> extractstates, files, states
      wbgindx = where(states.exposure gt mean(states.exposure) $
                      , Nscans, complement = complement, Ncomplement = Ncomplement)
      wbgstates = states[wbgindx]
      wbgfiles = files[wbgindx]

      ;; All the per-tuning files and states
      pertuningfiles = files[complement]
      pertuningstates = states[complement]

      ;; Unique tuning states
      utunindx = uniq(pertuningstates.fpi_state, sort(pertuningstates.fpi_state))
      sortindx = sort(pertuningstates[utunindx].tun_wavelength)
      ufpi_states = pertuningstates[utunindx[sortindx]].fpi_state
      utunwavelength = pertuningstates[utunindx[sortindx]].tun_wavelength

      wav = utunwavelength
      my_prefilters = pertuningstates[utunindx[sortindx]].prefilter
      Nwav = n_elements(utunindx)

      ;; Unique nb prefilters
      unbprefindx = uniq(pertuningstates[utunindx].prefilter, sort(pertuningstates[utunindx].prefilter))
      Nnbprefs = n_elements(unbprefindx)
      unbprefs = pertuningstates[utunindx[unbprefindx]].prefilter
      unbprefsref = dblarr(Nnbprefs)
      for inbpref = 0L, Nnbprefs-1 do begin
        ;; This is the reference point of the fine tuning for this prefilter:
        unbprefsref[inbpref] = double((strsplit(pertuningstates[utunindx[unbprefindx[inbpref]]].tuning $
                                                , '_', /extract))[0])
      endfor
      unbprefsref *= 1e-10      ; [m]

      if ~keyword_set(scans_only) then begin
        ;; Look for time-series calib file
        csearch = self.out_dir + '/calib_tseries/tseries_' + prefilters[ipref] $
                  + '_' + datestamp + '*_calib.sav'
        cfiles = file_search(csearch, count = Ncfiles)
        case Ncfiles of
          0: begin
            print, inam + ' : Could not find calibration file: ' + csearch
            print, inam + ' : Try executing red::polish_tseries on this dataset first!'
            return
          end
          1: cfile = cfiles[0]
          else: begin
            repeat begin
              tmp = red_select_subset(cfiles $
                                      , qstring = inam + ' : Select calibration file (scan subset).' $
                                      , count = Ncfileselect, indx = cindx, default = '-')
            endrep until Ncfileselect eq 1
            cfile = cfiles[cindx[0]]
          end
        endcase

        print, inam + ' : Loading calibration file -> '+file_basename(cfile)
        restore, cfile
        if(n_elements(ff) eq 5) then full = 1 else full = 0

        ;;tmean = mean(tmean) / tmean

        ;; Get the scan selection from wfiles (from the sav file)
        wbgfiles = wfiles
        self -> extractstates, wbgfiles, wbgstates
        uscans = wbgstates.scannumber
        Nscans = n_elements(uscans)
      endif else begin
        full = 0
        uscans = wbgstates[uniq(wbgstates.scannumber, sort(wbgstates.scannumber))].scannumber
        Nscans = n_elements(uscans)
        tmean = replicate(1.0, Nscans) ; Dummy time correction
      endelse

      ;; Per-tuning files, wb and nb, only for selected scans
      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                           , scan = uscans $
                           , cam = wbcamera $
                           , sel = wbindx, count = Nwb
      wbstates = pertuningstates[wbindx]
      wbfiles = pertuningfiles[wbindx]
      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                           , scan = uscans $
                           , cam = nbcamera $
                           , sel = nbindx, count = Nnb
      nbstates = pertuningstates[nbindx]
      nbfiles = pertuningfiles[nbindx]

      ;; Prepare for making output file names
      midpart = prefilters[ipref] + '_' + datestamp + '_scans=' $ 
                + red_collapserange(uscans, ld = '', rd = '')

      ;; Load prefilters
      for inbpref = 0L, Nnbprefs-1 do begin
        pfile = self.out_dir + '/prefilter_fits/chromis_'+unbprefs[inbpref]+'_prefilter.idlsave'
        if ~file_test(pfile) then begin
          print, inam + ' : prefilter file not found: '+pfile
          return
        endif
        
        restore, pfile          ; Restores variable prf which is a struct
        idxpref = where(my_prefilters eq unbprefs[inbpref], count)

        if inbpref eq 0 then begin
          prefilter_curve = [0.d0]
          prefilter_wav = [0.0d0]
          prefilter_wb = [0.0d0]
          units = prf.units
        endif else begin
          if units ne prf.units then begin
            print, inam + ' : Units in ' + pfile + ' do not match those in earlier read files.'
            print, inam + ' : Please rerun the prefilterfit step for these data.'
            retall
          endif
        endelse
        
        if(count eq 1) then begin
          prefilter_curve = [prefilter_curve, prf.pref]
          prefilter_wav = [prefilter_wav, prf.wav]
          prefilter_wb = [prefilter_wb, prf.wbint]
        endif else begin
          me = median(prf.wav)
          prefilter_curve = [prefilter_curve, red_intepf(prf.wav-me, prf.pref, wav[idxpref]*1.d10-me)]
          prefilter_wav = [prefilter_wav, wav[idxpref]*1.d10]
          prefilter_wb = [prefilter_wb, replicate(prf.wbint, count)]
        endelse
      endfor                    ; inbpref
      
      rpref = 1.d0/prefilter_curve[1:*]
      prefilter_wav = prefilter_wav[1:*]
      prefilter_wb = prefilter_wb[1:*]
      prefilter_curve = prefilter_curve[1:*]

      ;; Do WB correction?
      if(Nwb eq Nnb) then wbcor = 1B else wbcor = 0B


      
      ;; if(~keyword_set(noflats)) then begin
      ;;   ;; Load clean flats and gains
      ;;   for ii = 0L, Nwav-1 do begin
      
      ;;     tff = self.out_dir + 'flats/' + self.camttag + '.' + prefilters[ipref] $
      ;;           + '.'+st.uwav[ii]+'.unpol.flat'
      ;;     rff = self.out_dir + 'flats/' + self.camrtag + '.' + prefilters[ipref] $
      ;;           + '.'+st.uwav[ii]+'.unpol.flat'
      ;;     tgg = self.out_dir + 'gaintables/' + self.camttag + '.' + prefilters[ipref] $
      ;;           + '.' + st.uwav[ii] + '.lc4.gain'
      ;;     rgg = self.out_dir + 'gaintables/' + self.camrtag + '.'+prefilters[ipref] $
      ;;           + '.' + st.uwav[ii] + '.lc4.gain'
      ;;     if(ii eq 0) then print, inam + ' : Loading: '
      ;;     print,' -> '+tff
      ;;     print,' -> '+rff
      ;;     print,' -> '+tgg
      ;;     print,' -> '+rgg
      
      ;;     if ~file_test(tff) $
      ;;        or ~file_test(rff) $
      ;;        or ~file_test(tgg) $
      ;;        or ~file_test(rgg) then begin
      ;;       print, inam + ' : ERROR -> Flat/gain files not found'
      ;;       return
      ;;     endif
      
      ;;     if(ii eq 0) then begin
      ;;       dim = size(f0(tff),/dimen)
      ;;       tratio = fltarr(dim[0], dim[1], nwav)
      ;;       rratio = fltarr(dim[0], dim[1], nwav)
      ;;     endif 
      
      ;;     tmp = f0(tff)
      ;;     tmp1 = f0(tgg)
      ;;     idx = where(tmp1 gt 0.0001, count, complement = idx1)
      ;;     mask = bytarr(dim) + 1B
      ;;     tmp[idx] = mean(tmp[idx]) / tmp[idx]
      ;;     tmp[idx] = tmp[idx] / tmp1[idx]
      ;;     if(n_elements(idx1) gt 0) then begin
      ;;       tmp[idx1] = 0.0d0
      ;;       mask[idx1] = 0B
      ;;     endif
      ;;     tmp = red_fillpix(tmp, mask=red_cleanmask(mask),nthreads=6)
      ;;     tratio[*,*,ii] = temporary(tmp)
      
      ;;     tmp = f0(rff)
      ;;     tmp1 = f0(rgg)
      ;;     idx = where(tmp1 gt 0.0001, count, complement = idx1)
      ;;     mask = bytarr(dim) + 1B
      ;;     tmp[idx] = mean(tmp[idx]) / tmp[idx]
      ;;     tmp[idx] = tmp[idx] / tmp1[idx]
      ;;     if(n_elements(idx1) gt 0) then begin
      ;;       tmp[idx1] = 0.0d0
      ;;       mask[idx1] = 0B
      ;;     endif
      ;;     tmp = red_fillpix(tmp, mask = red_cleanmask(mask),nthreads=6)
      ;;     rratio[*,*,ii] = temporary(tmp)
      ;;   endfor
      ;; endif                     ; ~noflats

      ;; Load WB image and define the image border
      tmp = red_readdata(wbgfiles[0])

      if keyword_set(scans_only) then begin
        dimim = red_getborder(tmp, x0, x1, y0, y1, square=square)
      endif else begin
        dimim = [x1-x0+1, y1-y0+1]
      endelse
      
      if full then begin
        dimim[0] = nd[0]
        dimim[1] = nd[1]
      endif
      Nx = dimim[0]
      Ny = dimim[1]

      ;; Create temporary cube and open output file
      d = fltarr(dimim[0], dimim[1], Nwav)  

      if(n_elements(odir) eq 0) then odir = self.out_dir + '/crispex/' + timestamp + '/'
      file_mkdir, odir

      
      if keyword_set(float) then extent = '.fcube' else extent = '.icube'

      if(keyword_set(scans_only)) then begin
        head = red_unpol_lpheader(dimim[0], dimim[1], Nwav, float = float)
      endif else begin
        head = red_unpol_lpheader(dimim[0], dimim[1], Nwav*Nscans, float = float)

        ;; Open assoc file for output of multi-scan data cube.
        ofile = 'crispex_'+midpart+'_time-corrected'+extent

        if file_test(odir + '/' + ofile) then begin
          if keyword_set(overwrite) then begin
            print, 'Overwriting existing data cube:'
            print, odir + '/' + ofile
          endif else begin
            print, 'This data cube exists already:'
            print, odir + '/' + ofile
            return
          endelse
        endif
        
        openw, lun, odir + '/' + ofile, /get_lun
        writeu, lun, head
        point_lun, lun, 0L
        print, inam+' assoc file -> ',  odir + '/' + file_basename(ofile,extent)+'.assoc.pro'
        openw, lunf, odir + '/' + file_basename(ofile,extent)+'.assoc.pro', /get_lun
        printf,lunf, 'nx=', dimim[0]
        printf,lunf, 'ny=', dimim[1]
        printf,lunf, 'nw=', Nwav
        printf,lunf, 'nt=', Nscans
        printf,lunf, "openr,lun,'"+ofile+"', /get_lun"
        if keyword_set(float) then begin
          dat = assoc(lun, fltarr(dimim[0], dimim[1], nwav, /nozero), 512)
          printf,lunf, "dat = assoc(lun, fltarr(nx,ny,nw,/nozer), 512)"
        endif else begin
          dat = assoc(lun, intarr(dimim[0], dimim[1], nwav, /nozero), 512)
          printf,lunf, "dat = assoc(lun, intarr(nx,ny,nw,/nozer), 512)"
        endelse
        free_lun, lunf
      endelse

      ;; Start processing data
      if(~keyword_set(tiles) OR (~keyword_set(clips))) then begin
        tiles = [8, 16, 32, 64, 128]
        clips = [8, 4, 2, 1, 1]
      endif

      
      if keyword_set(aligncont) then begin

        ;; Get shifts based on continuum vs wideband alignment.

        aligndir = self.out_dir + '/align/' + timestamp $
                   + '/' + prefilters[ipref] + '/'
        
        nname = aligndir+'scan_numbers.fz'
        sname = aligndir+'continuum_shifts_smoothed.fz'
        
        if ~file_test(nname) or ~file_test(sname) then begin
          print, inam + ' : At least one file missing for aligncont option:'
          print, nname
          print, sname
          retall
        endif
        
        ;; Read the shifts for the continuum images
        align_scannumbers = f0(nname)
        align_shifts = f0(sname)

        ;; Use interpolation to get the shifts for the selected scans.
        nb_shifts = fltarr(2, Nscans)
        for bb=0, Nscans-1 do begin
          pos = where(align_scannumbers eq uscans[bb], cccc)
          if(cccc eq 1) then nb_shifts[*, bb] = align_shifts[*, pos] else begin
            nb_shifts[0, *] = interpol([reform(align_shifts[0, *])], [float(align_scannumbers)], [float(uscans)])
            nb_shifts[1, *] = interpol([reform(align_shifts[1, *])], [float(align_scannumbers)], [float(uscans)])
          endelse
        endfor
        pos = where(~finite(nb_shifts), cccc)
        if(cccc gt 0) then nb_shifts[pos] = 0
      endif


      iprogress = 0
      Nprogress = Nscans*Nwav
      for iscan = 0L, Nscans-1 do begin

        if(n_elements(selscan) gt 0) then if selscan ne strtrim(uscans[iscan], 2) then continue
;        print, inam + ' : processing scan -> '+strtrim(uscans[iscan], 2)

        ;; Save the wavelength points in a separate file, common to
        ;; all the scans.
        if(iscan eq 0) then begin
                                ;         wav = scan_nbstates.tun_wavelength
          fzwrite, wav, odir + '/' + 'wav_' + prefilters[ipref] +'.f0',' '
        endif

        ;; The NB files in this scan, sorted in tuning wavelength order.
        self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                             , cam = nbcamera, scan = uscans[iscan] $
                             , sel = scan_nbindx, count = count
        scan_nbfiles = pertuningfiles[scan_nbindx]
        scan_nbstates = pertuningstates[scan_nbindx]
        sortindx = sort(scan_nbstates.tun_wavelength)
        scan_nbfiles = scan_nbfiles[sortindx]
        scan_nbstates = scan_nbstates[sortindx]

        ;; The WB files in this scan, sorted as the NB files
        self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                             , cam = wbcamera, scan = uscans[iscan] $
                             , sel = scan_wbindx, count = count
        scan_wbfiles = pertuningfiles[scan_wbindx]
        scan_wbstates = pertuningstates[scan_wbindx]
        match2, scan_nbstates.fpi_state, scan_wbstates.fpi_state, sortindx
;        sortindx = sort(scan_wbstates.tun_wavelength)
        scan_wbfiles = scan_wbfiles[sortindx]
        scan_wbstates = scan_wbstates[sortindx]
        

        
        if(keyword_set(scans_only)) then begin
          ofile = 'crispex_' + prefilters[ipref] + '_' + datestamp + '_scan=' $
                  + string(uscans[iscan], format = '(i05)') + extent
          ofilewb = 'wb_' + prefilters[ipref] + '_' + datestamp + '_scan=' $
                    + string(uscans[iscan], format = '(i05)') + '.fz' 
          if file_test(odir + '/' + ofile) then begin
            if keyword_set(overwrite) then begin
              print, 'Overwriting existing data cube:'
              print, odir + '/' + ofile
            endif else begin
              print, 'Skip to next scan, this one exists already:'
              print, odir + '/' + ofile
              continue          ; Skip to next iteration of "for iscan ..." loop.
            endelse
          endif
        endif
        
        wb = (red_readdata(wbgfiles[iscan]))[x0:x1, y0:y1]

        if keyword_set(aligncont) then begin
          
          ;; Interpolate to get the shifts for all wavelengths for
          ;; this scan.

          ;; Get these from somewhere else (scan_wbstates does not
          ;; have the WB tuning or prefilter!):
          lambdaW = 3950e-10
          lambdaC = 3998.640d-10 + 1.258d-10
          
          xshifts = interpol([0., nb_shifts[0, iscan]], [lambdaW, lambdaC]*1e7 $
                             , scan_nbstates.tun_wavelength*1e7)
          yshifts = interpol([0., nb_shifts[1, iscan]], [lambdaW, lambdaC]*1e7 $
                             , scan_nbstates.tun_wavelength*1e7)
          
        endif

        for iwav = 0L, Nwav - 1 do begin 
          ;; state = strjoin((strsplit(file_basename(st.ofiles[iwav,iscan]),'.',/extract))[1:*],'.')
          state = ufpi_states[iwav]

          ;; Collect info about this frame here.
          
          nbhead = red_readhead(scan_nbfiles[iwav])

          red_progressbar, iprogress, Nprogress $
                           , /predict $
                           , 'Processing scan ' $
                           + strtrim(uscans[iscan], 2) + ' state=' + state 

          ;; Get destretch to anchor camera (residual seeing)
          if(wbcor) then begin
            wwi = (red_readdata(scan_wbfiles[iwav]))[x0:x1, y0:y1]
            grid1 = red_dsgridnest(wb, wwi, tiles, clips)
          endif

          if 0 then wwc = red_stretch(wwi, grid1) ; test

          ;; Read image and apply prefilter curve
          tmp = (red_readdata(scan_nbfiles[iwav]))[x0:x1, y0:y1] * rpref[iwav]

;          if (self.filetype eq 'ANA') then begin
;            tmp0 = (f0(ttf))[x0:x1, y0:y1] * tpref[iwav]
;            tmp1 = (f0(rrf))[x0:x1, y0:y1] * rpref[iwav]
;          endif else begin
;            tmp_raw0 = momfbd_read(ttf)
;            tmp_raw1 = momfbd_read(rrf)
;
;            tmp0 = (red_mozaic(tmp_raw0))[x0:x1, y0:y1] * tpref[iwav]
;            tmp1 = (red_mozaic(tmp_raw1))[x0:x1, y0:y1] * rpref[iwav]


;; The following part is not ported yet:
;; ***********************************************************
;            ;; Apply flat ratio after convolving with the PSF of the
;            ;; patch: red_img2momfbd
;            if(~keyword_set(noflats)) then begin
;              trat = (red_mozaic(red_img2momfbd(tmp_raw0, tratio[*,*,iwav])))[x0:x1, y0:y1]
;              rrat = (red_mozaic(red_img2momfbd(tmp_raw1, rratio[*,*,iwav])))[x0:x1, y0:y1]
;              
;              tmp0 = temporary(tmp0) * temporary(trat) 
;              tmp1 = temporary(tmp1) * temporary(rrat) 
;            endif
;; ***********************************************************



;          endelse 

;          ;; Combine cameras, compute scale factor avoiding borders...
;          dim = size(tmp0,/dim)
;          xx0 = round(dim[0] * 0.15)
;          xx1 = round(dim[0] * 0.85)
;          yy0 = round(dim[1] * 0.15)
;          yy1 = round(dim[1] * 0.85)
;          
;          me = median(tmp0[xx0:xx1,yy0:yy1] + tmp1[xx0:xx1,yy0:yy1]) * 0.5
;          sclt = me / (median(tmp0[xx0:xx1,yy0:yy1]))
;          sclr = me / (median(tmp1[xx0:xx1,yy0:yy1]))
          
;          tmp = (temporary(tmp0) * sclt + temporary(tmp1) * sclr) 
          
          if keyword_set(aligncont) then begin
            
            tmp = red_shift_sub(tmp, -xshifts[iwav], -yshifts[iwav])

;            ;; This is the continuum point for a Ca scan, has to be
;            ;; different for a Hb scan:
;            continnumpoint = scan_nbstates[iwav].fpi_state eq '3999_4000_+0'
;            if continnumpoint then begin
;              ;; Stretch the nb cont image to its wb image
;              gridx = red_dsgridnest(wwi, tmp, tiles, clips)
;              tmp = red_stretch(temporary(tmp), gridx)
;            endif
          endif

          ;; Apply destretch to anchor camera and prefilter correction
          if(wbcor) then tmp = red_stretch(temporary(tmp), grid1)
          
          if(~keyword_set(scans_only)) then begin
            ;; Apply derot, align, dewarp based on the output from
            ;; polish_tseries
            if(full) then begin
              bla = red_rotation(temporary(tmp), ang[iscan], $
                                 shift[0,iscan], shift[1,iscan], full=ff)
            endif else begin
              bla = red_rotation(temporary(tmp), ang[iscan], $
                                 shift[0,iscan], shift[1,iscan])
            endelse
            if(~keyword_set(nostretch)) then $
               bla = red_stretch(temporary(bla), reform(grid[iscan,*,*,*]))
            d[*,*,iwav] = rotate(temporary(bla), rot_dir) 

          endif else d[*,*,iwav] = rotate( temporary(tmp), rot_dir)

          iprogress++           ; update progress counter
          
        endfor                  ; iwav

        
        if n_elements(imean) eq 0 then begin 
          imean = fltarr(nwav)
          for ii = 0L, nwav-1 do imean[ii] = median(d[*,*,ii])
          ;;cscl = 4.0                    ; 32768 / 4096
          ;; if(keyword_set(scans_only)) then cscl = 1.0
          norm_spect = imean / 1.0 ;/ max(imean)
          norm_factor = 1.0
          spect_pos = wav *1.d10 ;+ double(prefilters[ipref])
;          print, inam + ' : saving -> '+odir + '/spectfile.'+prefilters[ipref]+'.idlsave'
          save, file = odir + '/spectfile.' + prefilters[ipref] + '.idlsave' $
                , norm_spect, norm_factor, spect_pos
        endif

        if(~keyword_set(scans_only)) then begin
          ;; Write this scan's data cube to assoc file
          if keyword_set(no_timecor) then tscl = 1 else tscl = mean(prefilter_wb) / tmean[iscan]
          if(keyword_set(float)) then begin
            dat[iscan] = d*tscl
          endif else begin
            d1 = round(d*tscl)
            dat[iscan] = fix(d1)
          endelse
          if(keyword_set(verbose)) then begin
            print, inam +'scan=',iscan,', max=', max(d1)            
          endif
        endif else begin
          ;; Write this scan's data cube as an individual file.
          print, inam + ' : saving to '+ odir + '/' + ofile
          openw, lun, odir + '/' + ofile, /get_lun
          writeu, lun, head
;          if(keyword_set(float)) then dat[iscan] = d else writeu, lun, fix(d + 0.5)
          if(keyword_set(float)) then writeu, lun, d else writeu, lun, fix(d + 0.5)
          free_lun, lun
          if keyword_set(wbwrite) then begin
            print, inam + ' : saving to '+ odir + '/' + ofilewb
            wbhead = red_mkhdr(wb) ; Just for now...
            red_writedata, odir + '/' + ofilewb, wb, head = wbhead $
                           , filetype = 'ANA', overwrite = overwrite
;              fzwrite, wb, odir + '/' + ofilewb, ' '
          endif
        endelse
      endfor                    ; iscan
      
      if(~keyword_set(scans_only)) then begin

        ;; Close assoc file for output of multi-scan data cube.
        free_lun, lun
        print, inam + ' : done'
        print, inam + ' : result saved to -> '+odir+'/'+ofile 
        if keyword_set(float) then begin
          red_flipthecube_unpol, odir+'/'+ofile, nt = Nscans, nw = Nwav $
                                 , Ny_limit = Ny_limit
        endif else begin
          red_flipthecube_unpol, odir + '/' + ofile, /icube, nt = Nscans, nw = Nwav $
                                 , Ny_limit = Ny_limit
        endelse
        ;;     make_crispex_sp_cube, odir+'/'+ofile, nwav, Nscans
      
      endif

    endfor                      ; ipref
    
  endfor                        ; itimestamp


end
