; docformat = 'rst'

;+
; Make a de-rotated and de-stretched time-series FITS data cube with
; momfbd-restored narrow-band images.
;
; It reads all temporal de-rotation and de-stretching information from
; the wideband cube produced by companion method make_wb_cube.
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
; :Returns:
; 
; 
; :Params:
; 
;     wcfile : in, type=string
; 
;       The name of the corrected WB cube file.
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
;    2017-08-17 : MGL. First version, based on code from
;                 chromis::make_crispex. 
;
;    2017-09-07 : MGL. Changed red_fitsaddpar --> red_addfitskeyword. 
; 
; 
; 
; 
;-
pro chromis::make_nb_cube, wcfile $
                           , aligncont = aligncont $
                           , clips_cont = clips_cont $
                           , momfbddir = momfbddir $
                           , no_timecor = no_timecor $
                           , nostretch = nostretch $
                           , overwrite = overwrite $
;                           , rot_dir = rot_dir $
                           , scans_only = scans_only $
                           , selscan = selscan $
                           , tiles_cont = tiles_cont $
                           , verbose = verbose $
                           , wbwrite = wbwrite 
;                           , noflats=noflats $
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Make prpara
  if n_elements(aligncont   ) ne 0 then red_make_prpara, prpara, 'aligncont'    , aligncont 
  if n_elements(blur        ) ne 0 then red_make_prpara, prpara, 'blur'         , blur         
  if n_elements(clips_cont       ) ne 0 then red_make_prpara, prpara, 'clips_cont'        , clips_cont         
  if n_elements(float       ) ne 0 then red_make_prpara, prpara, 'float'        , float  
  if n_elements(momfbddir   ) ne 0 then red_make_prpara, prpara, 'momfbddir'    , momfbddir    
  if n_elements(no_timecor  ) ne 0 then red_make_prpara, prpara, 'no_timecor'   , no_timecor 
  if n_elements(nostretch   ) ne 0 then red_make_prpara, prpara, 'nostretch'    , nostretch 
  if n_elements(np          ) ne 0 then red_make_prpara, prpara, 'np'           , np           
  if n_elements(overwrite   ) ne 0 then red_make_prpara, prpara, 'overwrite'    , overwrite
;  if n_elements(rot_dir     ) ne 0 then red_make_prpara, prpara, 'rot_dir'      , rot_dir         
  if n_elements(scans_only  ) ne 0 then red_make_prpara, prpara, 'scans_only'   , scans_only          
  if n_elements(selscan     ) ne 0 then red_make_prpara, prpara, 'selscan'      , selscan 
  if n_elements(tiles_cont       ) ne 0 then red_make_prpara, prpara, 'tiles_cont'        , tiles_cont        
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

  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0
  red_logdata, self.isodate, time_pig, pig = metadata_pig, rsun = rsun

  
  ;; Read the header from the corrected WB cube. Variables begin with
  ;; WC for Wideband Cube. 
  if ~file_test(wcfile) then begin
    print, 'WB cube missing, please run make_wb_cube.'
    print, wcfile
    retall
  endif
  wchead = red_readhead(wcfile)
  ;; Read parameters from the WB cube
  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
  fxbreadm, bunit, row = 1 $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
            , ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01
  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
  ;; wbgfiles (WideBand Global).
  fxbread, bunit, wbgfiles, 'WFILES', 1
  fxbclose, bunit

  x0 = wcX01Y01[0]
  x1 = wcX01Y01[1]
  y0 = wcX01Y01[2]
  y1 = wcX01Y01[3]
  
  self -> extractstates, wbgfiles, wbgstates
  prefilter = wbgstates[0].prefilter
  
  wchdr0 = red_readhead(wbgfiles[0])
  datestamp = strtrim(fxpar(wchdr0, 'STARTOBS'), 2)
  timestamp = (strsplit(datestamp, 'T', /extract))[1]
  
  
  
  

;
;  
;  ;; Find prefilter subdirs
;  search_dir = self.out_dir +'/'+momfbddir+'/'+timestamp+'/'
;  prefilters = file_basename(file_search(search_dir + '*' $
;                                         , count = Nprefs, /test_dir))
;  if Nprefs eq 0 then begin
;    print, inam + ' : No prefilter sub-directories found in: ' + search_dir
;    continue                    ; Next timestamp
;  endif
;  
;  ;; Select prefilter folders
;  selectionlist = strtrim(indgen(Nprefs), 2)+ '  -> ' + prefilters
;  tmp = red_select_subset(selectionlist $
;                          , qstring = inam + ' : Select prefilter directory ID:' $
;                          , count = Nprefs, indx = sindx)
;  if Nprefs eq 0 then begin
;    print, inam + ' : No prefilter sub-folders selected.'
;    continue                    ; Go to next timestamp
;  endif
;  prefilters = prefilters[sindx]
;  print, inam + ' : Selected -> '+ strjoin(prefilters, ', ')
;
;  ;; Loop over WB prefilters
;  for ipref = 0L, Nprefs-1 do begin
;

  search_dir = file_dirname(wbgfiles[0])+'/'
  extension = (strsplit(wbgfiles[0],'.',/extract))[-1]

  files = file_search(search_dir + '*.'+extension, count = Nfiles)      
  
  ;; Find all nb and wb per tuning files by excluding the global WB images 
  self -> selectfiles, files = files, states = states $
                       , cam = wbcamera, ustat = '' $
                       , sel = wbgindx, count = Nscans $
                       , complement = complement, Ncomplement = Ncomplement
  ;; We have no special state (or absence of state) to identify
  ;; the global WB images but we do know that their exposure times
  ;; are much larger than the ones corresponding to the individual
  ;; NB states.
  wbindx = where(states.exposure gt mean(states.exposure)*1.5 $
                  , Nscans, complement = complement, Ncomplement = Ncomplement)
;  wbstates = states[wbindx]
;  wbfiles = files[wbindx]

;  self -> selectfiles, files = files, states = states, cam = nbcamera, sel = sel
  
  ;; All the per-tuning files and states
  pertuningfiles = files[complement]
  pertuningstates = states[complement]

  ;; Unique tuning states
  utunindx = uniq(pertuningstates.fpi_state, sort(pertuningstates.fpi_state))
  Nwav = n_elements(utunindx)
  sortindx = sort(pertuningstates[utunindx].tun_wavelength)
  ufpi_states = pertuningstates[utunindx[sortindx]].fpi_state
  utunwavelength = pertuningstates[utunindx[sortindx]].tun_wavelength

  wav = utunwavelength
  my_prefilters = pertuningstates[utunindx[sortindx]].prefilter

  ;; Unique nb prefilters
  unbprefindx = uniq(pertuningstates[utunindx].prefilter, sort(pertuningstates[utunindx].prefilter))
  Nnbprefs = n_elements(unbprefindx)
  unbprefs = pertuningstates[utunindx[unbprefindx]].prefilter
  unbprefsref = dblarr(Nnbprefs)

  for inbpref = 0L, Nnbprefs-1 do begin
    ;; This is the reference point of the fine tuning for this prefilter:
    unbprefsref[inbpref] = double((strsplit(pertuningstates[utunindx[unbprefindx[inbpref]]].tuning $
                                            , '_', /extract))[0])
  endfor                        ; inbpref
  
  unbprefsref *= 1e-10          ; [m]
  
;  if ~keyword_set(scans_only) then begin
  ;; Look for time-series calib file
;  csearch = self.out_dir + '/calib_tseries/tseries_' + prefilters[ipref] $
;            + '_' + datestamp + '*_calib.sav'
;  cfiles = file_search(csearch, count = Ncfiles)
;  case Ncfiles of
;    0: begin
;      print, inam + ' : Could not find calibration file: ' + csearch
;      print, inam + ' : Try executing make_wb_cube on this dataset first!'
;      return
;    end
;    1: cfile = cfiles[0]
;    else: begin
;      repeat begin
;        tmp = red_select_subset(cfiles $
;                                , qstring = inam + ' : Select calibration file (scan subset).' $
;                                , count = Ncfileselect, indx = cindx, default = '-')
;      endrep until Ncfileselect eq 1
;      cfile = cfiles[cindx[0]]
;    end
;  endcase
;
;        print, inam + ' : Loading calibration file -> '+file_basename(cfile)
;        restore, cfile

  ;; Get the scan selection from wfiles (from the sav file)
;  wbgfiles = wfiles
  self -> extractstates, wbgfiles, wbgstates
  uscans = wbgstates.scannumber
  Nscans = n_elements(uscans)
;  endif else begin
;    full = 0
;    uscans = wbgstates[uniq(wbgstates.scannumber, sort(wbgstates.scannumber))].scannumber
;    Nscans = n_elements(uscans)
;    tmean = replicate(1.0, Nscans) ; Dummy time correction
;  endelse

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
  midpart = prefilter + '_' + datestamp + '_scans=' $ 
            + red_collapserange(uscans, ld = '', rd = '')

  ;; Load prefilters
  for inbpref = 0L, Nnbprefs-1 do begin
    pfile = self.out_dir + '/prefilter_fits/chromis_'+unbprefs[inbpref]+'_prefilter.idlsave'
    if ~file_test(pfile) then begin
      print, inam + ' : prefilter file not found: '+pfile
      return
    endif
      
    restore, pfile              ; Restores variable prf which is a struct
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
  endfor                        ; inbpref
    
  rpref = 1.d0/prefilter_curve[1:*]
  prefilter_wav = prefilter_wav[1:*]
  prefilter_wb = prefilter_wb[1:*]
  prefilter_curve = prefilter_curve[1:*]

  ;; Do WB correction?
  if Nwb eq Nnb then wbcor = 1B else wbcor = 0B

  ;; Load WB image and define the image border
  tmp = red_readdata(wbgfiles[0])

;  if full then begin
  Nx = wcND[0]
  Ny = wcND[1]
;  endif else begin
;    Nx = x1-x0+1
;    Ny = y1-y0+1
;  endelse
  
  ;; Create temporary cube and open output file
  d = fltarr(Nx, Ny, Nwav)  
  
  if(n_elements(odir) eq 0) then odir = self.out_dir + '/nb_cubes/' 
  file_mkdir, odir
  
;  if keyword_set(fitsoutput) then begin

  ofile = 'nb_'+midpart+'_corrected.fits'

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

  ;; Make FITS header for the NB cube
;  if ~keyword_set(scans_only) then begin
  hdr = wchead                                           ; Start with the WB cube header
  check_fits, d, hdr, /update                            ; Get the type right
  red_fitsaddkeyword, hdr, 'DATE', red_timestamp(/iso) $ ; DATE with time
                      , 'Creation UTC date of FITS header' ;
  if keyword_set(blur) then begin
    red_fitsaddkeyword, hdr, before='DATE', 'COMMENT', 'Intentionally blurred version'
  endif
  anchor = 'DATE' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'FILENAME', ofile ; New file name
  
  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'Prepare NB science data cube' $
                              , prpara = prpara $
                              , prproc = inam
  
  ;; Add info to headers
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'
  
  ;; Open fits file, dat is name of assoc variable
  dims = [Nx, Ny, Nwav, 1, Nscans] ;
  self -> fitscube_initialize, odir + ofile, hdr, lun, fileassoc, dims 
  
;  endif else begin
;    ;; Make headers for nb data as well as (possibly) for wb data
;    if keyword_set(float) then type = 4 else type = 2
;    red_mkhead, wbhead, type, [dimim[0], dimim[1]]
;    naxisx = [dimim[0], dimim[1], Nwav]
;    red_mkhead, nbhead, type, naxisx
;  endelse

  
;  endif else begin
;    
;    if keyword_set(float) then extent = '.fcube' else extent = '.icube'
;
;    if(keyword_set(scans_only)) then begin
;      head = red_unpol_lpheader(Nx, Ny, Nwav, float = float)
;    endif else begin
;      head = red_unpol_lpheader(Nx, Ny, Nwav*Nscans, float = float)
;
;      ;; Open assoc file for output of multi-scan data cube.
;      ofile = 'crispex_'+midpart+'_time-corrected'+extent
;
;      if file_test(odir + '/' + ofile) then begin
;        if keyword_set(overwrite) then begin
;          print, 'Overwriting existing data cube:'
;          print, odir + '/' + ofile
;        endif else begin
;          print, 'This data cube exists already:'
;          print, odir + '/' + ofile
;          return
;        endelse
;      endif
;      
;      openw, lun, odir + '/' + ofile, /get_lun
;      writeu, lun, head
;      point_lun, lun, 0L
;      print, inam+' assoc file -> ',  odir + '/' + file_basename(ofile,extent)+'.assoc.pro'
;      openw, lunf, odir + '/' + file_basename(ofile,extent)+'.assoc.pro', /get_lun
;      printf,lunf, 'nx=', Nx
;      printf,lunf, 'ny=', Ny
;      printf,lunf, 'nw=', Nwav
;      printf,lunf, 'nt=', Nscans
;      printf,lunf, "openr,lun,'"+ofile+"', /get_lun"
;      if keyword_set(float) then begin
;        dat = assoc(lun, fltarr(Nx, Ny, nwav, /nozero), 512)
;        printf,lunf, "dat = assoc(lun, fltarr(nx,ny,nw,/nozer), 512)"
;      endif else begin
;        dat = assoc(lun, intarr(Nx, Ny, nwav, /nozero), 512)
;        printf,lunf, "dat = assoc(lun, intarr(nx,ny,nw,/nozer), 512)"
;      endelse
;      free_lun, lunf
;    endelse
;  endelse

  ;; Set up for collecting time and wavelength data
  tbeg_array = dblarr(Nwav, Nscans)         ; Time beginning for state
  tend_array = dblarr(Nwav, Nscans)         ; Time end for state
  tavg_array = dblarr(Nwav, Nscans)         ; Time average for state
  w_array = dblarr(Nwav, Nscans)            ; Wavelengths
  hpln_array = dblarr(2, 2, Nwav, Nscans)   ; HPLN position for corners of FOV
  hplt_array = dblarr(2, 2, Nwav, Nscans)   ; HPLT position for corners of FOV
  
  ;; Start processing data
  if(~keyword_set(tiles_cont) OR (~keyword_set(clips_cont))) then begin
    tiles_cont = [8, 16, 32, 64, 128]
    clips_cont = [8, 4, 2, 1, 1]
  endif

    
  if keyword_set(aligncont) and prefilter eq '3950' then begin
    
    ;; Get shifts based on continuum vs wideband alignment.
    
    aligndir = self.out_dir + '/align/' + timestamp $
               + '/' + prefilter + '/'
    
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
;    if(iscan eq 0) then begin
;                                ;         wav = scan_nbstates.tun_wavelength
;      fzwrite, wav, odir + '/' + 'wav_' + prefilter +'.f0',' '
;    endif

    ;; The files in this scan, sorted in tuning wavelength order.
    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                         , cam = wbcamera, scan = uscans[iscan] $
                         , sel = scan_wbindx, count = count
    scan_wbfiles = pertuningfiles[scan_wbindx]
    scan_wbstates =  pertuningstates[scan_wbindx]
    sortindx = sort(scan_wbstates.tun_wavelength)
    scan_wbfiles = scan_wbfiles[sortindx]
    scan_wbstates = scan_wbstates[sortindx]
    
    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                         , cam = nbcamera, scan = uscans[iscan] $
                         , sel = scan_nbindx, count = count
    scan_nbfiles = pertuningfiles[scan_nbindx]
    scan_nbstates = pertuningstates[scan_nbindx]
    sortindx = sort(scan_nbstates.tun_wavelength)
    scan_nbfiles = scan_nbfiles[sortindx]
    scan_nbstates = scan_nbstates[sortindx]
    
    if(keyword_set(scans_only)) then begin
      if keyword_set(fitsoutput) then begin
        ofile = 'crispex_' + prefilter + '_' + datestamp + '_scan=' $
                + string(uscans[iscan], format = '(i05)') + '.fits' 
        ofilewb = 'wb_' + prefilter + '_' + datestamp + '_scan=' $
                  + string(uscans[iscan], format = '(i05)') + '.fits' 
      endif else begin
        ofile = 'crispex_' + prefilter + '_' + datestamp + '_scan=' $
                + string(uscans[iscan], format = '(i05)') + extent
        ofilewb = 'wb_' + prefilter + '_' + datestamp + '_scan=' $
                  + string(uscans[iscan], format = '(i05)') + '.fz' 
      endelse
      if file_test(odir + '/' + ofile) then begin
        if keyword_set(overwrite) then begin
          print, 'Overwriting existing data cube:'
          print, odir + '/' + ofile
        endif else begin
          print, 'Skip to next scan, this one exists already:'
          print, odir + '/' + ofile
          continue              ; Skip to next iteration of "for iscan ..." loop.
        endelse
      endif
    endif

    ;; Read global WB file to use as reference when destretching
    ;; pertuning wb files and then the corresponding nb files.
    wb = (red_readdata(wbgfiles[iscan]))[x0:x1, y0:y1]

    if keyword_set(aligncont) and prefilter eq '3950' then begin
      
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

      ;; DATE-??? keywords
      red_fitspar_getdates, nbhead $
                            , date_beg = date_beg $
                            , date_end = date_end $
                            , date_avg = date_avg
      tbeg_array[iwav, iscan] = red_time2double((strsplit(date_beg,'T',/extract))[1])
      tend_array[iwav, iscan] = red_time2double((strsplit(date_end,'T',/extract))[1])
      tavg_array[iwav, iscan] = red_time2double((strsplit(date_avg,'T',/extract))[1])
;          tbeg_array[iwav, iscan] = red_time2double((strsplit(fxpar(nbhead, 'DATE-BEG'),'T',/extract))[1])
;          tend_array[iwav, iscan] = red_time2double((strsplit(fxpar(nbhead, 'DATE-END'),'T',/extract))[1])
;          tavg_array[iwav, iscan] = red_time2double((strsplit(fxpar(nbhead, 'DATE-AVG'),'T',/extract))[1])
      ;;(tbeg_array[iwav, iscan]+tend_array[iwav, iscan])/2d

      ;; Wavelength
      w_array[iwav, iscan] = scan_nbstates[iwav].tun_wavelength

      red_wcs_hpl_coords, tavg_array[iwav, iscan], metadata_pig, time_pig $
                          , Nx, Ny, self.image_scale $
                          , hpln, hplt
      
      hpln_array[*, *, iwav, iscan] = hpln
      hplt_array[*, *, iwav, iscan] = hplt
      
;      ;; Pointing (Maybe we should average between tbeg and tend
;      ;; if there are several points in the interval? Depends on
;      ;; how noisy the positions are.)
;      hpln = interpol(metadata_pig[0, *], time_pig, tavg_array[iwav, iscan])
;      hplt = interpol(metadata_pig[1, *], time_pig, tavg_array[iwav, iscan])
;      ;; (hpln, hplt) are coordinates for the center of the FOV.
;      ;; Now tabulate the corner coordinates, assuming the FOV is
;      ;; aligned to solar coordinates. [The distance between the
;      ;; center of the FOV and the centers of the corner pixels is
;      ;; pixelsize*(Nx-1)/2 and pixelsize*(Ny-1), resp.]
;      hpln_array[0, 0, iwav, iscan] = hpln - double(self.image_scale) * (Nx-1)/2.d
;      hpln_array[1, 0, iwav, iscan] = hpln + double(self.image_scale) * (Nx-1)/2.d  
;      hpln_array[0, 1, iwav, iscan] = hpln - double(self.image_scale) * (Nx-1)/2.d 
;      hpln_array[1, 1, iwav, iscan] = hpln + double(self.image_scale) * (Nx-1)/2.d 
;      hplt_array[0, 0, iwav, iscan] = hplt - double(self.image_scale) * (Ny-1)/2.d
;      hplt_array[1, 0, iwav, iscan] = hplt - double(self.image_scale) * (Ny-1)/2.d 
;      hplt_array[0, 1, iwav, iscan] = hplt + double(self.image_scale) * (Ny-1)/2.d 
;      hplt_array[1, 1, iwav, iscan] = hplt + double(self.image_scale) * (Ny-1)/2.d 
;stop
      ;; Collect some more info for this frame here, like file
      ;; name, exposure time, detector gain, etc. Put this in FITS
      ;; tabulated keywords later.

      red_progressbar, iprogress, Nprogress $
                       , clock = clock, /predict $
                       , 'Processing scan ' $
                       + strtrim(uscans[iscan], 2) + ' state=' + state 

      ;; Get destretch to anchor camera (residual seeing)
      if(wbcor) then begin
        wwi = (red_readdata(scan_wbfiles[iwav]))[x0:x1, y0:y1]
        grid1 = red_dsgridnest(wb, wwi, tiles_cont, clips_cont)
      endif

      if 0 then wwc = red_stretch(wwi, grid1) ; test

      
      ;; Read image and apply prefilter curve
      nbim = (red_readdata(scan_nbfiles[iwav]))[x0:x1, y0:y1] * rpref[iwav]

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
      
      if keyword_set(aligncont) and prefilter eq '3950' then begin
        
        nbim = red_shift_sub(nbim, -xshifts[iwav], -yshifts[iwav])

;            ;; This is the continuum point for a Ca scan, has to be
;            ;; different for a Hb scan:
;            continnumpoint = scan_nbstates[iwav].fpi_state eq '3999_4000_+0'
;            if continnumpoint then begin
;              ;; Stretch the nb cont image to its wb image
;              gridx = red_dsgridnest(wwi, nbim, tiles_cont, clips_cont)
;              nbim = red_stretch(temporary(nbim), gridx)
;            endif
      endif

      ;; Apply destretch to anchor camera and prefilter correction
      if wbcor then nbim = red_stretch(temporary(nbim), grid1)
      
;      if(~keyword_set(scans_only)) then begin

      ;; Apply derot, align, dewarp based on the output from
      ;; polish_tseries

      nbim = red_rotation(temporary(nbim), ang[iscan], $
                          wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
      if(~keyword_set(nostretch)) then $
         nbim = red_stretch(temporary(nbim), reform(wcGRID[iscan,*,*,*]))
    
      d[*,*,iwav] = rotate(temporary(nbim), rot_dir) 

;      endif else d[*,*,iwav] = rotate( temporary(nbim), rot_dir)

      iprogress++               ; update progress counter
      
    endfor                      ; iwav

    
;    if n_elements(imean) eq 0 then begin 
;      imean = fltarr(nwav)
;      for ii = 0L, nwav-1 do imean[ii] = median(d[*,*,ii])
;      ;;cscl = 4.0                    ; 32768 / 4096
;      ;; if(keyword_set(scans_only)) then cscl = 1.0
;      norm_spect = imean / 1.0  ;/ max(imean)
;      norm_factor = 1.0
;      spect_pos = wav *1.d10 ;+ double(prefilters[ipref])
;;          print, inam + ' : saving -> '+odir + '/spectfile.'+prefilters[ipref]+'.idlsave'
;      save, file = odir + '/spectfile.' + prefilter + '.idlsave' $
;            , norm_spect, norm_factor, spect_pos
;    endif

;    if keyword_set(blur) then d = smooth(d, [29, 29, 1], /edge_wrap)

;    if(~keyword_set(scans_only)) then begin
    ;; Write this scan's data cube to assoc file
    if keyword_set(no_timecor) then tscl = 1 else tscl = mean(prefilter_wb) / wcTMEAN[iscan]
;    if keyword_set(fitsoutput) then begin
      for iwav = 0, Nwav-1 do $
         self -> fitscube_addframe, fileassoc, d[*, *, iwav] * tscl $
                                    , Nscan = Nscans, Ntuning = Nwav $
                                    , iscan = iscan, ituning = iwav 
;    endif else begin
;      if(keyword_set(float)) then begin
;        dat[iscan] = d*tscl
;      endif else begin
;        d1 = round(d*tscl)
;        dat[iscan] = fix(d1)
;      endelse
;    endelse
    if(keyword_set(verbose)) then begin
      print, inam +'scan=',iscan,', max=', max(d1)            
    endif
;     endif else begin
;       ;; Write this scan's data cube as an individual file.
;       if keyword_set(fitsoutput) then begin
;       endif else begin
;         print, inam + ' : saving to '+ odir + '/' + ofile
;         openw, lun, odir + '/' + ofile, /get_lun
;         writeu, lun, head
; ;          if(keyword_set(float)) then dat[iscan] = d else writeu, lun, fix(d + 0.5)
;         if(keyword_set(float)) then writeu, lun, d else writeu, lun, fix(d + 0.5)
;         free_lun, lun
;       endelse
;       if keyword_set(wbwrite) then begin
;         print, inam + ' : saving to '+ odir + '/' + ofilewb
;         wbhead = red_mkhdr(wb)  ; Just for now...
;         if keyword_set(fitsoutput) then begin
;           red_writedata, odir + '/' + ofilewb, wb, head = wbhead $
;                          , filetype = 'FITS', overwrite = overwrite
;         endif else begin
;           red_writedata, odir + '/' + ofilewb, wb, head = wbhead $
;                          , filetype = 'ANA', overwrite = overwrite
; ;              fzwrite, wb, odir + '/' + ofilewb, ' '
;         endelse
;       endif
;     endelse
  endfor                        ; iscan
    
;  if(~keyword_set(scans_only)) then begin

;  if keyword_set(fitsoutput) then begin

  ;; Close fits file.
  free_lun, lun
  print, inam + ' : Wrote file '+odir + ofile
  
  ;; Add any extensions.
  self -> fitscube_addwcs, odir + ofile, hpln_array, hplt_array $
                           , transpose(rebin(w_array, Nwav, Nscans, 2, 2, /sample) $
                                       , [2, 3, 0, 1]) $
                           , transpose(rebin(tavg_array, Nwav, Nscans, 2, 2, /sample) $
                                       , [2, 3, 0, 1]) 
  
  ;; Modify some headers
  value = fxpar(nbhead, 'CAMERA', comment = comment, count = count)
  if count eq 1 then fxhmodify, odir + ofile, 'CAMERA', value, comment
  value = fxpar(nbhead, 'DETECTOR', comment = comment, count = count)
  if count eq 1 then fxhmodify, odir + ofile, 'DETECTOR', value, comment
  fxhmodify, odir + ofile, 'DATE-BEG', self.isodate + 'T' + red_timestring(min(tbeg_array))
  fxhmodify, odir + ofile, 'DATE-AVG', self.isodate + 'T' + red_timestring(mean(tavg_array))
  fxhmodify, odir + ofile, 'DATE-END', self.isodate + 'T' + red_timestring(max(tend_array))

stop
  
  ;; Make a flipped cube
  
;  endif else begin
;    ;; Close assoc file for output of multi-scan data cube.
;    free_lun, lun
;    print, inam + ' : done'
;    print, inam + ' : result saved to -> '+odir+'/'+ofile 
;    if keyword_set(float) then begin
;      red_flipthecube_unpol, odir+'/'+ofile, nt = Nscans, nw = Nwav
;    endif else begin
;      red_flipthecube_unpol, odir + '/' + ofile, /icube, nt = Nscans, nw = Nwav
;    endelse
;    ;;     make_crispex_sp_cube, odir+'/'+ofile, nwav, Nscans
;  endelse
;endif

end
