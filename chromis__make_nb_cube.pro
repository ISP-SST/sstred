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
;    2017-09-07 : MGL. Changed red_fitsaddpar --> red_fitsaddkeyword. 
; 
;    2017-09-08 : MGL. Copy variable-keywords from the WB cube. 
; 
;    2017-09-28 : MGL. Add more variable-keywords. Make a flipped cube
;                 and copy variable keywords to it. 
; 
;    2017-10-20 : MGL. Add a WCS extension.
; 
;    2017-10-27 : MGL. New keyword noaligncont.
; 
;-
pro chromis::make_nb_cube, wcfile $
                           , noaligncont = noaligncont $
                           , clips_cont = clips_cont $
                           , momfbddir = momfbddir $
                           , no_timecor = no_timecor $
                           , nostretch = nostretch $
                           , overwrite = overwrite $
;                           , rot_dir = rot_dir $
;                           , scans_only = scans_only $
                           , selscan = selscan $
                           , tiles_cont = tiles_cont $
                           , verbose = verbose $
                           , wbwrite = wbwrite 
;                           , noflats=noflats $
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Make prpara
  if n_elements(noaligncont ) ne 0 then red_make_prpara, prpara, 'noaligncont'  , noaligncont 
  if n_elements(blur        ) ne 0 then red_make_prpara, prpara, 'blur'         , blur         
  if n_elements(clips_cont  ) ne 0 then red_make_prpara, prpara, 'clips_cont'   , clips_cont         
  if n_elements(float       ) ne 0 then red_make_prpara, prpara, 'float'        , float  
  if n_elements(momfbddir   ) ne 0 then red_make_prpara, prpara, 'momfbddir'    , momfbddir    
  if n_elements(no_timecor  ) ne 0 then red_make_prpara, prpara, 'no_timecor'   , no_timecor 
  if n_elements(nostretch   ) ne 0 then red_make_prpara, prpara, 'nostretch'    , nostretch 
  if n_elements(np          ) ne 0 then red_make_prpara, prpara, 'np'           , np           
  if n_elements(overwrite   ) ne 0 then red_make_prpara, prpara, 'overwrite'    , overwrite
;  if n_elements(rot_dir     ) ne 0 then red_make_prpara, prpara, 'rot_dir'      , rot_dir         
;  if n_elements(scans_only  ) ne 0 then red_make_prpara, prpara, 'scans_only'   , scans_only          
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

  ;; Read wcs extension of wb file to get pointing info
  fxbopen, wlun, wcfile, 'WCS-TAB', wbdr
  ttype1 = fxpar(wbdr, 'TTYPE1')
;  ttype2 = fxpar(wbdr, 'TTYPE2')
;  ttype3 = fxpar(wbdr, 'TTYPE3')
  fxbread, wlun, wwcs, ttype1
;  fxbread, wlun, whpln_index, ttype2
;  fxbread, wlun, whplt_index, ttype3
  fxbclose, wlun

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

  ofile = 'nb_'+midpart+'_corrected_im.fits'
  filename = odir+ofile
  
  if file_test(filename) then begin
    if keyword_set(overwrite) then begin
      print, 'Overwriting existing data cube:'
      print, filename
    endif else begin
      print, 'This data cube exists already:'
      print, filename
      return
    endelse
  endif

  ;; Make FITS header for the NB cube
;  if ~keyword_set(scans_only) then begin
  hdr = wchead                                             ; Start with the WB cube header
  check_fits, d, hdr, /update                              ; Get the type right
  red_fitsaddkeyword, hdr, 'DATE', red_timestamp(/iso) $   ; DATE with time
                      , 'Creation UTC date of FITS header' ;
  red_fitsdelkeyword, hdr, 'VAR_KEYS'                                ; Start fresh with variable-keywords. 

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

  ;; WB and NB data come from different cameras.
  red_fitsaddkeyword, hdr, 'CAMERA',   nbstates[0].camera
  ;; Get DETGAIN, DETOFFS, DETMODEL, DETFIRM from .fitsheader file,
  ;; i.e., red_readhead(.momfbd file). This would have to be handled
  ;; differently with CRISP because each Stokes image is a mix of two
  ;; detectors. 
  mhd=red_readhead(nbfiles[0])  ; Header of momfbd output file
  red_fitsaddkeyword, hdr, 'DETECTOR', red_fitsgetkeyword(mhd, 'DETECTOR', comment = dcomment), dcomment
  red_fitsaddkeyword, hdr, 'DETGAIN',  red_fitsgetkeyword(mhd, 'DETGAIN',  comment = dcomment), dcomment
  red_fitsaddkeyword, hdr, 'DETOFFS',  red_fitsgetkeyword(mhd, 'DETOFFS',  comment = dcomment), dcomment
  red_fitsaddkeyword, hdr, 'DETMODEL', red_fitsgetkeyword(mhd, 'DETMODEL', comment = dcomment), dcomment
  red_fitsaddkeyword, hdr, 'DETFIRM',  red_fitsgetkeyword(mhd, 'DETFIRM',  comment = dcomment), dcomment

  ;; Open fits file, dat is name of assoc variable
  dims = [Nx, Ny, Nwav, 1, Nscans] ;
  self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 
  
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
  tbeg_array = dblarr(Nwav, Nscans)     ; Time beginning for state
  tend_array = dblarr(Nwav, Nscans)     ; Time end for state
  tavg_array = dblarr(Nwav, Nscans)     ; Time average for state
  date_beg_array = strarr(Nwav, Nscans) ; DATE-BEG for state
  date_end_array = strarr(Nwav, Nscans) ; DATE-END for state
  date_avg_array = strarr(Nwav, Nscans) ; DATE-AVG for state
  exp_array = fltarr(Nwav, Nscans)      ; Total exposure time

  wcs = replicate({  wave:dblarr(2,2) $
                   , hplt:dblarr(2,2) $
                   , hpln:dblarr(2,2) $
                   , time:dblarr(2,2) $
                  }, Nwav, Nscans)

  ;; The narrowband cube is aligned to the wideband cube and all
  ;; narrowband scan positions are aligned to each other. So get hpln
  ;; and hplt from the wideband cube wcs coordinates.
  for iscan = 0L, Nscans-1 do begin
    for iwav = 0, Nwav-1 do begin
      ;; We rely here on hpln and hplt being the first two tabulated
      ;; coordinates. To make this more general, we should get the
      ;; actual indices from the headers. Maybe later...
      wcs[iwav, iscan].hpln = reform(wwcs[0,*,*,iscan])
      wcs[iwav, iscan].hplt = reform(wwcs[1,*,*,iscan])
    endfor                      ; iwav
  endfor                        ; iscan
  
  
  ;; Start processing data
  if(~keyword_set(tiles_cont) OR (~keyword_set(clips_cont))) then begin
    tiles_cont = [8, 16, 32, 64, 128]
    clips_cont = [8, 4, 2, 1, 1]
  endif

    
  if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
    
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
    
;    if(keyword_set(scans_only)) then begin
;      if keyword_set(fitsoutput) then begin
;        ofile = 'crispex_' + prefilter + '_' + datestamp + '_scan=' $
;                + string(uscans[iscan], format = '(i05)') + '.fits' 
;        ofilewb = 'wb_' + prefilter + '_' + datestamp + '_scan=' $
;                  + string(uscans[iscan], format = '(i05)') + '.fits' 
;      endif else begin
;        ofile = 'crispex_' + prefilter + '_' + datestamp + '_scan=' $
;                + string(uscans[iscan], format = '(i05)') + extent
;        ofilewb = 'wb_' + prefilter + '_' + datestamp + '_scan=' $
;                  + string(uscans[iscan], format = '(i05)') + '.fz' 
;      endelse
;      if file_test(odir + '/' + ofile) then begin
;        if keyword_set(overwrite) then begin
;          print, 'Overwriting existing data cube:'
;          print, odir + '/' + ofile
;        endif else begin
;          print, 'Skip to next scan, this one exists already:'
;          print, odir + '/' + ofile
;          continue              ; Skip to next iteration of "for iscan ..." loop.
;        endelse
;      endif
;    endif

    ;; Read global WB file to use as reference when destretching
    ;; pertuning wb files and then the corresponding nb files.
    wb = (red_readdata(wbgfiles[iscan]))[x0:x1, y0:y1]

    if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
      ;; Interpolate to get the shifts for all wavelengths for
      ;; this scan.

;      ;; Get these from somewhere else (scan_wbstates does not
;      ;; have the WB tuning or prefilter!):
;      lambdaW = 3950e-10
;      lambdaC = 3998.640d-10 + 1.258d-10
      
      icont = where(scan_nbstates.prefilter eq '3999')
      xshifts = interpol([0., nb_shifts[0, iscan]] $
;                         , [lambdaW, lambdaC]*1e7 $
                         , [scan_wbstates[icont].tun_wavelength $
                            , scan_nbstates[icont].tun_wavelength]*1e7 $
                         , scan_nbstates.tun_wavelength*1e7)
      yshifts = interpol([0., nb_shifts[1, iscan]] $
                         , [scan_wbstates[icont].tun_wavelength $
                            , scan_nbstates[icont].tun_wavelength]*1e7 $
;                         , [lambdaW, lambdaC]*1e7 $
                         , scan_nbstates.tun_wavelength*1e7)
    endif

    for iwav = 0L, Nwav - 1 do begin 

      state = ufpi_states[iwav]

      ;; Collect info about this frame here.
      
      nbhead = red_readhead(scan_nbfiles[iwav])

      ;; DATE-??? keywords
      red_fitspar_getdates, nbhead $
                            , date_beg = date_beg $
                            , date_end = date_end $
                            , date_avg = date_avg
      date_beg_array[iwav, iscan] = date_beg
      date_end_array[iwav, iscan] = date_end
      date_avg_array[iwav, iscan] = date_avg
      tbeg_array[iwav, iscan] = red_time2double((strsplit(date_beg,'T',/extract))[1])
      tend_array[iwav, iscan] = red_time2double((strsplit(date_end,'T',/extract))[1])
      tavg_array[iwav, iscan] = red_time2double((strsplit(date_avg,'T',/extract))[1])

      ;; Wavelength and time
      wcs[iwav, iscan].wave = scan_nbstates[iwav].tun_wavelength*1d9
      wcs[iwav, iscan].time = tavg_array[iwav, iscan]

      ;; Exposure time
      exp_array[iwav, iscan] = scan_nbstates[iwav].exposure
      
      red_progressbar, iprogress, Nprogress $
                       , /predict $
                       , 'Processing scan ' $
                       + strtrim(uscans[iscan], 2) + ' state=' + state 

      ;; Get destretch to anchor camera (residual seeing)
      if(wbcor) then begin
        wwi = (red_readdata(scan_wbfiles[iwav]))[x0:x1, y0:y1]
        grid1 = red_dsgridnest(wb, wwi, tiles_cont, clips_cont)
      endif

      ;; Read image and apply prefilter curve
      nbim = (red_readdata(scan_nbfiles[iwav]))[x0:x1, y0:y1] * rpref[iwav]

      if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
        ;; Apply alignment to compensate for time-variable chromatic
        ;; aberrations.
        nbim = red_shift_sub(nbim, -xshifts[iwav], -yshifts[iwav])
      endif

      ;; Apply destretch to anchor camera and prefilter correction
      if wbcor then nbim = red_stretch(temporary(nbim), grid1)

      ;; Apply derot, align, dewarp based on the output from
      ;; polish_tseries

      nbim = red_rotation(temporary(nbim), ang[iscan], $
                          wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
      if(~keyword_set(nostretch)) then $
         nbim = red_stretch(temporary(nbim), reform(wcGRID[iscan,*,*,*]))
    
      d[*,*,iwav] = rotate(temporary(nbim), rot_dir) 

      iprogress++               ; update progress counter
      
    endfor                      ; iwav
    
    if keyword_set(no_timecor) then tscl = 1 else tscl = mean(prefilter_wb) / wcTMEAN[iscan]
      for iwav = 0, Nwav-1 do $
         self -> fitscube_addframe, fileassoc, d[*, *, iwav] * tscl $
                                    , Nscan = Nscans, Ntuning = Nwav $
                                    , iscan = iscan, ituning = iwav 
    if(keyword_set(verbose)) then begin
      print, inam +'scan=',iscan,', max=', max(d1)            
    endif
  endfor                        ; iscan

  ;; Close fits file and make a flipped version.
  self -> fitscube_finish, lun, flipfile = flipfile, wcs = wcs

  self -> fitscube_addvarkeyword, filename, 'DATE-BEG', date_beg_array $
                                  , comment = 'Beginning of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(min(tbeg_array)) $
                                  , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, filename, 'DATE-END', date_end_array $
                                  , comment = 'End time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(max(tend_array)) $
                                  , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, filename, 'DATE-AVG', date_avg_array $
                                  , comment = 'Average time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(mean(tavg_array)) $
                                  , axis_numbers = [3, 5] 

  ;; Modify some headers (should do in flipfile as well) Or should
  ;; they be var-keys?
;  fxhmodify, filename, 'DATE-BEG', self.isodate + 'T' + red_timestring(min(tbeg_array))
;  fxhmodify, filename, 'DATE-AVG', self.isodate + 'T' + red_timestring(mean(tavg_array))
;  fxhmodify, filename, 'DATE-END', self.isodate + 'T' + red_timestring(max(tend_array))

  ;; Copy variable-keywords from wb cube file.
  self -> fitscube_addvarkeyword, filename, 'SCANNUM',  old_filename = wcfile
  self -> fitscube_addvarkeyword, filename, 'ATMOS_R0', old_filename = wcfile

  ;; Add also XPOSURE but based on NB data
  self -> fitscube_addvarkeyword, filename $
                                  , 'XPOSURE', comment = 'Summed exposure times' $
                                  , tunit = 'ms' $
                                  , exp_array, keyword_value = mean(exp_array) $
                                  , axis_numbers = [3, 5] 


  ;; Copy some variable-keywords from the ordinary nb cube to the
  ;; flipped version.
  self -> fitscube_addvarkeyword, flipfile, 'SCANNUM',  old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'ATMOS_R0', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-BEG', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-AVG', old_filename = filename, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-END', old_filename = filename, /flipped

  if 0 then begin

    im=readfits(filename)  
    sp=readfits(flipfile)

    him=headfits(filename)  
    hsp=headfits(flipfile)

    scn = red_fitsgetkeyword(filename, 'SCANNUM', comment = comment, variable_values = scn_values)
    xps = red_fitsgetkeyword(filename, 'XPOSURE', comment = comment, variable_values = xps_values)
    print, scn, xps

    r0 = red_fitsgetkeyword(filename, 'ATMOS_R0', comment = comment, variable_values = r0_values)
    help, r0_values
  
    
  endif
  
end
