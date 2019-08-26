; docformat = 'rst'

;+
; Fit prefilter model to CRISP data.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz, ISP
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
;    shift : in, optional, type=float
;
;        Initial value for the line vs. atlas shift.
;
;    init : in, optional, type=array
;
;        Initialize the parameter array to this.
;
;    fixcav : in, optional, type=float, default='Not fixed'
;
;        Fix the number of cavities to this number.
;
;    cwl : in, optional, type=float, default="Varies depending on prefilter"
;
;        The central wavelength of the line in Ångström. SHould be
;        used only if keyword pref is also used.
;
;    dir : in, optional, type=string
;
;        Set this to the time-stamp directory to use and bypass the
;        selection dialogue. 
;
;    hints : in, optional, type=boolean
;
;        If set, various hints will aid in the selection of time-stamp
;        directory. The selection dialogue will have more info than
;        just the time stamps and the FOV of the directories will be
;        displayed.
;
;    mask : in, optional, type=boolean
;
;        If set, will allow the user to mask out spectral positions
;        from the fit.
;
;    nasym : in, optional, type=integer, default=0
;
;        The number of asymmetry terms to fit (0, 1, or 2).
;
;    noabsunits : in, optional, type=boolean
;
;        If set, skip calibrations to establish absolute intensity
;        units. 
;
;    pref : in, optional, type=string
;
;        Select prefilter to fit.
;
;    scan : in, optional, type=integer, default=0
;
;        Use data from this single scan only.
;
;    stretch : in, optional, type=boolean
;
;        Allow the wavelength scale to stretch.
;
;    useflats : in, optional, type=boolean
;
;        Select between flats directories rather than science data. 
; 
; :History:
; 
;   2016-11-28 : MGL. Moved helper functions to their own files and
;                added header. Make and save a final plot of the fit.
;                Prevent user from setting both /cgs and /si keywords.
; 
;   2016-12-04 : JdlCR. Allow to mask regions of the mean spectra.
;                Sometimes there is no real quiet-sun and the line
;                center must be masked.
;
;   2017-02-14 : JdlCR. Allow to also mask a section of the FOV. Many
;                observers forget to take quiet-Sun data for
;                calibration.
;
;   2017-04-07 : MGL. Use XROI GUI to select area. Added progress
;                bars. 
;
;   2017-04-18 : MGL. Remove si and cgs keywords, always use SI units.  
;
;   2017-06-05 : MGL. Construct the units string as specified in the
;                FITS standard 3.0.
;
;   2017-07-06 : THI. Use rdx_readdata to also support compressed data.
;
;   2017-10-05 : MGL. Selection list now includes mu and number of
;                files in directories and has a default selection
;                based on those numbers. Also display mosaic of FOV
;                for the different directories. New keyword dir.
;
;   2017-11-28 : MGL. Add legends and color to diagnostic plot.
;
;   2017-12-01 : MGL. New keyword hints. Reorganize hints
;                calculations. New keyword useflats.
;
;   2017-12-04 : MGL. New keyword noabsunits.
;
;   2018-11-30 : MGL. New version based on chromis::fitprefilter.
;
;   2018-12-17 : MGL. New keywords pref and cwl.
;
;   2018-12-17 : MGL. New keywords shift, init, fixcav, nasym,
;                stretch. 
;
;   2019-06-14 : MGL. Do backscatter correction for 8542 and 7772
;                (this actually changes the measured intensity) and
;                gain correction of images.
; 
;-
pro crisp::fitprefilter, cwl = cwl $
                         , dir = dir $
                         , hints = hints $
                         , fixcav = fixcav $
                         , init = init $
                         , mask = mask $
                         , nasym = nasym $
                         , noabsunits = noabsunits $
                         , pref = pref_keyword $
                         , scan = scan $
                         , shift = shift $
                         , stretch = stretch $
                         , useflats = useflats 
;                           , time = time 
  
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])
  
  camWB = 'Crisp-W'
  camsNB = ['Crisp-T', 'Crisp-R']
  camNB = camsNB[0]
  
  ;; For now! We may be able to work around this later!
  noabsunits = keyword_set(useflats)
  
  if keyword_set(noabsunits) then begin
    units = 'dn'                ; "Digital number"
    unitscalib = 0
  endif else begin
    ;; units = 'Watt/(m2 Hz ster)' ; SI units
    units = 'W m^-2 Hz^-1 sr^-1' ; SI units
    unitscalib = 1
  endelse

  if n_elements(nasym) eq 0 then nasym = 0
  
  if n_elements(dir) eq 0 then begin

    ;; Directory not provided, user has to choose one.
    
    if keyword_set(useflats) then begin

      if ~ptr_valid(self.flat_dir) then begin
        print, inam+' : ERROR : undefined flat_dir'
        return
      endif
      dirs = *self.flat_dir

      ;; In flats directories WB and NB data are often separated. So
      ;; use only directories where there actually are NB data.
      indx = where(file_test(dirs+'/'+camsNB[0]), cnt)
      if cnt eq 0 then begin
        print, inam + ' : No flats directories with NB data.'
        print, dirs
        retall
      endif
      dirs = dirs[indx]

    endif else begin
      
      if ~ptr_valid(self.data_dirs) then begin
        print, inam+' : ERROR : undefined data_dir'
        return
      endif
      dirs = *self.data_dirs

    endelse
    Ndirs = n_elements(dirs)
    if( Ndirs eq 0) then begin
      print, inam+' : ERROR : no directories defined'
      return
    endif else begin
      if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
      else dirstr = dirs[0]
    endelse

    if keyword_set(hints) then begin

      ;; Find some info about the directories
      
      prefs = strarr(Ndirs)
      fnames = strarr(Ndirs)
      times = dblarr(Ndirs)
      Nfiles = lonarr(Ndirs)
      contr = fltarr(Ndirs)

      ;; First get mu and zenith angle
      timeregex = '[0-2][0-9]:[0-5][0-9]:[0-6][0-9]'
      for idir = 0, Ndirs-1 do begin
        times[idir] = red_time2double(stregex(dirs[idir], timeregex, /extract))
      endfor                    ; idir
      red_logdata, self.isodate, times, mu = mu, zenithangle = za

      
      for idir = 0, Ndirs-1 do begin

        print, dirs[idir]
        fnamesN = file_search(dirs[idir]+'/'+camNB+'/*',   count = NfilesN)
        Nfiles[idir] = NfilesN
        if keyword_set(unitscalib) then begin
          fnamesW = file_search(dirs[idir]+'/'+camWB+'/*', count = NfilesW)
        endif else begin
          NfilesW = 0
        endelse

        if n_elements(ims) eq 0 then begin
          hdr = red_readhead(fnamesN[0])
          xs = red_fitsgetkeyword(hdr, 'NAXIS1')
          ys = red_fitsgetkeyword(hdr, 'NAXIS2')
          ims = fltarr(xs, ys, Ndirs)
          xds = xs/4
          yds = ys/4
          mos = fltarr(xds*Ndirs, yds)
        endif

        if keyword_set(unitscalib) && NfilesW gt 0 then begin
          im = red_readdata(fnamesW[0], /silent)
          ;;ims[0, 0, idir] = total(red_readdata(fnamesW[0], /silent), 3)
        endif else begin
          im = red_readdata(fnamesN[0], /silent)
          ;;ims[0, 0, idir] = total(red_readdata(fnamesN[0], /silent), 3)
        endelse
        if size(im, /n_dim) eq 3 then im = total(im, 3)
        ims[0, 0, idir] = im
        mos[idir*xds:(idir+1)*xds-1, *] $
           = rebin(ims[*, *, idir], xds, yds)/median(ims[*, *, idir])
        ;; Get more hints only for potentially interesting dirs
        if NfilesN gt 0 && NfilesN lt 1000 && mu[idir] gt 0.9 then begin

          self -> extractstates, fnamesN, sts
          prefs[idir] = ', prefs='+strjoin(sts[uniq(sts.prefilter,sort(sts.prefilter))].prefilter, ',')

          contr[idir] = stddev(ims[20:-20, 20:-20, idir])/mean(ims[20:-20, 20:-20, idir])

;          red_fitspar_getdates, hdr, date_beg = date_beg 
;          times[idir] = red_time2double((strsplit(date_beg, 'T', /extract))[1])

        endif
      endfor                    ; idir

      ;; Select data folder containing a quiet-Sun-disk-center dataset

      ;; Default based on some heuristics of unknown value. Deselect
      ;; directories with many files, far from disc center, with large
      ;; contrast (possibly spots?).
      findx = where(Nfiles gt 0 and Nfiles lt 400 $
                    and mu gt 0.9 $
                    and contr lt median(contr) $
                    , Nwhere)
;      findx = where(Nfiles lt 500 and mu gt 0.9, Nwhere)

      case Nwhere of
        0 : tmp = min(za, default) ; Default default is close to local noon
        1 : default = findx[0]
        else : begin            ; Have to pick one of findx
          tmp = min(za[findx], ml) 
          default = findx[ml]
        end
      endcase
      
      selectionlist = file_basename(dirs) $
                      + ' (µ=' +string(mu,format='(f4.2)') $
                      + ', #files='+strtrim(Nfiles,2) $
                      + prefs $
                      + ')'
      
      hintlist = file_basename(dirs) $
                 + ' (µ=' +string(mu,format='(f4.2)') $
                 + ', #files='+strtrim(Nfiles,2) $
                 + ')'
      

      ;; Display some visual hints
      scrollwindow, xs = xds*Ndirs, ys = yds
      tv, bytscl(mos, .6, 1.4)
      for idir = 0, Ndirs-1 do $
         cgtext, 5+idir*xds, 5, align = 0, /device, color = default eq idir?'cyan':'green' $
                 , strtrim(idir, 2)+' : '+red_strreplace(hintlist[idir], 'µ', '$\mu$')
      
      qstring = 'Select data to be processed'

    endif else begin

      ;; No hints
      selectionlist = file_basename(dirs)
      qstring = 'Select data to be processed (call with /hints to get more info or with dir=(timestamp) to bypass)'
      default = 0
      
    endelse

    
    ;; Do the selection
    tmp = red_select_subset(selectionlist $
                            , qstring = qstring $
                            , count = Nselect, indx = sindx, default = default)
    idx = sindx[0]
    print, inam + ' : Will process the following data:'
    print, selectionlist[idx], format = '(a0)'
    dirs = dirs[idx]
    
  endif else begin
    
    if file_test(dir) then begin
      dirs = dir
    endif else begin
      ;; Maybe just the timestamp dir?
      td = file_dirname((*self.data_dirs)[0])
      dirs = td + '/' + dir
    endelse
    
  endelse
  
  ;; Get files and states

  ;; Get one scan (scan 0 by default)
  if n_elements(scan) eq 0 then scan = 0

  for icam = 0, 1 do begin

    camNB = camsNB[icam]
    
    filesNBall = file_search(dirs+'/'+camNB+'/*', count=nfilesNB)
    self->extractstates, filesNBall, statesNBall
    
    idx=where(statesNBall.scannumber eq scan, ct)
    if ct eq 0 then begin
      print, inam+' : ERROR, invalid scan number'
      return
    endif
    statesNBall = statesNBall[idx]
    filesNBall = filesNBall[idx]
    
    if keyword_set(unitscalib) then begin ; Get the corresponding WB files
      filesWBall = file_search(dirs+'/'+camWB+'/*', count=nfilesWB)
      self -> extractstates, filesWBall, statesWBall

      idx=where(statesWBall.scannumber eq scan, ct)
      if ct eq 0 then begin
        print, inam+' : ERROR, invalid scan number'
        return
      endif
      statesWBall = statesWBall[idx]
      filesWBall = filesWBall[idx]
      
      ;; NEED NB and WB file lists to be synched!
      if min(statesWBall.framenumber eq statesNBall.framenumber) eq 0 then stop ; Will this happen?
      
    endif

    upref = statesNBall[uniq(statesNBall.prefilter,sort(statesNBall.prefilter))].prefilter

    ;; Did we specify prefilter(s)?
    if n_elements(pref_keyword) ne 0 then begin
      if max(upref eq pref_keyword) ne 1 then begin
        print, inam + ' : Keyword pref does not match prefilters of available data.'
        help, pref_keyword, upref
        retall
      endif
      upref = pref_keyword
    endif
    Npref = n_elements(upref)
    
    ;; Loop prefilters
    for ipref = 0L, Npref-1 do begin

      print, inam + ' : Process prefilter ' + upref[ipref]

      if n_elements(cwl) eq 0 then begin
        ;; Default line wavelengths (in Ångström) for different
        ;; prefilters. (Maybe we should base this on the "line" part
        ;; of the state, rather than the "prefilter" part?
        ;; (state=prefilter_line_[+-]tuning))
        case upref[ipref] of
          '5173' : cwl = 5172.70   
          '5876' : cwl = 5876.28   
          '5896' : cwl = 5895.93
          '6173' : cwl = 6173.34
          '6302' : cwl = 6302.50
          '6563' : cwl = 6562.82
          '7772' : cwl = 7772.00
          '8542' : cwl = 8542.13
          else : cwl = double(upref[ipref])
        endcase
      endif
      
      self -> selectfiles, files = filesNBall, states = statesNBall $
                           , prefilter = upref[ipref] $
                           , selected = sel, count = Nsel
      
      statesNB = statesNBall[sel]
      
      ;; Sort selected states and files
      idx1 = sort(statesNB.tun_wavelength)
      statesNB = statesNB[idx1]
      
      if keyword_set(unitscalib) then begin
        ;; Select the corresponding WB frames
        statesWB = statesWBall[sel[idx1]]
      endif
      
      ustates = statesNB[uniq(statesNB.tun_wavelength, sort(statesNB.tun_wavelength))]
      Nwav = n_elements(ustates)

      ;; Load data and compute mean spectrum
      time_avgs = dblarr(Nwav)
      spec      = dblarr(Nwav)
      wav       = dblarr(Nwav)
      pref      = strarr(Nwav)
      specwb    = dblarr(Nwav)

      if self.dodescatter and (statesNB[0].prefilter eq '8542' $
                               or statesNB[0].prefilter eq '7772') then begin
        self -> loadbackscatter, statesNB[0].detector $
                                 , statesNB[0].prefilter, bgainn, bpsfn
        if keyword_set(unitscalib) then begin
          self -> loadbackscatter, statesWB[0].detector $
                                 , statesWB[0].prefilter, bgainw, bpsfw
        endif
      endif 
      
      for istate = 0L, Nwav-1 do begin

        red_progressbar, istate, Nwav, 'Process fpi_state '+ustates[istate].fpi_state, /predict

        ;; Let's not assume that all images for one state must be in the
        ;; same file... just in case.
        
        pos = where(statesNB.fpi_state eq ustates[istate].fpi_state, count)
        
        ;; Get darks for this camera
        self -> get_calib, statesNB[pos[0]], darkdata=darkN, gaindata=gainN, status = status
        if status ne 0 then begin
          print, inam+' : ERROR, cannot find dark file for NB'
          stop
        endif
        dim = size(darkN, /dim)
        
        if keyword_set(unitscalib) then begin
          self -> get_calib, statesWB[pos[0]], darkdata=darkW, gaindata=gainW, status = status
          if status ne 0 then begin
            print, inam+' : ERROR, cannot find dark file for WB'
            stop
          endif
        endif
        
        ;; Sum files with same tuning, checking for outliers. The returned
        ;; "sum" is the average of the summed frames, so still in counts
        ;; for a single frame. Also correct for dark.
        imN = rdx_sumfiles(statesNB[pos].filename, /check, nthreads = 4, time_avg = time_avg) - darkN
        if self.dodescatter and (statesNB[pos[0]].prefilter eq '8542' $
                                 or statesNB[pos[0]].prefilter eq '7772') then begin
          imN = rdx_descatter(temporary(imN), bgainn, bpsfn, nthreads = nthread)
        endif
        imN *= rdx_fillpix(gainN)

        if keyword_set(unitscalib) then begin
          imW = rdx_sumfiles(statesWB[pos].filename, /check, nthreads = 4) - darkW
          if self.dodescatter and (statesWB[pos[0]].prefilter eq '8542' $
                                   or statesWB[pos[0]].prefilter eq '7772') then begin
            imW = rdx_descatter(temporary(imW), bgainw, bpsfw, nthreads = nthread)
          endif
          imW *= rdx_fillpix(gainW)
        endif
        
        ;; Get the spectrum point in counts
        if keyword_set(mask) then begin
          if istate eq 0 then begin
            mmask = red_select_area(imN, /noedge, /xroi)
            ind = where(mmask gt 0)
          endif 
          spec[istate] = median(double(imN[ind]))
          if keyword_set(unitscalib) then $
             specWB[istate] = median(double(imW[ind]))
        endif else begin
          dx = round(dim[0]*0.12)
          dy = round(dim[1]*0.12)
          spec[istate] = median(double(imN[dx:dim[0]-dx-1,dy:dim[1]-dy-1]))
          if keyword_set(unitscalib) then $
             specWB[istate] = median(double(imW[dx:dim[0]-dx-1,dy:dim[1]-dy-1]))
        endelse

        time_avgs[istate] = red_time2double(time_avg)
        wav[istate]       = statesNB[pos[0]].tun_wavelength*1.d10 ; [Å] Tuning wavelength
        pref[istate]      = statesNB[pos[0]].prefilter
        
        wav[istate]       += cwl - pref[istate] ; Adjust wavelength scale
        
      endfor                    ; istate

      ;; We use in spec[] rdx_sumfiles() of the data files, so the
      ;; effective exposure time is just 1 x XPOSURE.
      hdrN = red_readhead(statesNB[pos[0]].filename)
      xposure = fxpar(hdrN, 'XPOSURE')

      ;; copy spectra for each prefilter
      
      idx = where(pref eq upref[ipref], nwav)
      lambda = wav[idx]
      spectrum = spec[idx]
      if keyword_set(unitscalib) then wbint = mean(specwb) else wbint = 1.
      
      ;; Load satlas
      ;; red_satlas, wav[0]-5.1, wav[-1]+5.1, xl, yl, /si, cont = cont 
      red_satlas, wav[0]-1, wav[-1]+1, xl, yl, /si, cont = cont 
      
      ;; Make FPI transmission profile
      dw = xl[1] - xl[0]
      np = round((0.080 * 8) / dw)
      np = long((max(xl) - min(xl)) / dw) - 2
      if np/2*2 eq np then np -=1
      tw = (dindgen(np)-np/2)*dw                                               ;+ double(upref[ipref])
      tr = self -> fpi_profile(tw, upref[ipref], erh=-0.01d, /offset_correction) ;
      tr /= total(tr)
      
      ;; Convolve the spectrum with the FPI profile
      yl1 = fftconvol(yl, tr)     
      
      ;; Prepdata
      
      if keyword_set(mask) then w = red_maskprefilter(wav, spec) else w = dblarr(n_elements(wav)) + 1.0d0
      
      dat = {xl:xl, yl:yl1, spectrum:spectrum, lambda:lambda, pref:cwl, w:w}

      ;; Init guess model
      par = dblarr(8)
      ;; par = [max(spectrum) * 2d0 / cont[0], 0.01d0, 0.01d0, 3.3d0, 3.0d0, -0.01d0, -0.01d0, 1d]
      if n_elements(init) gt 0 then par[0] = init else begin
        par[0] =  max(spectrum) * 2d0 / cont[0] ; Scale factor
        par[1] = -0.001d0                       ; Line shift (satlas-obs)
        par[2] = 0d                             ; Pref. shift (-0.5d in master?)
        par[3] = 6d                             ; Pref. FWHM
        par[4] = 2d                             ; Prefilter number of cavities)
        par[5] = 0d                             ; asymmetry term 1
        par[6] = 0d                             ; asymmetry term 2
        par[7] = 1d                             ; Wavelength stretch
      endelse

      ;; Pars = {fts_scal, fts_shift, pref_w0, pref_dw}
      fitpars = replicate({mpside:2, limited:[0,0], limits:[0.0d, 0.0d], fixed:0, step:1.d-5}, 8)

      ;; Scale factor
      fitpars[0].limited[*] = [1,0]
      fitpars[0].limits[*] = [0.d0, 0.0d0]

      ;; Line shift (satlas-obs)
      if n_elements(shift) gt 0 then begin
        par[1] = shift
      endif else begin
        fitpars[1].limited[*] = [1,1]
        fitpars[1].limits[*] = [-1.0,1.0]
      endelse
      
      ;; Pref. shift
      fitpars[2].limited[*] = [1,1]
      fitpars[2].limits[*] = [-3.0d0,+3.0d0]

      ;; Pref. FWHM
      fitpars[3].limited[*] = [1,1]
      case upref[ipref] of
        '8542' : begin
          ;; 8542 is much wider than the other CRISP prefilters
          par[3] = 8.5d
          fitpars[3].limits[*] = [8d, 9d]          
        end
        else : begin
          par[3] = 4.d
          fitpars[3].limits[*] = [2d, 6d]
        end
      endcase
      
      ;; Prefilter number of cavities)
      if n_elements(fixcav) ne 0 then begin
        par[4] = fixcav
        fitpars[4].fixed = 1
      endif else begin
        fitpars[4].limited[*] = [1,1]
        fitpars[4].limits[*] = [1.4d0, 2.4d0]
        fitpars[4].fixed = 1
      endelse

      ;; Asymmetry term 1
      if nasym gt 0 then begin
        fitpars[5].limited[*] = [1,1]
        fitpars[5].limits[*] = [-1.d0, 1.d0]
        par[5] = 0.0001d0       ; asymmetry term 1
      endif else begin
        fitpars[5].fixed = 1
      endelse

      ;; Asymmetry term 2
      if nasym gt 1 then begin
        fitpars[6].limited[*] = [1,1]
        fitpars[6].limits[*] = [-1.d0, 1.d0]
        par[6] = 0.0001d0       ; asymmetry term 2
      endif else begin
        fitpars[6].fixed = 1
      endelse

      ;; Wavelength stretch
      if ~keyword_set(stretch) then begin
        fitpars[7].fixed = 1
      end
      
      ;; Now call mpfit
      par = mpfit('red_prefilterfit', par, xtol=1.e-4, functar=dat, parinfo=fitpars, ERRMSG=errmsg) ; <------
      prefilter = red_prefilter(par, dat.lambda, dat.pref)                                          ; <------

      print, inam + ' : p[0] -> ', par[0], ' (scale factor)'
      print, inam + ' : p[1] -> ', par[1], ' (solar atlas shift)'
      print, inam + ' : p[2] -> ', par[2], ' (prefilter shift)'
      print, inam + ' : p[3] -> ', par[3], ' (prefilter FWHM)'
      print, inam + ' : p[4] -> ', par[4], ' (prefilter number of cavities)'
      print, inam + ' : p[5] -> ', par[5], ' (asymmetry term 1)'            
      print, inam + ' : p[6] -> ', par[6], ' (asymmetry term 2)'            
      print, inam + ' : p[7] -> ', par[7], ' (Wavelength stretch)'          



      
      ;; save curve
      file_mkdir, self.out_dir+'/prefilter_fits/'
      prf = {wav:lambda $
             , pref:prefilter $
             , spec:spectrum $
             , wbint:wbint $
             , reg:upref[ipref]$
             , fitpars:par $
             , fts_model:interpol(yl1, xl+par[1], wav)*prefilter $
             , units:units $
             , time_avg:mean(time_avgs) $
             , xposure:xposure $
            }

      ;; Save the fit
      save, prf $
            , file = self.out_dir + '/prefilter_fits/' $
            + camNB + '_' + upref[ipref] + '_prefilter.idlsave'

      cgwindow
      mx = max([spectrum, interpol(yl1, xl+par[1], wav)*prefilter, prefilter/par[0] * max(spectrum)]) * 1.05
      
      colors = ['blue', 'red', 'black']
      lines = [0, 2, 0]
      psyms = [16, -3, -3]
      prefilter_plot = red_prefilter(par, xl, dat.pref) ; <------
      cgplot, /add, lambda/10., spectrum, line = lines[0], color = colors[0] $
              , xtitle = '$\lambda$ / 1 nm', psym = psyms[0], yrange = [0, mx] $
              , title = camNB + ' ' + upref[ipref]
      cgplot, /add,/over,(xl+par[1])/10,yl1*prefilter_plot   $
              , color = colors[1], line = lines[1], psym = psyms[1]
      cgplot, /add, /over,(xl+par[1])/10, prefilter_plot/par[0] * max(spectrum) $
              , color = colors[2], line = lines[2], psym = psyms[2]
      
      
      
;    cgplot, /add, /over, wav/10., interpol(yl1, xl+par[1], wav)*prefilter $
;    , color = colors[1], line = lines[1], psym = psyms[1]
;    cgplot, /add, /over, wav/10., prefilter/par[0] * max(spectrum), color = colors[2], line = lines[2], psym = psyms[2]
      
      cglegend, /add, align = 3, /data $
                , location = [!x.crange[0] + (!x.crange[1]-!x.crange[0])*0.1, mean(!y.crange)*.02] $
                , title = ['obs scan'], color = colors[0], psym = psyms[0], length = 0.0
      cglegend, /add, align = 5, /data, location = [mean(!x.crange), mean(!y.crange)*.02] $
                , title = ['filtered spectrum'], line = lines[1], color = colors[1], length = 0.05
      cglegend, /add, align = 2, /data $
                , location = [!x.crange[1] - (!x.crange[1]-!x.crange[0])*0.01, mean(!y.crange)*.02] $
                , title = ['fitted prefilter'], line = lines[2], color = colors[2], length = 0.05
      
;      file_mkdir, self.out_dir+'/prefilter_fits/'
      
      cgcontrol, output = self.out_dir + '/prefilter_fits/'+camNB+'_'+upref[ipref]+'_prefilter.pdf'
      
    endfor                      ; ipref
    
  endfor                        ; icam

end
