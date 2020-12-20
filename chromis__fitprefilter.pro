; docformat = 'rst'

;+
; Fit prefilter model to CHROMIS data.
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
; :Keywords:
;
;    cwl : in, optional, type=double
;
;      The central wavelength (in Ångström) of the spectral line to
;      fit. Should be used only if keyword pref is also used.
;
;    dir : in, optional, type=string
;
;      Set this to the time-stamp directory to use and bypass the
;      selection dialogue. 
;
;    fit_fwhm : in, optional, type=boolean
;
;      Keep the filter profile FWHM parameter to the initial value
;      fixed during the fit.
;
;    fit_ncav : in, optional, type=boolean
;
;      Keep the number of cavities parameter fixed to the initial
;      value during the fit.
;
;    fit_parameters  : in, optional, type="boolean array[8]", default=[1,1,1,1,1,1,1,0]
;
;      What parameters to include in the fit, starting from the
;      initial value, in the following order: [scale factor, solar
;      atlas shift, prefilter shift, prefilter FWHM, prefilter number
;      of cavities, asymmetry term 1, asymmetry term 2, wavelength
;      stretch]. If nasym is used, its value will overrule whatever is
;      set for the asymmetry terms in this keyword or its default. If
;      a shorter array is given, the default will be used for the
;      unspecified items. Also the fit_xxx keywords for individual
;      parameters override this keyword.
;
;    fit_prshift : in, optional, type=boolean
;
;      Keep the filter profile shift parameter fixed to the initial
;      value during the fit.
;
;    fit_shift : in, optional, type=boolean
;
;      Keep the spectrum shift parameter fixed to the initial value
;      during the fit.
;
;    fit_stretch : in, optional, type=boolean
;
;      Keep the spectrum stretch parameter fixed to the initial value
;      during the fit.
;
;    hints : in, optional, type=boolean
;
;      If set, various hints will aid in the selection of time-stamp
;      directory. The selection dialogue will have more info than just
;      the time stamps and the FOV of the directories will be
;      displayed.
;
;    mask : in, optional, type=boolean
;
;      If set, will allow the user to mask out active areas of the FOV
;      as well as spectral positions from the fit.
;
;    nasym : in, optional, type=integer, default=0
;
;      Number of asymmetry terms to include in the fit (0, 1, or 2).
;
;    noabsunits : in, optional, type=boolean
;
;      If set, skip calibrations to establish absolute intensity
;      units.
;
;    pref : in, optional, type=string
;
;      The four-digit tag of the prefilter to calibrate.
;
;    scan : in, optional, type=integer, default=0
;
;      Use data from this single scan only.
;
;    useflats : in, optional, type=boolean
;
;      Select between flats directories rather than science data.
; 
;    value_fwhm : in, optional, type=double, default="design value"
;
;      Initial value in Ångström for the filter profile FWHM
;      parameter.
;
;    value_ncav : in, optional, type=double, default=3
;
;      Initial value for the number of cavities parameter. 
;
;    value_parameters : in, optional, type=dblarr(8)
;
;      Initial values for the fitted parameters.
;
;    value_prshift : in, optional, type=double, default=0
;
;      Initial value in Ångström for the filter profile shift
;      parameter.
;
;    value_shift : in, optional, type=double, default=0
;
;      Initial value in Ångström for the spectrum shift parameter.
;
;    value_stretch : in, optional, type=double, default=1
;
;      Initial value for the  spectrum stretch parameter.
;
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
;   2018-08-30 : MGL. Change units from 'W s^-1 m^-2 Hz^-1 sr^-1' to
;                'W m^-2 Hz^-1 sr^-1' 
; 
;   2019-12-09 : MGL. Changed correction of HRE reflectivity from
;                -0.07 to -0.09.
; 
;   2020-11-18 : MGL. New keywords cwl, nasym, pref, fit_parameters,
;                fit_fwhm, fit_shift, fit_prshift, fit_ncav,
;                fit_stretch, value_parameters, value_fwhm,
;                value_shift, value_prshift, value_ncav,
;                value_stretch. Default cwl and fwhm values
;                individually for each filter.
;
;-
pro chromis::fitprefilter, cwl = cwl_keyword $
                           , dir = dir $
                           , fit_fwhm = fit_fwhm $
                           , fit_ncav = fit_ncav $
                           , fit_parameters = fit_parameters_keyword $ 
                           , fit_prshift = fit_prshift $
                           , fit_shift = fit_shift $
                           , fit_stretch = fit_stretch $
                           , hints = hints $
                           , mask = mask $
                           , nasym = nasym $
                           , noabsunits = noabsunits $
                           , pref = pref_keyword $
                           , scan = scan $
                           , useflats = useflats $
                           , value_fwhm = value_fwhm $
                           , value_ncav = value_ncav $
                           , value_parameters = value_parameters $
                           , value_prshift = value_prshift $
                           , value_shift = value_shift $
                           , value_stretch = value_stretch

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  camWB = 'Chromis-W'
  camNB = 'Chromis-N'
  

  ;; For now! We may be able to work around this later! If there are
  ;; WB DC data throughout the observations! Can we check from here?
  noabsunits = keyword_set(useflats)

  if keyword_set(noabsunits) then begin
    units = 'dn'                ; "Digital number"
    unitscalib = 0
  endif else begin
    ;; units = 'Watt/(m2 Hz ster)' ; SI units
    units = 'W m^-2 Hz^-1 sr^-1' ; SI units
    unitscalib = 1
  endelse

  if n_elements(value_ncav) eq 0 then value_ncav = 3.d
  if n_elements(value_stretch) eq 0 then value_stretch = 1.d
  if n_elements(value_shift) eq 0 then value_shift = -0.01d
  if n_elements(value_prshift) eq 0 then value_prshift = -0.01d

  fit_parameters = [1,1,1,1,1,1,1,0] ; default
  ;; Check n_elements() rather than keyword_set() so we don't
  ;; set to false just because a keyword is not used.
  if n_elements(fit_parameters_keyword) gt 0 then begin
    fit_parameters[0] = fit_parameters_keyword ; set with keyword
  endif
  if n_elements(fit_shift)   gt 0 then fit_parameters[1] = fit_shift
  if n_elements(fit_prshift) gt 0 then fit_parameters[2] = fit_prshift
  if n_elements(fit_fwhm)    gt 0 then fit_parameters[3] = fit_fwhm
  if n_elements(fit_ncav)    gt 0 then fit_parameters[4] = fit_ncav
  if n_elements(nasym)       gt 0 then begin
    fit_parameters[5] = nasym ge 1
    fit_parameters[6] = nasym ge 2
  endif
  if n_elements(fit_stretch) gt 0 then fit_parameters[7] = fit_stretch
  
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
      indx = where(file_test(dirs+'/Chromis-N'), cnt)
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
          xds = xs/5
          yds = ys/5
          mos = fltarr(xds*Ndirs, yds)
        endif

        if keyword_set(unitscalib) && NfilesW gt 0 then begin
          im = red_readdata(fnamesW[0], /silent)
        endif else begin
          im = red_readdata(fnamesN[0], /silent)
        endelse
        if size(im, /n_dim) eq 3 then im = total(im, 3)
        ims[0, 0, idir] = im
;        if keyword_set(unitscalib) && NfilesW gt 0 then begin
;          ims[0, 0, idir] = total(red_readdata(fnamesW[0], /silent), 3)
;        endif else begin
;          ims[0, 0, idir] = total(red_readdata(fnamesN[0], /silent), 3)
;        endelse
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

;  detector = self->getdetector( camNB )
  filesNB = file_search(dirs+'/'+camNB+'/*.fits', count=nfilesNB)
  filesNB = red_sortfiles(filesNB)
  self->extractstates, filesNB, statesNB

  if keyword_set(unitscalib) then begin
    filesWB = file_search(dirs+'/'+camWB+'/*.fits', count=nfilesWB)
;    detectorwb = self->getdetector( camWB )
    filesWB = red_sortfiles(filesWB)
    self -> extractstates, filesWB, statesWB
  endif
  
  ;; Read one scan (scan 0 by default)
  
  if n_elements(scan) eq 0 then scan = 0
  idx=where(statesNB[*].scannumber eq scan, ct)
  if ct eq 0 then begin
    print, inam+' : ERROR, invalid scan number'
    return
  endif
  
  ;; Sort selected states and files
  
  idx1 = sort(statesNB[idx].tun_wavelength)
  ;;
  statesNB = statesNB[idx[idx1]]
  if keyword_set(unitscalib) then statesWB = statesWB[idx[idx1]]
  
  ustate = statesNB[uniq(statesNB[*].tun_wavelength, sort(statesNB[*].tun_wavelength))].fullstate  
  Nstates = n_elements(ustate)
  
  ;; Load data and compute mean spectrum
  time_avgs = dblarr(Nstates)
  spec      = dblarr(Nstates)
  wav       = dblarr(Nstates)
  pref      = strarr(Nstates)
  specwb    = dblarr(Nstates)

  for istate =0L, Nstates-1 do begin

    red_progressbar, istate, Nstates, 'Loop over states: '+ustate[istate], /predict

    self -> get_calib, statesNB[istate], darkdata=darkN, status = status
    if status ne 0 then begin
      print, inam+' : ERROR, cannot find dark file for NB'
      stop
    endif

    if keyword_set(unitscalib) then begin
      self -> get_calib, statesWB[istate], darkdata=darkW, status = status
      if status ne 0 then begin
        print, inam+' : ERROR, cannot find dark file for WB'
        stop
      endif
    endif

    ;; Let's not assume that all images for one state must be in the
    ;; same file... just in case.
    
    pos = where(statesNB[*].fullstate eq ustate[istate], count)
;    print, inam+'loading files for state -> '+ustate[istate]

    time_avg = dblarr(count)
    for kk=0L, count-1 do begin
      imsN = float(rdx_readdata(statesNB[pos[kk]].filename, header = hdrN))
      time_avg[kk] = red_time2double((strsplit(fxpar(hdrN, 'DATE-AVG'), 'T', /extract))[1])
      dim = size(imsN, /dim)
      if n_elements(dim) gt 2 then nsli = double(dim[2]) else nsli = 1d
      imsN = total(imsN, 3, /double) / nsli
      imsN -= darkN

      if keyword_set(unitscalib) then begin
        imsW = float(rdx_readdata(statesWB[pos[kk]].filename))
        dim = size(imsW, /dim)
        if n_elements(dim) gt 2 then nsli = double(dim[2]) else nsli = 1d
        imsW = total(imsW, 3, /double) / nsli
        imsW -= darkW
      endif
      
      if keyword_set(mask) then begin

        if kk eq 0 and istate eq 0 then begin
          mmask = red_select_area(imsN[*,*,0], /noedge, /xroi)
          nzero = where(mmask gt 0)
          bla = imsN[*,*,0]
          ind = array_indices(bla, nzero)
        endif
        
        imsN1 = double(imsN[reform(ind[0,*]),reform(ind[1,*])])
        if keyword_set(unitscalib) then imsW1 = double(imsW[reform(ind[0,*]),reform(ind[1,*])])
        
      endif else begin
        
        dx = round(dim[0]*0.12)
        dy = round(dim[1]*0.12)

        imsN1 = double(imsN[dx:dim[0]-dx-1,dy:dim[1]-dy-1])
        if keyword_set(unitscalib) then imsW1 = double(imsW[dx:dim[0]-dx-1,dy:dim[1]-dy-1])
        
      endelse ;; if mask

      spec[istate] += median(imsN1)
      if keyword_set(unitscalib) then specWB[istate] += median(imsW1)

    endfor                      ; kk

    ;; We use in spec[] the median of all frames (or average thereof,
    ;; if more than a single file), so the effective exposure time is
    ;; just 1 x XPOSURE.
    xposure = fxpar(hdrN, 'XPOSURE')
    
    spec[istate] /= count
    if keyword_set(unitscalib) then specwb[istate] /= count
    
    time_avgs[istate] = time_avg
    wav[istate]       = statesNB[pos[0]].tun_wavelength*1.d10
    pref[istate]      = statesNB[pos[0]].prefilter

  endfor                        ; istate



  upref = pref[uniq(pref, sort(pref))] 

  ;; Did we specify prefilter(s)?
  if n_elements(pref_keyword) ne 0 then begin
    if max(upref eq pref_keyword) ne 1 then begin
      print, inam + ' : Keyword pref does not match prefilters of available data.'
      help, pref_keyword, upref
      stop
      retall
    endif
    upref = pref_keyword
  endif
  Npref = n_elements(upref)

  file_mkdir, self.out_dir+'/prefilter_fits/'

  ;; Loop prefilters
  for ipref = 0L, npref-1 do begin

    red_progressbar, ipref, Npref, 'Process prefilters: ' + upref[ipref], /predict
    
    if n_elements(cwl_keyword) eq 0 then begin
      ;; Default line wavelengths (in Ångström) for different
      ;; prefilters. (Maybe we should base this on the "line" part
      ;; of the state, rather than the "prefilter" part?
      ;; (state=prefilter_line_[+-]tuning))
      case upref[ipref] of
        '3934' : cwl = 3933.66  ; Ca II K core
        '3969' : cwl = 3968.47  ; Ca II H core 
        '4862' : cwl = 4861.85  ; H-beta core
        else : cwl = double(upref[ipref])     
      endcase
    endif else cwl = cwl_keyword
      
    if n_elements(value_fwhm) eq 0 then begin
      ;; Default prefilter FWHM (in Ångström).
      case upref[ipref] of
        '3934' : value_fwhm = 4.1d ; Ca II K core
        '3969' : value_fwhm = 4.1d ; Ca II H core 
        '4862' : value_fwhm = 4.8d ; H-beta core
        else   : value_fwhm = 4.1d
      endcase
    endif

    ;; copy spectra for each prefilter
    
    idx = where(pref eq upref[ipref], nwav)
    measured_lambda = wav[idx] + (cwl - pref[idx]) ; Adjust wavelength scale
    measured_spectrum = spec[idx]
    if keyword_set(unitscalib) then wbint = mean(specwb[idx]) else wbint = 1.
    
    ;; Load satlas
    red_satlas, measured_lambda[0]-5.1, measured_lambda[-1]+5.1 $
                , atlas_lambda, atlas_spectrum $
                , /si, cont = cont 

    dw = atlas_lambda[1] - atlas_lambda[0]
    np = round((0.080 * 8) / dw)
    if np/2*2 eq np then np -=1
    tw = (dindgen(np)-np/2)*dw + double(upref[ipref])
    tr = chromis_profile(tw, erh=-0.09d0)
    tr /= total(tr)
    
    ;; Convolve the spectrum with the FPI profile
    atlas_spectrum_convolved = fftconvol(atlas_spectrum, tr)

    ;; Prepdata
    
    if Nwav gt 1 then begin

      if keyword_set(mask) then w = red_maskprefilter(measured_lambda, measured_spectrum) $
      else w = dblarr(n_elements(measured_lambda)) + 1.0d0
      dat = {xl:atlas_lambda, yl:atlas_spectrum_convolved $
             , spectrum:measured_spectrum, lambda:measured_lambda, pref:cwl, w:w}

      par = dblarr(8)
      if n_elements(value_parameters) gt 0 then par[0] = value_parameters 
      
      ;; Pars = {fts_scal, fts_shift, pref_w0, pref_dw}
      fitpars = replicate({mpside:2, limited:[0,0], limits:[0.0d, 0.0d], fixed:0, step:1.d-5}, 8)
      fitpars.fixed = ~fit_parameters
      
      ;; Scale factor
      par[0] = max(measured_spectrum) * 2d0 / cont[0]
      if ~fitpars[0].fixed then begin
        fitpars[0].limited[*] = [1,0]
        fitpars[0].limits[*]  = [0.0d0, 0.0d0]
      endif
      
      ;; Line shift (satlas-obs)
      if n_elements(value_shift) gt 0 then begin
        par[1] = value_shift
      endif
      if ~fitpars[1].fixed then begin
        fitpars[1].limited[*] = [1,1]
        fitpars[1].limits[*]  = [-1.0,1.0]
      endif
      
      ;; Pref. shift
      if n_elements(value_prshift) gt 0 then begin
        par[2] = value_prshift
      endif
      if ~fitpars[2].fixed then begin
        fitpars[2].limited[*] = [1,1]
        fitpars[2].limits[*]  = [-3.0d0,+3.0d0]
      endif
      
      ;; Pref. FWHM
      if n_elements(value_fwhm) ne 0 then begin
        par[3] = value_fwhm
      endif
      if ~fitpars[3].fixed then begin
        fitpars[3].limited[*] = [1,1]
        fitpars[3].limits[*]  = par[3] + [-0.5, 0.5] 
      endif
      
      ;; Prefilter number of cavities)
      if n_elements(value_ncav) ne 0 then begin
        par[4] = value_ncav
      endif
      if ~fitpars[4].fixed then begin
        fitpars[4].limited[*] = [1,1]
        fitpars[4].limits[*]  = [2.0d0, 3.5d0]
      endif
    
      ;; Asymmetry term 1
      if ~fitpars[5].fixed then begin
        fitpars[5].limited[*] = [1,1]
        fitpars[5].limits[*]  = [-1.d0, 1.d0]
        par[5] = 0.01d0
      endif else begin
        par[5] = 0.0d0
      endelse
      
      ;; Asymmetry term 2
      if ~fitpars[6].fixed then begin
        fitpars[6].limited[*] = [1,1]
        fitpars[6].limits[*] = [-1.d0, 1.d0]
        par[6] = 0.01d0
      endif else begin
        par[6] = 0.0d0
      endelse
      
      ;; Wavelength stretch
      if n_elements(value_stretch) ne 0 then begin
        par[7] = value_stretch
      endif
      if ~fitpars[7].fixed then begin
        fitpars[7].limited[*] = [1,1]
        fitpars[7].limits[*]  = [0.9d0, 1.1d0]
      endif
      
      ;; Now call mpfit
      par = mpfit('red_prefilterfit', par, xtol=1.e-4, functar=dat, parinfo=fitpars, ERRMSG=errmsg)
      prefilter = red_prefilter(par, dat.lambda, dat.pref) ; <------ cwl!?

      fit_fix = replicate('Fit', 8)
      for i = 0, 7 do if fitpars[i].fixed then fit_fix[i] = 'Fix'
      
      print, inam + ' : p[0] -> ', par[0], ' '+fit_fix[0]+' (scale factor)'
      print, inam + ' : p[1] -> ', par[1], ' '+fit_fix[1]+' (solar atlas shift)'
      print, inam + ' : p[2] -> ', par[2], ' '+fit_fix[2]+' (prefilter shift)'
      print, inam + ' : p[3] -> ', par[3], ' '+fit_fix[3]+' (prefilter FWHM)'
      print, inam + ' : p[4] -> ', par[4], ' '+fit_fix[4]+' (prefilter number of cavities)'
      print, inam + ' : p[5] -> ', par[5], ' '+fit_fix[5]+' (asymmetry term 1)'            
      print, inam + ' : p[6] -> ', par[6], ' '+fit_fix[6]+' (asymmetry term 2)'            
      print, inam + ' : p[7] -> ', par[7], ' '+fit_fix[7]+' (wavelength stretch)'          


      ;; save curve
      
      prf = {wav:measured_lambda $
             , pref:prefilter $
             , spec:measured_spectrum $
             , wbint:wbint $
             , reg:upref[ipref] $
             , fitpars:par $
             , fts_model:interpol(atlas_spectrum_convolved, atlas_lambda+par[1], measured_lambda)*prefilter $
             , units:units $
             , time_avg:mean(time_avgs) $
             , xposure:xposure $
            }

      cgwindow
      colors = ['blue', 'red', 'black']
      lines = [0, 2, 0]
      psyms = [16, -3, -3]
      prefilter_plot = red_prefilter(par, atlas_lambda+par[1], dat.pref) 
      mx = max([measured_spectrum $
                , atlas_spectrum_convolved*prefilter_plot $
                , prefilter_plot/par[0] * max(measured_spectrum) $
               ]) * 1.05
      
      ;; Plot measured spectrum
      cgplot, /add, measured_lambda/10., measured_spectrum, line = lines[0], color = colors[0] $
              , xtitle = '$\lambda$ / 1 nm', psym = psyms[0], yrange = [0, mx] $
              , title = camNB + ' ' + upref[ipref]
      ;; Plot atlas spectrum times prefilter profile
      cgplot,/add,/over,(atlas_lambda+par[1])/10,atlas_spectrum_convolved*prefilter_plot   $
             , color = colors[1], line = lines[1], psym = psyms[1]
      ;; Plot prefilter profile
      cgplot, /add, /over,(atlas_lambda+par[1])/10, prefilter_plot/par[0] * max(measured_spectrum) $
              , color = colors[2], line = lines[2], psym = psyms[2]
      
      cglegend, /add, align = 3, /data $
                , location = [!x.crange[0] + (!x.crange[1]-!x.crange[0])*0.1, mean(!y.crange)*.02] $
                , title = ['obs scan'], color = colors[0], psym = psyms[0], length = 0.0
      cglegend, /add, align = 5, /data, location = [mean(!x.crange), mean(!y.crange)*.02] $
                , title = ['filtered spectrum'], line = lines[1], color = colors[1], length = 0.05
      cglegend, /add, align = 2, /data $
                , location = [!x.crange[1] - (!x.crange[1]-!x.crange[0])*0.01, mean(!y.crange)*.02] $
                , title = ['fitted prefilter'], line = lines[2], color = colors[2], length = 0.05


      
      cgcontrol, output = self.out_dir + '/prefilter_fits/chromis_'+upref[ipref]+'_prefilter.pdf'
      
    endif else begin

      y1 = interpol(atlas_spectrum, atlas_lambda, measured_lambda)
      prefilter = [measured_spectrum/atlas_spectrum_convolved]
      prf = {wav:measured_lambda $
             , pref:prefilter $
             , spec:measured_spectrum $
             , wbint:wbint $
             , reg:upref[ipref] $
             , fitpars:prefilter $
             , fts_model:y1 $
             , units:units $
             , time_avg:mean(time_avgs) $
             , xposure:xposure $
            }

    endelse
    
    save, file=self.out_dir + '/prefilter_fits/chromis_'+upref[ipref]+'_prefilter.idlsave', prf
    
  endfor                        ; ipref
  
  
end
