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
;    noabsunits : in, optional, type=boolean
;
;        If set, skip calibrations to establish absolute intensity
;        units. 
;
;    scan : in, optional, type=integer, default=0
;
;        Use data from this single scan only.
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
; 
;-
pro chromis::fitprefilter, dir = dir $
                           , useflats = useflats $
                           , noabsunits = noabsunits $
                           , hints = hints $
                           , mask = mask $
;                           , pref = pref $
                           , scan = scan ;$
;                           , time = time 
  

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; For now! We may be able to work around this later!
  noabsunits = keyword_set(useflats)
  
  if keyword_set(noabsunits) then begin
    units = 'dn'                ; "Digital number"
    unitscalib = 0
  endif else begin
    ;; units = 'Watt/(s m2 Hz ster)' ; SI units
    units = 'W s^-1 m^-2 Hz^-1 sr^-1' ; SI units
    unitscalib = 1
  endelse
  
  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

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
        
        fnamesN = file_search(dirs[idir]+'/Chromis-N/*fits', count = NfilesN)
        Nfiles[idir] = NfilesN
        if keyword_set(unitscalib) then begin
          fnamesW = file_search(dirs[idir]+'/Chromis-W/*fits', count = NfilesW)
        endif else begin
          NfilesW = 0
        endelse

        if n_elements(ims) eq 0 then begin
          hdr = red_readhead(fnamesN[0])
          xs = red_fitsgetkeyword(hdr, 'NAXIS1')
          ys = red_fitsgetkeyword(hdr, 'NAXIS2')
          ims = fltarr(xs, ys, Ndirs)
          mos = fltarr(xs/5*Ndirs, ys/5)
        endif

        if keyword_set(unitscalib) && NfilesW gt 0 then begin
          ims[0, 0, idir] = total(red_readdata(fnamesW[0], /silent), 3)
        endif else begin
          ims[0, 0, idir] = total(red_readdata(fnamesN[0], /silent), 3)
        endelse
        mos[idir*xs/5:(idir+1)*xs/5-1, *] $
           = rebin(ims[*, *, idir], xs/5, ys/5)/median(ims[*, *, idir])
        
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
      scrollwindow, xs = xs/5*Ndirs, ys = ys/5
      tv, bytscl(mos, .6, 1.4)
      for idir = 0, Ndirs-1 do $
         cgtext, 5+idir*xs/5, 5, align = 0, /device, color = default eq idir?'cyan':'green' $
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
  
  cams = *self.cameras

  detector = self->getdetector( cams[-1] )
  filesNB = file_search(dirs+'/'+cams[-1]+'/*.fits', count=nfilesNB)
  filesNB = red_sortfiles(filesNB)
  self->extractstates, filesNB, statesNB

  if keyword_set(unitscalib) then begin
    filesWB = file_search(dirs+'/'+cams[0]+'/*.fits', count=nfilesWB)
    detectorwb = self->getdetector( cams[0] )
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
  
  ;; load data and compute mean spectrum

  spec   = dblarr(Nstates)
  wav    = dblarr(Nstates)
  pref   = strarr(Nstates)
  specwb = dblarr(Nstates)

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

;    darkfile = file_search(self.out_dir +'/darks/'+detector+'_'+statesNB[istate].cam_settings+'.dark.fits', count=ct)
;    darkfilewb = file_search(self.out_dir +'/darks/'+detectorwb+'_'+statesWB[istate].cam_settings+'.dark.fits', count=ct)

;    darkN = red_readdata(darkfile, /silent)
;    darkW = red_readdata(darkfilewb, /silent)

    ;; Let's not assume that all images for one state must be in the
    ;; same file... just in case.
    
    pos = where(statesNB[*].fullstate eq ustate[istate], count)
;    print, inam+'loading files for state -> '+ustate[istate]

    for kk=0L, count-1 do begin
      
      imsN = float(rdx_readdata(statesNB[pos[kk]].filename))
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

    spec[istate] /= count
    if keyword_set(unitscalib) then specwb[istate] /= count
    
    wav[istate]  = statesNB[pos[0]].tun_wavelength*1.d10
    pref[istate] = statesNB[pos[0]].prefilter
    
  endfor                        ; istate


  ;; Loop prefilters

  file_mkdir, self.out_dir+'/prefilter_fits/'
  upref = pref[uniq(pref, sort(pref))] 
  npref = n_elements(upref)

  for ipref = 0L, npref-1 do begin

    red_progressbar, ipref, Npref, 'Loop over prefilters: ' + upref[ipref], /predict

    ;; copy spectra for each prefilter
    
    idx = where(pref eq upref[ipref], nwav)
    iwav = wav[idx]
    ispec = spec[idx]
    if keyword_set(unitscalib) then wbint = mean(specwb[idx]) else wbint = 1.
    
    ;; Load satlas
    red_satlas, iwav[0]-5.1, iwav[-1]+5.1, xl, yl, /si, cont = cont 
    dw = xl[1] - xl[0]
    np = round((0.080 * 8) / dw)
    if np/2*2 eq np then np -=1
    tw = (dindgen(np)-np/2)*dw + double(upref[ipref])
    tr = chromis_profile(tw, erh=-0.07d0)
    tr /= total(tr)
    yl1 = fftconvol(yl, tr)

    ;; Prepdata
    
    if Nwav gt 1 then begin

      if keyword_set(mask) then w = red_maskprefilter(iwav, ispec) else w = dblarr(n_elements(iwav)) + 1.0d0
      dat = {xl:xl, yl:yl1, ispec:ispec, iwav:iwav, pref:double(upref[ipref]), w:w}

      ;; Pars = {fts_scal, fts_shift, pref_w0, pref_dw}
      fitpars = replicate({mpside:2, limited:[0,0], limits:[0.0d, 0.0d], fixed:0, step:1.d-5}, 7)
      
      fitpars[0].limited[*] = [1,0]
      fitpars[0].limits[*] = [0.d0, 0.0d0]
      
      fitpars[1].limited[*] = [1,1]
      fitpars[1].limits[*] = [-1.0,1.0]
      
      fitpars[2].limited[*] = [1,1]
      fitpars[2].limits[*] = [-3.0d0,+3.0d0]
      
      fitpars[3].limited[*] = [1,1]
      fitpars[3].limits[*] = [2.0d0, 7.5d0]
      
      fitpars[4].limited[*] = [1,1]
      fitpars[4].limits[*] = [2.5d0, 3.5d0]
      fitpars[4].fixed = 1
      
      fitpars[5].limited[*] = [1,1]
      fitpars[5].limits[*] = [-1.d0, 1.d0]
      
      fitpars[6].limited[*] = [1,1]
      fitpars[6].limits[*] = [-1.d0, 1.d0]
      
      
      ;; Now call mpfit
      
      par = [max(ispec) * 2d0 / cont[0], 0.01d0, 0.01d0, 3.3d0, 3.0d0, -0.01d0, -0.01d0]
      par = mpfit('chromis_prefilterfit', par, xtol=1.e-4, functar=dat, parinfo=fitpars, ERRMSG=errmsg)
      prefilter = chromis_prefilter(par, dat.iwav, dat.pref)
      
      ;; save curve
      
      prf = {wav:iwav, pref:prefilter, spec:ispec, wbint:wbint, reg:upref[ipref], $
             fitpars:par, fts_model:interpol(yl1, xl+par[1], iwav)*prefilter, units:units}
      
      cgwindow
      mx = max([ispec, interpol(yl1, xl+par[1], iwav)*prefilter, prefilter/par[0] * max(ispec)]) * 1.05
      
      colors = ['blue', 'red', 'black']
      lines = [0, 2, 0]
      psyms = [16, -3, -3]
      cgplot, /add, iwav/10., ispec, line = lines[0], color = colors[0] $
              , xtitle = '$\lambda$ / 1 nm', psym = psyms[0], yrange = [0, mx]
      cgplot, /add, /over, iwav/10., interpol(yl1, xl+par[1], iwav)*prefilter $
              , color = colors[1], line = lines[1], psym = psyms[1]
      cgplot, /add, /over, iwav/10., prefilter/par[0] * max(ispec), color = colors[2], line = lines[2], psym = psyms[2]
      
      cglegend, /add, align = 3, /data $
                , location = [!x.crange[0] + (!x.crange[1]-!x.crange[0])*0.1, mean(!y.crange)*.02] $
                , title = ['obs scan'], color = colors[0], psym = psyms[0], length = 0.0
      cglegend, /add, align = 5, /data, location = [mean(!x.crange), mean(!y.crange)*.02] $
                , title = ['model scan'], line = lines[1], color = colors[1], length = 0.05
      cglegend, /add, align = 2, /data $
                , location = [!x.crange[1] - (!x.crange[1]-!x.crange[0])*0.01, mean(!y.crange)*.02] $
                , title = ['fitted prefilter'], line = lines[2], color = colors[2], length = 0.05
      
      cgcontrol, output = self.out_dir + '/prefilter_fits/chromis_'+upref[ipref]+'_prefilter.pdf'
      
    endif else begin

      y1 = interpol(yl, xl, iwav)
      prefilter = [ispec/yl1]
      prf = {wav:iwav, pref:prefilter, spec:ispec, wbint:wbint, reg:upref[ipref], $
            fitpars:prefilter, fts_model:y1, units:units}
      
    endelse
    
    save, file=self.out_dir + '/prefilter_fits/chromis_'+upref[ipref]+'_prefilter.idlsave', prf
    
  endfor                        ; ipref
  
  
end
