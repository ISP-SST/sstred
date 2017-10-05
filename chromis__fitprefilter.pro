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
;    mask : in, optional, type=boolean
;
;        If selected, will allow the user to mask out spectral
;        positions from the fit.
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
;   2016-10-05 : MGL. Selection list now includes mu and number of
;                files in directories and has a default selection
;                based on those numbers.
; 
; 
;-
pro chromis::fitprefilter, time = time, scan = scan, pref = pref, mask = mask

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  if ~ptr_valid(self.data_dirs) then begin
    print, inam+' : ERROR : undefined data_dir'
    return
  endif
  dirs = *self.data_dirs

  Ndirs = n_elements(dirs)
  if( Ndirs eq 0) then begin
    print, inam+' : ERROR : no directories defined'
    return
  endif else begin
    if Ndirs gt 1 then dirstr = '['+ strjoin(dirs,';') + ']' $
    else dirstr = dirs[0]
  endelse

  ;; Find some info about the directories
  fnames = strarr(Ndirs)
  times = dblarr(Ndirs)
  Nfiles = lonarr(Ndirs)
  xs = 1920
  ys = 1200
  ims = fltarr(xs, ys, Ndirs)
  contr = fltarr(Ndirs)
  for idir = 0, Ndirs-1 do begin
    fnames[idir] = (file_search(dirs[idir]+'/Chromis-W/*fits', count = nf))[0]
    ims[0, 0, idir] = total(red_readdata(fnames[idir], /silent), 3)
    contr[idir] = stddev(ims[20:-20, 20:-20, idir])/mean(ims[20:-20, 20:-20, idir])
    Nfiles[idir] = nf
    hdr = red_readhead(fnames[idir])
    red_fitspar_getdates, hdr, date_beg = date_beg 
    times[idir] = red_time2double((strsplit(date_beg, 'T', /extract))[1])
  endfor                        ; idir

  red_logdata, self.isodate, times, mu = mu, zenithangle = za

  selectionlist = file_basename(dirs) $
                  + ' (µ=' +string(mu,format='(f4.2)') + ', #files='+strtrim(Nfiles,2) + ')'
  mos = fltarr(xs/5*Ndirs, ys/5)
  for idir = 0, Ndirs-1 do mos[idir*xs/5:(idir+1)*xs/5-1, *] $
     = rebin(ims[*, *, idir], xs/5, ys/5)/median(ims[*, *, idir])
                                   
  ;; Select data folder containing a quiet-Sun-disk-center dataset

  ;; Select default
  default = 0
  ;; Deselect directories with many files, far from disc center, with
  ;; large contrast (possibly spots?).
  findx = where(Nfiles lt 150 and mu gt 0.9 and contr lt median(contr), Nwhere)
  if Nwhere gt 0 then begin
    ;; Have to pick one!
    tmp = min(abs(times[findx]-red_time2double('13:00:00')), ml) ; Near local noon good?
    default = findx[ml]
  endif

  ;; Display some visual hints
  scrollwindow, xs = xs/5*Ndirs, ys = ys/5
  tv, bytscl(mos, .6, 1.4)
  for idir = 0, Ndirs-1 do $
     cgtext, 5+idir*xs/5, 5, align = 0, /device, color = default eq idir?'cyan':'green' $
             , strtrim(idir, 2)+' : '+red_strreplace(selectionlist[idir], 'µ', '$\mu$')

  ;; Do the selection
  tmp = red_select_subset(selectionlist $
                          , qstring = 'Select data to be processed' $
                          , count = Nselect, indx = sindx, default = default)
  idx = sindx[0]
  print, inam + ' : Will process the following data:'
  print, selectionlist[idx], format = '(a0)'
  dirs = dirs[idx]

  ;; Get files and states
  
  cams = *self.cameras
  detector = self->getdetector( cams[-1] )
  detectorwb = self->getdetector( cams[0] )

  files = file_search(dirs+'/'+cams[-1]+'/*.fits', count=nfiles)
  files1 = file_search(dirs+'/'+cams[0]+'/*.fits', count=nfiles1)
  ;;
  files = red_sortfiles(files)
  files1= red_sortfiles(files1)
  
  self->extractstates, files, states
  self->extractstates, files1, states1
  
  ;; Read one scan (scan 0 by default)
  
  if(n_elements(scan) eq 0) then scan = 0
  idx=where(states[*].scannumber eq scan, ct)
  if(ct eq 0) then begin
    print, inam+' : ERROR, invalid scan number'
    return
  endif
  
  ;; Sort selected states and files
  
  idx1 = sort(states[idx].tun_wavelength)
  ;;
  states = states[idx[idx1]]
  states1 = states1[idx[idx1]]
  
  ustate = states[uniq(states[*].tun_wavelength, sort(states[*].tun_wavelength))].fullstate  
  Nstates = n_elements(ustate)
  
  ;; load data and compute mean spectrum

  spec = dblarr(Nstates)
  wav   = dblarr(Nstates)
  pref  = strarr(Nstates)
  specwb = dblarr(Nstates)

  for istate =0L, Nstates-1 do begin

    red_progressbar, istate, Nstates, 'Loop over states: '+ustate[istate], clock = clock, /predict

    darkfile = file_search(self.out_dir +'/darks/'+detector+'_'+states[istate].cam_settings+'.dark.fits', count=ct)
    darkfilewb = file_search(self.out_dir +'/darks/'+detectorwb+'_'+states1[istate].cam_settings+'.dark.fits', count=ct)

    if(ct ne 1) then begin
      print, inam+' : ERROR, cannot find dark file'
      stop
    endif
    dd = red_readdata(darkfile, /silent)
    ddwb = red_readdata(darkfilewb, /silent)

    ;; Let's not assume that all images for one state must be in the
    ;; same file... just in case.
    
    pos = where(states[*].fullstate eq ustate[istate], count)
;    print, inam+'loading files for state -> '+ustate[istate]
    
    for kk=0L, count-1 do begin
      tmp = float(rdx_readdata(states[pos[kk]].filename))
      tmpwb = float(rdx_readdata(states1[pos[kk]].filename))
      
      dim = size(tmp,/dim)
      dx = round(dim[0]*0.12)
      dy = round(dim[1]*0.12)
      
      nsli = 1
      if(n_elements(dim) gt 2) then begin
         nsli = dim[2]
         tmp   = total(  tmp,3, /double) / double(nsli)
         tmpwb = total(tmpwb,3, /double) / double(nsli)
      endif

      tmp   -= dd
      tmpwb -= ddwb
      
      if(keyword_set(mask)) then begin

         if(kk eq 0 and istate eq 0) then begin
            mmask = red_select_area(tmp[*,*,0], /noedge, /xroi)
            nzero = where(mmask gt 0)
            bla = tmp[*,*,0]
            ind = array_indices(bla, nzero)
         endif
         
 
         tmp1 = double(tmp[reform(ind[0,*]),reform(ind[1,*])])
         tmpwb1 = double(tmpwb[reform(ind[0,*]),reform(ind[1,*])])
         
      endif else begin
 
         tmp1 = double(tmp[dx:dim[0]-dx-1,dy:dim[1]-dy-1])
         tmpwb1 = double(tmpwb[dx:dim[0]-dx-1,dy:dim[1]-dy-1])
        
      endelse ;; if mask

      spec[istate] += median(tmp1)
      specwb[istate] += median(tmpwb1)

    endfor                      ; kk

    spec[istate] /= count
    specwb[istate] /= count
    
    wav[istate] = states[pos[0]].tun_wavelength*1.d10
    pref[istate] = states[pos[0]].prefilter
  endfor                        ; istate


  ;; loop prefilters

  file_mkdir, self.out_dir+'/prefilter_fits/'
  upref = pref[uniq(pref, sort(pref))] 
  npref = n_elements(upref)

  for ipref=0, npref-1 do begin

    red_progressbar, ipref, Npref, 'Loop over prefilters: ' + upref[ipref], clock = clock, /predict

    ;; copy spectra for each prefilter
    
    idx = where(pref eq upref[ipref], nwav)
    iwav = wav[idx]
    ispec = spec[idx]
    wbint = mean(specwb[idx])
    
    
    ;; Load satlas
    red_satlas, iwav[0]-5.1, iwav[-1]+5.1, xl, yl, /si, cont = cont 
    dw = xl[1] - xl[0]
    np = round((0.080 * 8) / dw)
    if(np/2*2 eq np) then np -=1
    tw = (dindgen(np)-np/2)*dw + double(upref[ipref])
    tr = chromis_profile(tw, erh=-0.07d0)
    tr/=total(tr)
    yl1 = fftconvol(yl, tr)

    ;;units = 'Watt/(s m2 Hz ster)' ; SI units
    units = 'W s^-1 m^-2 Hz^-1 sr^-1' ; SI units
    
    ;; Prepdata
    
    if(nwav gt 1) then begin
      if(keyword_set(mask)) then w = red_maskprefilter(iwav, ispec) else w = dblarr(n_elements(iwav)) + 1.0d0
      dat = {xl:xl, yl:yl1, ispec:ispec, iwav:iwav, pref:double(upref[ipref]), w:w}

      ;; Pars = {fts_scal, fts_shift, pref_w0, pref_dw}
      fitpars = replicate({mpside:2, limited:[0,0], limits:[0.0d, 0.0d], fixed:0, step:1.d-5}, 7)
      ;;
      fitpars[0].limited[*] = [1,0]
      fitpars[0].limits[*] = [0.d0, 0.0d0]
      ;;
      fitpars[1].limited[*] = [1,1]
      fitpars[1].limits[*] = [-1.0,1.0]
      ;;
      fitpars[2].limited[*] = [1,1]
      fitpars[2].limits[*] = [-3.0d0,+3.0d0]
      ;;
      fitpars[3].limited[*] = [1,1]
      fitpars[3].limits[*] = [2.0d0, 7.5d0]
      ;;
      fitpars[4].limited[*] = [1,1]
      fitpars[4].limits[*] = [2.5d0, 3.5d0]
      fitpars[4].fixed = 1
      ;;
      fitpars[5].limited[*] = [1,1]
      fitpars[5].limits[*] = [-1.d0, 1.d0]
      ;;
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
      cgplot, /add, iwav, ispec, line = 1
      cgplot, /add, /over, iwav, interpol(yl1, xl+par[1], iwav)*prefilter
      cgplot, /add, /over, iwav, prefilter/par[0] * max(ispec), line=2
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
