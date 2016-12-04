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
;     cgs : in, optional, type=boolean
;   
;        Set this to use CGS units (default is to normalize with the
;        continuum). 
; 
;     si : in, optional, type=boolean
;   
;        Set this to use SI units (default is to normalize with the
;        continuum).
;
;    mask : in, optional, type=boolean
;
;           if selected, will allow the user to mask out spectral
;           positions from the fit.
; 
; 
; :History:
; 
;   2016-11-28 : MGL. Moved helper functions to their own files and
;                added header. Make and save a final plot of the fit.
;                Prevent user from setting both /cgs and /si keywords.
; 
;   2016-12-04 : JdlCR. Allow to mask regions of the mean
;                spectra. Sometimes there is no real quiet-sun and the line center
;                must be masked.
;
; 
; 
;-
pro chromis::fitprefilter, time = time, scan = scan, pref = pref, cgs = cgs, si = si, mask = mask

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])+': '

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  if keyword_set(si) and keyword_set(cgs) then begin
    print, inam+' : Please do not set both /si and /cgs.'
    retall
  endif

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

  ;; Select data folder containing a quiet-Sun-disk-center dataset

  idx = 0L
  if(Ndirs gt 1) then begin
    print, inam+'Time stamps found:'
    for ii = 0, Ndirs-1 do print, '   '+red_stri(ii)+' -> '+file_basename(dirs[ii])
    read, idx, prompt='Please select folder ID: '
  endif
  dirs = dirs[idx]
  print, inam+'Selected -> '+dirs
  
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
    print, inam+'ERROR, invalid scan number'
    return
  endif
  
  ;; Sort selected states and files
  
  idx1 = sort(states[idx].tun_wavelength)
  ;;
  states = states[idx[idx1]]
  states1 = states1[idx[idx1]]
  
  ustate = states[uniq(states[*].tun_wavelength, sort(states[*].tun_wavelength))].fullstate  
  nstate = n_elements(ustate)
  
  ;; load data and compute mean spectrum

  spec = dblarr(nstate)
  wav   = dblarr(nstate)
  pref  = strarr(nstate)
  specwb = dblarr(nstate)

  for ii =0L, nstate-1 do begin
    
    darkfile = file_search(self.out_dir +'/darks/'+detector+'_'+states[ii].cam_settings+'.dark.fits', count=ct)
    darkfilewb = file_search(self.out_dir +'/darks/'+detectorwb+'_'+states1[ii].cam_settings+'.dark.fits', count=ct)

    if(ct ne 1) then begin
      print, inam+'ERROR, cannot find dark file'
      stop
    endif
    dd = red_readdata(darkfile)
    ddwb = red_readdata(darkfilewb)

    ;; Let's not assume that all images for one state must be
    ;; in the same file... just in case.
    
    pos = where(states[*].fullstate eq ustate[ii], count)
    print, inam+'loading files for state -> '+ustate[ii]
    
    for kk=0L, count -1 do begin
      tmp = red_readdata(states[pos[kk]].filename)
      tmpwb = red_readdata(states1[pos[kk]].filename)
      
      dim = size(tmp,/dim)
      dx = round(dim[0]*0.12)
      dy = round(dim[1]*0.12)
      
      nsli = 1
      if(n_elements(dim) gt 2) then nsli = dim[2]
      
      if(nsli gt 1) then begin
        tmp1 = total(tmp[dx:dim[0]-dx-1,dy:dim[1]-dy-1,*],3, /double) / double(nsli)
        tmpwb1 = total(tmpwb[dx:dim[0]-dx-1,dy:dim[1]-dy-1,*],3, /double) / double(nsli)
      endif else begin
        tmp1 = double(tmp[dx:dim[0]-dx-1,dy:dim[1]-dy-1])
        tmpwb1 = double(tmpwb[dx:dim[0]-dx-1,dy:dim[1]-dy-1])
      endelse

      spec[ii] += median(tmp1-dd)
      specwb[ii] += median(tmpwb1-ddwb)

    endfor                      ; kk

    spec[ii] /= count
    specwb[ii] /= count
    
    wav[ii] = states[pos[0]].tun_wavelength*1.d10
    pref[ii] = states[pos[0]].prefilter
  endfor                        ; ii


  ;; loop prefilters

  file_mkdir, self.out_dir+'/prefilter_fits/'
  upref = pref[uniq(pref, sort(pref))] 
  npref = n_elements(upref)

  for pp=0, npref-1 do begin

    ;; copy spectra for each prefilter
    
    idx = where(pref eq upref[pp], nwav)
    iwav = wav[idx]
    ispec = spec[idx]
    wbint = mean(specwb[idx])
    
    
    ;; Load satlas
    
    red_satlas, iwav[0]-0.1, iwav[-1]+0.1, xl, yl, cgs=cgs, si=si, cont = cont
    dw = xl[1] - xl[0]
    np = round((0.080 * 8) / dw)
    if(np/2*2 eq np) then np -=1
    tw = (dindgen(np)-np/2)*dw + double(upref[pp])
    tr = chromis_profile(tw, erh=-0.02d0)
    tr/=total(tr)
    yl1 = fftconvol(yl, tr)
    
    
    ;; Prepdata
    
    if(nwav gt 1) then begin
       w = red_maskprefilter(iwav, ispec)
       dat = {xl:xl, yl:yl1, ispec:ispec, iwav:iwav, pref:double(upref[pp]), w:w}

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
      fitpars[3].limits[*] = [2.0d0, 5.0d0]
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

      par = [max(ispec) * 10d0 / cont[0], 0.01d0, 0.01d0, 3.3d0, 3.0d0, -0.01d0, -0.01d0]
      par = mpfit('chromis_prefilterfit', par, xtol=1.e-4, functar=dat, parinfo=fitpars, ERRMSG=errmsg)
      prefilter = chromis_prefilter(par, dat.iwav, dat.pref)
      
      ;; save curve

      prf = {wav:iwav, pref:prefilter, spec:ispec, wbint:wbint, reg:upref[pp], $
             fitpars:par, fts_model:interpol(yl1, xl+par[1], iwav)*prefilter}

      cgwindow
      cgplot, /add, iwav, ispec, line = 1
      cgplot, /add, /over, iwav, interpol(yl1, xl+par[1], iwav)*prefilter
      cgplot, /add, /over, iwav, chromis_prefilter(par, iwav, pref)/par[0] * max(ispec), line=2
      cgcontrol, output = self.out_dir + '/prefilter_fits/chromis_'+upref[pp]+'_prefilter.pdf'

    endif else begin

      y1 = interpol(yl, xl, iwav)
      prefilter = [ispec/yl1]
      prf = {wav:iwav, pref:prefilter, spec:ispec, wbint:wbint, reg:upref[pp], $
             fitpars:prefilter, fts_model:y1}

    endelse

    save, file=self.out_dir + '/prefilter_fits/chromis_'+upref[pp]+'_prefilter.idlsave', prf

  endfor                        ; pp
  
  
end
