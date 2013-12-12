; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    timecor  : 
;   
;   
;   
;    out_dir  : 
;   
;   
;   
;    rot_dir  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
; 
; 
;-
pro red::make_stokes_crispex, timecor = timecor, out_dir = out_dir, rot_dir = rot_dir
  inam = 'red::make_stokes_crispex : '
  if(~keyword_set(rot_dir)) then rot_dir = 0
                                ;
                                ; Search reduction folders
                                ;
  fold = file_search(self.out_dir + '/momfbd/*', /test_directory, count = ct)
  if(ct eq 0) then begin
     print, inam + 'Error, no subfolders were found in '+self.out_dir+'/momfbd/'
     return
  endif
                                ;
                                ; select one of them if ct >1
                                ;
  iread = 0
  if(ct gt 1) then begin
     print, inam + 'reduction subfolders:'
     for ii = 0L, ct-1 do print, string(ii,format='(I2)') +' -> '+ file_basename(fold[ii])
     read, iread, prompt = 'Please, choose one state (enter the ID number): '
  endif
  fold = fold[iread]
  time = file_basename(fold)

  iread = 0L
  fold = file_search(fold+'/*', /test_directory, count = ct)
  if(ct gt 1) then begin
     print, inam + 'reduction subfolders:'
     for ii = 0L, ct-1 do print, string(ii,format='(I2)') +' -> '+ file_basename(fold[ii])
     read, iread, prompt = 'Please, choose one state (enter the ID number): '
  endif
  fold = fold[iread]
  
  pref = file_basename(fold)
  print, inam + 'selected state -> '+ pref
                                ;
                                ; search files
                                ;
  ss = fold + '/cfg/results/stokes/stokesIQUV.*f0'
  files = file_search(ss, count = ct)
  if(ct lt 1) then begin
     print, inam + 'Error, no files found in '+ file_basename(ss)
     return
  endif
                                ;
                                ; sort states in form (nwav, nscan) -> returns a structure
                                ; but the files are in st.ofiles.
                                ;
  st = red_get_stkstates(files)
  nwav = st.nwav
  nscan = st.nscan              ;< 140
  print, inam + 'found '+ red_stri(nscan) + ' scans'
                                ;
                                ; get residual cross-talk from Stokes I
                                ;
  k = 0L
                                ;
                                ; Load first scan 
                                ;
  imean = fltarr(nwav)
  for ii = 0L, nwav-1 do imean[ii] = median((f0(st.ofiles[ii,0]))[*,*,0])
                                ;
  ppc = red_select_spoints(st.uiwav, imean)
  
  cstk = fltarr(4, n_elements(ppc), nscan)
  ntot = 100. / (n_elements(ppc) * nscan - 1.0)
                                ;
  for ss = 0L, nscan - 1 do begin
     for ww = 0L, n_elements(ppc)-1 do begin
        if(st.flag[ww,ss]) then continue ; no image in this state
        if(k eq 0) then dim = size(f0(st.ofiles[ppc[ww],ss]), /dim)
        cstk[*,ww,ss] = red_get_cross(temporary(f0(st.ofiles[ppc[ww],ss])))
        print, string(13B), 'measuring crosstalk from I -> Q,U,V: ', k*ntot,'%', format='(A,A,F5.1,A,$)'
        k++
     endfor
  endfor
  print, ' '
                                ;
                                ; Average the wavelength-dependent crosstalk and remove outliers
                                ; using a median filter
                                ;
  !p.multi = [0,2,3]
  window, 0, xs = 600, ys = 800, title = 'Calibration details'
  cstk = total(cstk, 2) / n_elements(ppc)
  imean = reform(cstk[0,*])
  cstk = cstk[1:*,*]
                                ;
  for ii = 0, 2 do cstk[ii,*] = median(reform(cstk[ii,*]), 3)
  plot, cstk[0,*], ytitle = 'I -> Q', xtitle='tstep'
  plot, cstk[1,*], ytitle = 'I -> U', xtitle='tstep'
  plot, cstk[2,*], ytitle = 'I -> V', xtitle='tstep'
                                ; stop
                                ;
                                ; load prefilter
                                ;
  
  self->getcamtags, dir = self.data_dir
  ;; pr1 = self.out_dir + 'prefilter_fits/'+self.camttag+'.'+pref+'.prefilter.f0'
  ;; pr2 = self.out_dir + 'prefilter_fits/'+self.camrtag+'.'+pref+'.prefilter.f0'
  ;; ;
  ;; if(~file_test(pr1) AND ~file_test(pr2)) then begin
  ;;    print, inam + 'WARNING, could not find prefilter fit files, the prefilter will not be corrected!'
  ;;    prefilter = fltarr(nwav) + 1.0
  ;; endif else begin
  ;;    if(file_test(pr1) AND file_test(pr2)) then begin
  ;;       prefilter = 2.0 / ((f0(pr1) + f0(pr2)))
  ;;    endif else begin
  ;;       if(file_test(pr1)) then prefilter = 1.0 / f0(pr1) 
  ;;       if(file_test(pr2)) then prefilter = 1.0 / f0(pr2) 
  ;;    endelse
  ;; endelse
  ;;
  ;; Load prefilter
  ;;
  tpfile = self.out_dir + '/prefilter_fits/'+self.camttag+'.'+pref+'.prefilter.f0'
  tpwfile = self.out_dir + '/prefilter_fits/'+self.camttag+'.'+pref+'.prefilter_wav.f0'
  rpfile = self.out_dir + '/prefilter_fits/'+self.camrtag+'.'+pref+'.prefilter.f0'
  rpwfile = self.out_dir + '/prefilter_fits/'+self.camrtag+'.'+pref+'.prefilter_wav.f0'

  if(file_test(tpfile) AND file_test(tpwfile)) then begin
     print, inam + 'Loading:'
     print, '  -> ' + file_basename(tpfile)
     tpref = f0(tpfile)
     print, '  -> ' + file_basename(tpwfile)
     twav = f0(tpwfile)
  endif else begin
     print, inam + 'prefilter files not found!'
     stop
  endelse

  if(file_test(rpfile) AND file_test(rpwfile)) then begin
     print, inam + 'Loading:'
     print, '  -> ' + file_basename(rpfile)
     rpref = f0(tpfile)
     print, '  -> ' + file_basename(rpwfile)
     rwav = f0(tpwfile)
  endif else begin
     print, inam + 'prefilter files not found!'
     stop
  endelse
  
  prefilter = float(2.0 / (red_intepf(twav, tpref, st.uiwav*1.d-3) + red_intepf(rwav, rpref, st.uiwav*1.d-3)))

  plot, st.uiwav * 1.e-3, 1. / prefilter, ytitle = 'Prefilter', xtitle = 'Wavelength'
                                ;
                                ; load time-series-calibration data
                                ; vars: tstep, clip, tile, scale, ang, shift, grid, time, date, wfiles, tmean
                                ;
  polish = 0B
  pfile = self.out_dir + '/calib_tseries/tseries.'+pref+'.'+time+'.calib.sav'
  if(file_test(pfile)) then begin
     print, inam + 'Using time-series calibration in ' + pfile
     restore, pfile
     polish = 1B
  endif else print, inam + 'WARNING, no time-series calibration found!'
  
                                ;
                                ; Correct time intensity variations? Measure from WB
                                ;
  if(keyword_set(timecor) AND (polish)) then begin
                                ;tmean = fltarr(nscan)
                                ;for ii = 0L, nscan-1 do tmean[ii] = mean(cub[*,*,*,0,ii])
                                ;
     ntm = n_elements(tmean)
     ttime = dblarr(ntm)
     for ii = 0L, ntm-1 do ttime[ii] = red_time2double(time[ii])
     ttime -= ttime[0]
     cc = poly_fit(ttime[0:n_elements(tmean)-1], tmean, 2, tfit)
     plot, ttime / 60.d0, tmean, /line, ystyle =3
     oplot, ttime / 60.0d0, tfit
                                ;
                                ; Use flat to obtain NB variation
                                ;
     scl = self->count2diskcenter(pref) * (ttime^2. * cc[2] + ttime * cc[1] + cc[0])
                                ;
  endif else scl = fltarr(nscan) + median(imean)
                                ;
  scl = 10000. / scl
                                ;
                                ; create data cube as integer
                                ; Stokes = stokes * 1E4
                                ;  nx, ny, nwav * 4* nt
                                ;  
  k = 0L
  ntot = 100. / (ct - 1.0)
                                ;
  cub = intarr(dim[0], dim[1], nwav,4, nscan)
                                ;skip = [0,1,2,11,12,13]
  for ss = 0L, nscan - 1 do begin
     for ww = 0L, nwav - 1 do begin
        if(st.flag[ww,ss]) then continue ; no image in this state
                                ;
                                ; Load image
                                ;
        tmp = f0(st.ofiles[ww,ss])
                                ;
                                ; Correct Crosstalk
                                ;
                                ;dum =where(skip eq ww, count)
        for ii = 0, 2 do tmp[*,*,ii+1] -= cstk[ii,ss] * tmp[*,*,0]
                                ;
                                ; prefilter and scale (remember that prefilter in reality is 1/prefilter)
                                ;
        tmp *= total(prefilter[ww] * scl[ss])
                                ;
                                ; De-rotate, align and de-stretch into cube
                                ;
;     aa[ss] = round(rotate(stretch(red_rotation(temporary(cmap1), ang[ss], total(shift[0,ss]), total(shift[1,ss])), reform(grid[ss,*,*,*])), rot_dir))

        for ii = 0, 3 do begin
           cub[*,*,ww,ii,ss] = rotate(round(stretch(red_rotation(tmp[*,*,ii], ang[ss], total(shift[0,ss]), total(shift[1,ss])), reform(grid[ss,*,*,*]))), rot_dir)
        endfor
                                ;
        print, string(13B), inam + 'creating data cube -> ', k * ntot, '%', FORMAT='(A,A,F5.1,A,$)'
        k++
     endfor
  endfor
  print,' '
                                ;
                                ; Create lp header for image cube
                                ;
  extra = 'stokes=[I,Q,U,V], ns=4'
  typestring = '(integer)'
  datatype = 2
  header = 'stokes=[I,Q,U,V], ns=4 : ' +' datatype='+strtrim(datatype,2)+' '+typestring
  header = header + ', dims='+strtrim(3,2)
  header = header + ', nx='+strtrim(dim[0],2)
  header = header + ', ny='+strtrim(dim[1],2)
  header = header + ', nt='+strtrim(nscan*nwav*4L,2)
  if ((byte(1L, 0, 1))[0] eq 1) then endianstr = 'endian=l'  $ ; little endian
  else endianstr = 'endian=b'                                  ; big endian
  header = header + ', '+endianstr
                                ;
                                ; Save
                                ;
  if(~keyword_set(out_dir)) then out_dir = fold + '/cfg/results/crispex/'
  file_mkdir, out_dir
  ofile = 'crispex.'+pref+'.fullstokes.icube'
  openw, lun, out_dir+ofile, /get_lun
                                ;
                                ; Header
                                ;
  hh = bytarr(512)
  hh1 = byte(header)
  hh[0:n_elements(hh1)-1] = hh1
  writeu, lun, hh
                                ;
                                ; Data
                                ;
  print, inam + 'writting cube to file -> '+out_dir + ofile, format = '(A,$)'
  for ii = 0L, nscan - 1 do for jj =0L, 3 do writeu, lun, cub[*,*,*,jj,ii]
  free_lun, lun
  print, ' done'
                                ;
                                ; SP file
                                ;
  extra = 'stokes=[I,Q,U,V], ns=4'
  typestring = '(integer)'
  datatype = 2
  header = 'stokes=[I,Q,U,V], ns=4 : ' +' datatype='+strtrim(datatype,2)+' '+typestring
  header = header + ', dims='+strtrim(3,2)
  header = header + ', nx='+strtrim(nwav,2)
  header = header + ', ny='+strtrim(nscan,2)
  header = header + ', nt='+strtrim(dim[0]*dim[1]*4L,2)
  if ((byte(1L, 0, 1))[0] eq 1) then endianstr = 'endian=l'  $ ; little endian
  else endianstr = 'endian=b'                                  ; big endian
  header = header + ', '+endianstr
                                ;
  ofile = 'crispex.'+pref+'.fullstokes_sp.icube'
  openw, lun, out_dir+ofile, /get_lun
                                ;
                                ; Header
                                ;
  hh = bytarr(512)
  hh1 = byte(header)
  hh[0:n_elements(hh1)-1] = hh1
  writeu, lun, hh
                                ;
                                ; Data
                                ;
  print, inam + 'transposing cube dimensions ... ', format = '(A,$)'
  cub = transpose(temporary(cub), [2,4,3,0,1])
  print, inam + 'done'
  print, inam + 'writting spectral cube to file -> '+out_dir + ofile, format = '(A,$)'
  for jj = 0L, dim[1] -1 do for ii = 0L, dim[0]-1 do writeu,lun, cub[*,*,*,ii,jj]
  free_lun, lun
  print, ' done'
                                ;
  return
end
