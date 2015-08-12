pro red::init_data_destretch, pref = pref, scan = scan, min = min, max = max, smooth = smooth, bad = bad, $
                                 clip = clip, tile = tile, verbose=verbose, addwb = addwb, show = show, $
                                 t0 = t0, t1 = t1, scale = scale, tstep = tstep, norotation = norotation, $
                                 xbd = xbd, ybd = ybd, np = np
 ;; Get procedure name
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0]) + ': '


  ;; Defaults
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(smooth) eq 0) then smooth = 3.0
  if(n_elements(clip) eq 0) then clip = [72,  24, 8, 4, 2]
  if(n_elements(tile) eq 0) then tile = [6, 16,   32, 48, 64]
  if(n_elements(nthreads) eq 0) then nthreads = 6L
  if(keyword_set(verbose)) then verb=1 else verb=0
  if(~keyword_set(scale)) then scale = 1.0 / float(self.image_scale)
  if(~keyword_set(norotation)) then rot = 1 else rot = 0
  if(n_elements(xbd) eq 0) then xbd = 512
  if(n_elements(ybd) eq 0) then ybd = 512
  if(n_elements(np) eq 0)  then np = 3
  cubic = 0

  
  ;; Select folder
  dir = file_search(self.out_dir+'/data/*',/test_dir, count = ct)
  if(ct eq 0) then begin
     print, inam+'ERROR, no data found in -> '+ self.out_dir+'/data/'
     return
  endif
  if(ct gt 1) then begin
     print, inam+'Found folders:'
     for ii=0,ct-1 do print, string(ii,format='(I4)')+' -> '+file_basename(dir[ii])
     idx = 0L
     read, idx, prompt='Select folder ID: '
     dir = dir[idx]+'/'
  endif
  print, inam+'Selected folder: -> '+dir
  
  

  ;; Get files
  files = file_search(dir +'/'+ self.camwb+'/cam*', count = nf)
  print, inam + 'Found '+string(nf,format='(I0)')+' files'
  


  ;; Get states
  st = red_getstates(files, /links)
  if(n_elements(pref) gt 0) then begin
     pos = where(st.pref eq pref, count)
     if(count eq 0) then begin
        print, inam+'ERROR, pref '+pref+' does not exist'
        return
     endif
     files = files[pos]
  endif else begin
     upref = st.pref[uniq(st.pref,sort(st.pref))]
     if(n_elements(upref) gt 1) then begin
        print, inam+'Found prefilter:'
        for ii=0,n_elements(upref)-1 do print, string(ii,format='I4')+' -> '+upref[ii]
        pos  =0
        read,pos, prompt='Select prefilter ID: '
        pos  = where(st.pref eq upref[pos])
        pref = upref[pos]
        files = files[pos]
     endif else pref=upref
  endelse
  
  files = red_sortfiles(temporary(files))
  st = red_getstates(files, /links)
  

  ;; Get unique states
  uscan = st.scan[uniq(st.scan,sort(st.scan))]
  pos = uniq(st.wav,sort(st.dwav))
  uwav = st.wav[pos]
  udwav = st.dwav[pos]
  idx = sort(udwav)
  udwav = udwav[idx]
  uwav = uwav[idx]
  nscan = n_elements(uscan)
  nwav = n_elements(uwav)
  ulc = st.lc[uniq(st.lc,sort(st.lc))]
  nlc = n_elements(ulc)
  fstate = st.scan+'.'+st.state
  ufstate = fstate[uniq(fstate)]
  
  upref = st.pref[uniq(st.pref,sort(st.pref))]
  
  self->getcamtags, dir = dir
  cams = [self.camwbtag]
  cc = 0

  ;; Output data
  file_mkdir, self.out_dir + '/calib_tseries/'
  ofile = self.out_dir +'/calib_tseries/header.add.'+pref+'.'+file_basename(dir)+'.idlsave'

  ;; Load alignclips
  file = self.out_dir+'/calib/align_clips.'+pref+'.sav'
  if(~file_test(file)) then begin
     print, inam+'ERROR, you must perform the pinhole calibration'
     return
  endif else restore, file

  
  ;;get number of adquisitions
  numadd = intarr(n_elements(uwav), n_elements(ulc), nscan)
  times = strarr(nscan)
  tstate = st.scan+'.'+st.pref

  print, inam+'Reading headers/time stamps ... ', format='(A,$)'
  for ss=0, nscan-1 do begin
     for ll = 0, n_elements(ulc)-1 do for ww = 0,n_elements(uwav)-1 do begin
        istate = strjoin([uscan[ss], pref, uwav[ww], ulc[ll]], '.')
        
        pos = where(fstate eq istate, count)
        numadd[ww,ll,ss] = count
     endfor

     ;;
     ;; Get time stamps for each scan for derotation and polarimetry
     ;; 
     istate = strjoin([uscan[ss], pref], '.')
     pos = where(tstate eq istate, count)
     
     if(count eq 0) then begin
        print, inam+'Error, state not found -> ', istate
        return
     endif
     
    ;; fzread, dum, files[pos[count/2]], h
     times[ss] = fzhead(files[pos[count/2]])
     len = strpos(times[ss],'Te=')
     times[ss] = strmid(times[ss],len+14,12)
  endfor
  print, 'done'
  
  
  if(self.isodate eq '') then begin
     date = ''
     read, date, prompt=inam+'please enter the ISO-date of the observations (YYYY-MM-DD): '
  endif else date = self.isodate
  date = strjoin(strsplit(date,'-./', /extract),'-')
  
  ang = red_lp_angles(times, date)
  mang = mean(ang)
  ang -= mang

  maxadd = max(numadd)

  print, inam+'Max number of acquisitions -> '+string(maxadd,format='(I0)')
  d = fltarr(sx, sy, maxadd, n_elements(uwav), n_elements(ulc))
  if(n_elements(refs) eq 0) then refs = fltarr(sx,sy,nscan)
  if(n_elements(ifiles) eq 0) then ifiles = strarr(maxadd, nwav, nlc, nscan)
  


  ;;
  ;; load darks and flats
  ;;
  nam = self.out_dir + '/darks/'+cams+'.dark'
  if(~file_test(nam)) then begin
     print,inam+'ERROR, darks not found'
     return
  endif
  
  ;; Create variables and load flat and dark
  print, inam+'loading dark '+file_basename(nam)
  dd = red_clipim(f0(nam), cl[*,0])
  dim = size(dd,/dim)
  ffile = self.out_dir+'/flats/'+cams+'.'+pref+'.flat' 
  if(~file_test(ffile)) then begin
     print, unam+' ERROR, flat not found -> '+ffile
  endif
  print, inam+'loading flat '+ffile
  ff = red_flat2gain(red_clipim(f0(ffile), cl[*,0]), mi = min, ma = max, $
                     smooth = smooth, bad = bad, /preserve)     
  
  ;;
  ;; get reference images based on file size
  ;;

  for ss = 0, nscan-1 do begin
     ;; get file size
     pos = where(st.scan eq uscan[ss], count)
     size = double((file_info(files[pos])).size)
     
     ;; remove trends due to intensity increase
     dum = poly_fit(dindgen(count), size, 1, yfit = norm)
     norm = mean(norm) / norm
     size *= norm
     
     ;; Normalize to units ~ 1
     mi = min(size)
     ma = max(size)
     range = 1.d0 / (ma - mi)
     size = (size - mi) * range
     
     ;; take second best? Hopefully not outlier
     idx = sort(size)
     refs[*,*,ss] = red_fillpix((red_clipim(f0(files[pos[idx[(count-2)>0]]]), cl[*,0]) -dd) * ff, nthreads=nthreads)
     if(rot) then refs[*,*,ss] = red_rotation( refs[*,*,ss], ang[ss])
  endfor
  
  
  shifts_refs = red_aligncube(refs, np, xbd = xbd, ybd = ybd, cubic = cubic, /aligncube)
  


  for ss = 0, nscan - 1 do begin
     if(n_elements(scan) ne 0) then if(uscan[ss] ne string(scan, format='(I05)')) then continue
     print, inam+' scan='+string(ss,format='(I0)'), format='(A,$)'
     for ll = 0, n_elements(ulc)-1 do for ww = 0, n_elements(uwav)-1 do begin
;        pos = where(st.wav eq uwav[ww] AND st.lc eq ulc[ll] AND
;        st.scan eq uscan[ss], count)
        istate = strjoin([uscan[ss], pref, uwav[ww], ulc[ll]], '.')
        pos = where(fstate eq istate, count)

        ifiles[0:count-1,ww,ll,ss] = files[pos]
           
     endfor
     print, string(13B), format='(A,$)'
  endfor
  print,' '



  ;;
  ;; Polish time series
  ;;
  dt = red_time2double(times)
  dt = mean(dt[1:*] - dt[0:*])
  tstep = float(round(180. / dt)) +1 ;; 3 minutes!
  scale = 1.0 / self.image_scale

 ; refs1 = refs
  
  ;shifts_refs = red_aligncube(refs1, np, xbd = xbd, ybd = ybd, cubic = cubic, /aligncube)
  print, inam + "computing time-series polishing ... ", format='(A,$)'
  tgrid = red_destretch_tseries(refs, scale, tile, clip, tstep)
  print, 'done'
  

  print, inam+'saving '+file_basename(ofile)
  
  save, file = ofile, ifiles, files, refs, times, mang, ang, times, $
        st, shifts_refs, tile, clip, uwav, udwav, ulc, uscan, fstate,$
        nlc, nwav, nscan, numadd, ff, dd, pref, upref,  $
        size, date, cl, sx, sy, maxadd, tgrid, tstep, dt, scale
  

  
end
