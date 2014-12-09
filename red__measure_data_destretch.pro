pro red::measure_data_destretch, pref = pref, scan = scan, min = min, max = max, smooth = smooth, bad = bad, $
                                 clip = clip, tile = tile, verbose=verbose, addwb = addwb, show = show, $
                                 t0 = t0, t1 = t1, scale = scale, tstep = tstep, norotation = norotation
 ;; Get procedure name
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0]) + ': '


  ;; Defaults
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(smooth) eq 0) then smooth = 3.0
  if(n_elements(clip) eq 0) then clip = [32,  16, 8, 4, 2]
  if(n_elements(tile) eq 0) then tile = [6, 16, 32, 48, 64]
  if(n_elements(nthreads) eq 0) then nthreads = 6L
  if(keyword_set(verbose)) then verb=1 else verb=0
  if(~keyword_set(scale)) then scale = 1.0 / float(self.image_scale)
  if(~keyword_set(norotation)) then rot = 1 else rot = 0

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
  
  self.getcamtags, dir = dir
  cams = [self.camwbtag]
  cc = 0

  ;; Output data
  file_mkdir, self.out_dir + '/calib_tseries/'
  if(n_elements(done) eq 0) then done = bytarr(n_elements(uscan))
  ofile = self.out_dir +'/calib_tseries/tseries.add.'+pref+'.'+file_basename(dir)+'.idlsave'

  if(file_test(ofile)) then begin
     restore, ofile
     print, inam+'Resuming from file -> ', ofile
  endif else print, inam+'Calibration data to ', ofile
  
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
  if(n_elements(rms) eq 0) then rms = fltarr(maxadd, n_elements(uwav), n_elements(ulc))
  if(n_elements(refs) eq 0) then refs = fltarr(sx,sy,nscan)
  if(n_elements(ifiles) eq 0) then ifiles = strarr(maxadd, nwav, nlc, nscan)
  if(keyword_set(show)) then window, 0, xs = sx, ys=sy, title='Reference image for current scan'
  

  if(n_elements(t0) eq 0) then t0 = 0
  if(n_elements(t1) eq 0) then t1 = nscan-1
  if(n_elements(trefs) eq 0) then trefs = dblarr(nscan)


;; Loop scans
  for ss = t0,t1 do begin
     ;; Get darks and flats in the first time-step
     if(ss eq 0) then begin
        nam = self.out_dir + '/darks/'+cams+'.dark'
        if(~file_test(nam)) then begin
           print,inam+'ERROR, darks not found'
           return
        endif
        
        ;; Create variables and load flats and darks
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
     endif ;; (cc eq 0)

     if(n_elements(scan) ne 0) then if(uscan[ss] ne string(scan, format='(I05)')) then continue
     if(done[ss]) then begin
        print, inam+'Scan '+string(ss,format='(I0)')+' is already processed -> skipping'
        continue
     endif

     
     ;; Load images, get RMS and select reference
     rms[*] = 0.0
     mrms = 0.0
     print, inam+'loading WB data for scan '+string(ss, format='(I0)'), format ='(A,$)'
     for ll = 0, n_elements(ulc)-1 do for ww = 0, n_elements(uwav)-1 do begin
        pos = where(st.wav eq uwav[ww] AND st.lc eq ulc[ll] AND st.scan eq uscan[ss], count)
        for nn = 0L, count-1 do begin
           if(verb) then print, inam+'loading file -> '+file_basename(files[pos[nn]])
           ifiles[nn,ww,ll,ss] = files[pos[nn]]
           d[*,*,nn,ww,ll] = red_fillpix((red_clipim(f0(files[pos[nn]]), cl[*,0]) -dd) * ff, nthreads=nthreads)
           
           if(rot) then d[*,*,nn,ww,ll] = red_rotation( d[*,*,nn,ww,ll], ang[ss])
           
           rms[nn,ww,ll] = stdev(d[*,*,nn,ww,ll]) / mean(d[*,*,nn,ww,ll])
           if(rms[nn,ww,ll] gt mrms) then begin
              mrms = rms[nn,ww,ll]
              refs[*,*,ss] = d[*,*,nn,ww,ll]
           endif
        endfor
     endfor
     print, ' ... done'
     
     if(keyword_set(addwb)) then refs[*,*,ss] = total(total(total(d,3),3),3)/total(numadd[*,*,ss])
     if(keyword_set(show)) then tvscl, histo_opt(refs[*,*,ss])
     

     ;; Get in-scan destretch grids
     print, inam+'Computing de-stretch grids and shifts (this may take a while!) ... ', format='(A,$)'
     for ll = 0, n_elements(ulc)-1 do begin
        for ww = 0,n_elements(uwav)-1 do for nn=0, numadd[ww,ll,ss]-1 do begin
           if(n_elements(shifts) eq 0) then shifts = intarr([2,maxadd,n_elements(uwav), n_elements(ulc),nscan]) 
           
           ;; Align images to reference by shifting array (pixel accuracy
           ;; to avoid extra interpolations)
           shifts[*,nn,ww,ll,ss] = round(red_shc(refs[*,*,ss], d[*,*,nn,ww,ll], /filt))
           d[*,*,nn,ww,ll] = shift(d[*,*,nn,ww,ll], shifts[*,nn,ww,ll])
           
           tmp = red_dsgridnest(refs[*,*,ss], d[*,*,nn,ww,ll], tile, clip)
           if(n_elements(corrs) eq 0) then begin
              dim = size(tmp, /dim)
              corrs = fltarr([dim,maxadd,n_elements(uwav), n_elements(ulc),nscan]) 
           endif
           corrs[*,*,*,nn,ww,ll,ss] = temporary(tmp)
           d[*,*,nn,ww,ll] = stretch(d[*,*,nn,ww,ll], corrs[*,*,*,nn,ww,ll,ss])
        endfor
     endfor
     print, 'done'
     
     done[ss] = 1
     save, file=ofile, corrs, shifts, ifiles, files, st, rms, $
           refs, uscan, uwav, ulc, fstate, ufstate, nwav, nlc, $
           nscan, uscan, pref, upref, done, tile, clip, udwav, $
           numadd, cl, times, rot, ang, mang
  endfor   ;; cc (Scans)
  
  ;tgrid = red_destretch_tseries(refs, scale, tile, clip, tstep)

  return
end
