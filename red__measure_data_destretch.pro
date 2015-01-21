pro red::measure_data_destretch, clip = clips, tile = tiles, verbose=verbose, overwrite = overwrite, $
                                 t0 = t0, t1 = t1, norotation = norotation, notimecor = notime
                                 
 ;; Get procedure name
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0]) + ': '


  ;; Defaults
  if(n_elements(clip) eq 0) then clips = [32, 16,  8,  4,  2]
  if(n_elements(tile) eq 0) then tiles = [5,  12, 32, 48, 64]
  if(n_elements(nthreads) eq 0) then nthreads = 6L
  if(~keyword_set(norotation)) then rot = 1 else rot = 0
  
  cubic = 0

  
  ;; Select folder
  dir = file_search(self.out_dir+'/calib_tseries/header.add.*', count = ct)
  if(ct eq 0) then begin
     print, inam+'ERROR, no data found in -> '+ self.out_dir+'/calib_tseries/'
     return
  endif
  
  if(ct gt 1) then begin
     print, inam+'Found folders:'
     for ii=0,ct-1 do print, string(ii,format='(I4)')+' -> '+file_basename(dir[ii])
     idx = 0L
     read, idx, prompt='Select folder ID: '
     dir = dir[idx]+'/'
  endif
  print, inam+'loading '+dir
  restore, dir
  tobs = (strsplit(dir, '.',/extract))[-2]
  odir1 =  self.out_dir+'/calib_tseries/wb/'+tobs+'/'
  file_mkdir, odir1
  
  if(n_elements(t0) eq 0) then t0 = 0
  if(n_elements(t1) eq 0) then t1 = nscan-1

  sx = max(cl[0:1,0])  -  min(cl[0:1,0]) + 1
  sy = max(cl[2:3,0])  -  min(cl[2:3,0]) + 1

  ;; Correct time series?
  if(~keyword_Set(notime)) then begin
     for ss=t0, t1 do refs[*,*,ss] = red_stretch(refs[*,*,ss], reform(tgrid[ss,*,*,*]))
     tcor = 1
  endif else tcor = 0
  
;; Loop scans
  for ss = t0, t1 do begin

          
     ;; ofile?
     ofile = self.out_dir + '/calib_tseries/tstep.add.'+pref+'.'+tobs+'.'+string(ss,format='(I05)')+'.idlsave'
     if(file_test(ofile) AND ~keyword_set(overwrite)) then begin
        print, inam+'skipping tt='+string(ss, format='(I0)')+', file exists!'
        continue
     endif
     
     ;; define arrays to store data and shifts
     d = fltarr(sx, sy, maxadd, n_elements(uwav), n_elements(ulc))
     shifts = fltarr([2,maxadd,n_elements(uwav), n_elements(ulc)]) 

     nadd = 0L
     print, inam+'loading/processing for scan '+string(ss, format='(I0)')
     for ll = 0, n_elements(ulc)-1 do for ww = 0, n_elements(uwav)-1 do begin
        
        for nn = 0L, numadd[ww,ll,ss]-1 do begin
           d[*,*,nn,ww,ll] = red_fillpix((red_clipim(f0(ifiles[nn,ww,ll,ss]), cl[*,0]) -dd) * ff, nthreads=nthreads) 
           if(rot) then d[*,*,nn,ww,ll] = red_rotation( d[*,*,nn,ww,ll], ang[ss])
           
           shifts[*,nn,ww,ll] = round(red_shc(refs[256:767,256:767,ss], d[256:767,256:767,nn,ww,ll], /filt))
           d[*,*,nn,ww,ll] = red_shift_im(d[*,*,nn,ww,ll], (shifts[0,nn,ww,ll]), (shifts[1,nn,ww,ll]))
           
           tmp = red_dsgridnest(refs[*,*,ss], d[*,*,nn,ww,ll], tiles, clips)
           pos = where(~finite(tmp), ctf)
           if(ctf gt 0) then tmp[pos] = 0.0
           
           if(nn eq 0 AND ll eq 0 AND ww eq 0) then begin
              dim = size(tmp, /dim)
              corrs = fltarr([dim, maxadd, n_elements(uwav), n_elements(ulc)])
              avcorr = fltarr(dim)
           endif
           
           corrs[*,*,*,nn,ww,ll] = temporary(tmp)
           avcorr += corrs[*,*,*,nn,ww,ll]
           nadd++      
        endfor
     endfor
     
     avcorr /= float(nadd)
     print, 'avx = ', mean(avcorr[0,*,*])
     print, 'avy = ', mean(avcorr[1,*,*])
     
     ;; Remove average distortion from measurements
     wb = fltarr(sx, sy)
     for ll = 0, n_elements(ulc)-1 do begin
        for ww = 0,n_elements(uwav)-1 do for nn = 0, numadd[ww,ll,ss]-1 do begin
           corrs[*,*,*,nn,ww,ll] -= avcorr
           wb += red_stretch(d[*,*,nn,ww,ll], corrs[*,*,*,nn,ww,ll])
        endfor
     endfor
     wb /= float(nadd)
     writefits, odir1+'wb.'+pref+'.'+uscan[ss]+'.fits', wb, times[ss]
     print, ' done, saving '+file_basename(ofile)

     ;; save
     save, file=ofile, corrs, shifts, rot, tiles, clips, tcor
 
  endfor   ;; ss (Scans)

  
  
  return
end
