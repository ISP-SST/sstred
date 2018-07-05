
pro mencode, root, ofile, fps=fps, bitrate = bitrate, verbose = verbose, nomp4box=nomp4box
  inam = 'mencode : '
  ;
  if(~keyword_set(fps)) then fps = 25
  if(~keyword_set(bitrate)) then bitrate = 30000
  if(~keyword_set(nomp4box)) then vof = 'rawvideo' else vof='avi'
  ;
  bit = strcompress(string(long(bitrate)), /remove)
  fp = strcompress(string(fps), /remove)
  ;nt = strcompress(string(nthreads), /remove)
  ;
  print, inam+'Input  -> '+root
  print, inam+'Output -> '+ofile+'.mp4'
  print, inam+'bitrate -> '+bit+' Kbits/s'
  print, inam+'fps -> '+fp
  ;
  command = 'mencoder -noskip '+root+' -o /dev/null -ovc x264 -x264encopts subq=7:pass=1:bitrate='+bit+':frameref=5:bframes=1:me=umh:partitions=all:trellis=1:qp_step=4:qcomp=0.9:direct_pred=auto:keyint=30 -ofps '+fp
  ;
  print, inam+'processing pass 1'
  print,' '
  if(keyword_set(verbose)) then spawn, command else spawn, command, h
  ;
  command = 'mencoder -noskip '+root+' -o '+ofile+'.h264 -of '+vof+' -ovc x264 -x264encopts subq=7:pass=2:bitrate='+bit+':frameref=5:bframes=1:me=umh:partitions=all:trellis=1:qp_step=4:qcomp=0.9:direct_pred=auto:keyint=30 -ofps '+fp
  ;
  print, inam+'processing pass 2'
  print,' '
  if(keyword_set(verbose)) then spawn, command else spawn, command, h
  ;
  if(~keyword_set(nomp4box)) then begin
     ;
     spawn, 'rm '+ofile+'.mp4'
  ;
     print, inam+'Creating mp4 file'
     spawn, 'MP4Box -fps '+fp+' -add ./'+ofile+'.h264 '+ofile+'.mp4'
  ;
     spawn, 'rm '+ofile+'.h264'
     print, inam+'Output -> '+ofile+'.mp4'
  endif else begin
     spawn, 'mv '+ofile+'.h264 '+ofile+'.avi'
     print, inam+'Output -> '+ofile+'.avi'

  endelse
  ;
end

function red_quicklook_chromis_states, f

  ;; Create struct to store data
  
  nf = n_elements(f)
  st = {files:f, scan:strarr(nf), stat:strarr(nf), pref:strarr(nf), hr:strarr(nf), cam:strarr(nf)}

  for ii = 0L, nf-1 do begin
     tmp = strsplit(file_basename(f[ii],'.fits'), '_',/extract)
     st.scan[ii] = tmp[2]
     st.pref[ii] = tmp[4]
     st.hr[ii] = tmp[5]
     st.stat[ii] = tmp[4] +'_'+tmp[5]
     st.cam[ii] = tmp[1]
  endfor
  
  
  return, temporary(st)
end

function red_quicklook_chromis_extract, st, stat
  
  idx = where(st.stat eq stat, ct)
  res = {files:st.files[idx], scan:st.scan[idx], stat:stat[idx], pref:st.pref[idx], hr:st.hr[idx], cam:st.cam[idx]}

  return, res
end

function rms_select, tmp

  dim = size(tmp, /dim)
  if(n_elements(dim) eq 2) then return, tmp

  x0 = dim[0]/3
  x1 = dim[0]-dim[0]/3
  y0 = dim[1]/3
  y1 = dim[1]-dim[0]/3
  
  rms = fltarr(dim[2])
  for ii = 0, dim[2] -1 do begin
     rms[ii] = stdev(tmp[x0:x1,y0:y1,ii]) / median(tmp[x0:x1,y0:y1,ii])
  endfor

  dum = max(rms,p)

  
  return, tmp[*,*,p]
end



pro red_quicklook_chromis, root, dark=dark, flat=flat, cam=cam, $
                           fps = fps, bitrate = bitrate, clip=clip, $
                           statid = statid, foldid = foldid, getfolders = getfolders, $
                           getstates = getstates, states = states, folders = folders, $
                           rms_select = rms_select, rebin = rebin, idlcompress=idlcompress, no_normalize=no_normalize
  


  ;; Defaults
  
  if(n_elements(dark) eq 0) then dark = 0.0
  if(n_elements(flat) eq 0) then flat = 1.0
  if(n_elements(cam) eq 0) then cam = 'Chromis-N/'
  if(n_elements(fps) eq 0) then fps = 8
  if(n_elements(bitrate) eq 0) then bitrate = 40000
  if(n_elements(clip) ne 4) then clip = [0,0,0,0]
  if(n_elements(flat) gt 1) then dum = red_flat2gain(flat,gain_nozero=gain) else gain = 1.0/flat
  
 

  
  ;; Search folders
  
  idx = 0
  fold = file_search(root+'/*', /test_dir, count = nf)
  
  if(nf eq 0) then begin
     print, 'ERROR, could not find folders in '+root
     return
  endif

  if(keyword_set(getfolders)) then begin
     folders = fold
     return
  endif
  
  
  
  times = strarr(nf)
   if(nf gt 1 and (n_elements(foldid) eq 0)) then print, 'Found folders:'
  
  for ii = 0, nf - 1 do begin
     times[ii] = (strsplit(fold[ii], '/',/extract))[-1]
      if(nf gt 1 and (n_elements(foldid) eq 0)) then print, string(ii, format='(I5)')+' -> '+ times[ii]
  endfor
  
  if(nf gt 1 and (n_elements(foldid) eq 0)) then read, idx, prompt='Please type folder ID: '

  if(n_elements(foldid) eq 1) then begin
     if(foldid ge 0 AND foldid lt nf) then idx = foldid else begin
        print, 'keyword foldid out of range, exiting! -> foldid='+string(foldid)
        return
     endelse
  endif

  
  fold = fold[idx]+'/'+cam
  times = times[idx]
  print, 'Using -> '+fold

  f = file_search(fold+'/*', count = nf)



  


  ;; Parse states based on file names
  
  st = red_quicklook_chromis_states(f)
  ustat = st.stat[uniq(st.stat, sort(st.stat))]
  ns = n_elements(ustat)


  if(keyword_set(getstates)) then begin
     states = ustat
     return
  endif
  
  idx = 0
  if(ns gt 1 and (n_elements(statid) eq 0)) then begin
     print, 'Found states:'
     for ii=0, ns-1 do begin
        print, string(ii, format='(I5)') + ' -> ' + ustat[ii]
     endfor
     read, idx, prompt='Please select state ID: '
  endif

  if(n_elements(statid) eq 1) then if(statid gt 0 AND statid lt ns) then idx = statid

  
  ustat = ustat[idx]
  print, 'Using -> '+ustat
  




  
  ;; Extract ustat from st

  st = red_quicklook_chromis_extract(st, ustat)
  uscan = st.scan[uniq(st.scan, sort(st.scan))]
  nt = n_elements(uscan) 
  print,'Found scans -> '+strcompress(string(nt),/remove)
  
  for ii = 0,nt-1 do begin

     
     ;; Get files in that scan
     
     idx = where(st.scan eq uscan[ii], ct)
     if(ct gt 1) then idx = idx[ct/2-1]
     print,string(13B), 'loading -> '+ file_basename(st.files[idx]),format='(A,A,$)'

     if(~keyword_set(rms_select)) then begin
       if(ii eq 0) then tmp = red_readdata(st.files[idx], /silent) $
       else tmp = red_readdata(st.files[idx], /silent, nslice=0)
     endif else tmp = rms_select(red_readdata(st.files[idx],/silent))
     
     ;; init some stuff
     
     if(ii eq 0) then begin
        dim = size(tmp, /dim)
        
        x0 = clip[0]
        x1 = dim[0] - 1 - clip[1]
        y0 = clip[2]
        y1 = dim[1] - 1 - clip[3]

        nx = x1 - x0 + 1
        ny = y1 - y0 + 1

        cub = fltarr(nx, ny, nt)
        me = fltarr(nt)
     endif
          

     cub[*,*,ii] = ((tmp-dark)*gain)[x0:x1, y0:y1]
     me[ii] = median(cub[*,*,ii])
  endfor
  print, ''

  if(nt gt 3 and ~keyword_set(no_normalize)) then begin
     mm = mean(me)
     cc = poly_fit(findgen(nt), me/mm, 2, yfit=yfit)
     for ii = 0, nt-1 do cub[*,*,ii] /= yfit[ii]*mm
  endif
  
  cub = bytscl(histo_opt(cub))
  





  
  ;; Prep movie file

  pid = strcompress(string(call_external('libc.so.6', 'getpid')),/remove)
  ofile = 'ql_'+times+'_'+ustat
  file_mkdir,'quicklook/'
  odir = 'tmp_'+pid+'/'
  file_mkdir, odir
  spawn, 'rm '+odir+'*.png'
  oofile = 'tempo_'+pid

  
  if(keyword_set(idlcompress)) then begin
      oVid = IDLffVideoWrite('quicklook/'+ofile+'.mp4')
      vidStream = oVid.AddVideoStream(nx, ny, fps, PRESET='medium', BIT_RATE=bitrate)
   endif




  
  ;; compress images


  thisDevice = !D.Name
  Set_Plot, 'Z'
  Device, Set_Resolution=[nx, ny], Z_Buffer=0, Set_Pixel_Depth=24, Decomposed=0,set_font="Helvetica",/tt_font
  !P.font = 1
  
  k = 0L
  for ii=0, nt-1 do begin
     Erase
     loadct,0,/silent
     tv, cub[*,*,ii]
     
     ;loadct,3,/silent
     xyouts, [0.01], [0.97], 'scan='+string(ii,format='(I5)'), /normal, charsize=3., color=255, font=1

     snap = tvrd()
     tvlct,rr,gg,bb,/get

     
     if(keyword_set(idlcompress)) then !NULL=oVid.Put(vidStream, snap) $
     else write_png, odir+'img_'+string(k++, format='(I05)')+'.png', snap, rr, gg, bb
     
  endfor

 ; set_plot,'X'
  Device, decomposed=0

  if(~keyword_set(idlcompress)) then begin
     print, 'Saving -> '+'quicklook/'+ofile+'.mp4'
     mencode,'"mf://'+odir+'/*.png" -mf fps='+strcompress(string(fps), /remove), oofile, fps=fps, bitrate=bitrate
     spawn, 'rm '+odir+'*.png'
     spawn,'mv ./'+oofile+'.mp4 quicklook/'+ofile+'.mp4'
  endif
  
  oVid = 0
  print, ''
end
