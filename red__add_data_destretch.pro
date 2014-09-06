pro red::add_data_destretch, pref = pref, scan = scan, min = min, max = max, smooth = smooth, bad = bad

  ;; Get procedure name
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0]) + ': '


  ;; Defaults
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(smooth) eq 0) then smooth = 3.0



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
  files = file_search(dir + self.camwb+'/cam*', count = nf)
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
     st = red_getstates(files, /links)
  endif
  if(n_elements(scan) gt 0) then begin
     pos = where(st.scan eq scan, count)
     if(count eq 0) then begin
        print, inam+'ERROR, scan '+scan+' does not exist'
        return
     endif
     files = files[pos]
     st = red_getstates(files, /links)
  endif
 


  ;; Get unique states
  uscan = st.scan[uniq(st.scan,sort(st.scan))]
  pos = uniq(st.wav,sort(st.dwav))
  uwav = st.wav[pos]

  ulc = st.lc[uniq(st.lc,sort(st.lc))]
  upref = st.pref[uniq(st.pref,sort(st.pref))]
  
  self.getcamtags, dir = dir
  cams = [self.camwbtag,self.camttag,self.camrtag]


  ;; Loop filters
  for pp = 0, n_elements(upref)-1 do begin
     ;; Get darks and flats
     for cc = 0, 2 do begin
        nam = self.out_dir + '/darks/'+cams[cc]+'.dark'
        if(~file_test(nam)) then begin
           print,inam+'ERROR, darks not found'
           return
        endif

        print, inam+'loading dark '+file_basename(nam)
        tmp = f0(nam)
        if cc eq 0 then begin
           dim = size(tmp,/dim)
           dd = fltarr(dim[0],dim[1],3)
           ff = fltarr(dim[0],dim[1],n_elements(uwav),3)
        endif
        dd[*,*,cc] = temporary(tmp)
        

        
        ;; load FF data
        for ww = 0, n_elements(uwav)-1 do begin
           if(cc eq 0) then ffile = self.out_dir+'/flats/'+cams[cc]+'.'+pref[pp]+'.flat' $
           else ffile = self.out_dir+'/flats/'+strjoin([cams[cc],pref[pp],uwav[ww],'unpol.flat'],'.')

           if(~file_test(ffile)) then begin
              print, unam+' ERROR, flat not found -> '+ffile
           endif
           print, inam+'loading flat '+ffile
           ff[*,*,ww,cc] = red_flat2gain(f0(ffile), mi = min, ma = max, smooth = smooth, bad = bad, /preserve)
        endfor

        
        


     endfor ;; cc (cams)
  endfor ;; upref
  


  stop
  return
end
