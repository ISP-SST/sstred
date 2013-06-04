pro red::makegains, descatter = descatter, nthreads = nthreads, cam = cam, pref = pref, min = min, max = max, bad=bad, preserve=preserve, smoothsize = smoothsize
                                ;
  tosearch = self.out_dir+'/flats/*.flat'
  files = file_search(tosearch, count = ct)
                                ;
  if(ct eq 0) then begin
     print, 'red::makegains : No flats found in: ' + tosearch
  endif
                                ;
  firsttime = 1B
  for ii = 0L, ct -1 do begin
     tmp = strsplit(file_basename(files[ii]), '.', /extract)
     if(keyword_set(pref)) then begin
        if(tmp[1] ne pref) then begin
           print, 'red::makegains : skipping prefilter -> '+tmp[1]
           continue
        endif
     endif
     fzread, flat, files[ii], head
                                ;
                                ; Only one camera?
                                ;
     if n_elements(cam) ne 0 then if tmp[0] NE cam then continue
                                ;
     if(keyword_set(descatter)) then begin
        if((tmp[1] eq '8542' OR tmp[1] eq '7772') AND self.dodescatter) then begin
           psff = self.descatter_dir+'/'+tmp[0]+'.psf.f0'
           bgf = self.descatter_dir+'/'+tmp[0]+'.backgain.f0'
           if(file_test(psff) AND file_test(bgf)) then begin
              psf = f0(psff)
              bg = f0(bgf)
              flat = red_cdescatter(flat, bg, psf, nthreads = nthreads, verbose = 1)
           endif
        endif
     endif

     gain = red_flat2gain(flat, ma=max, mi=min, bad=bad, preserve=preserve, smoothsize=smoothsize)
     
     namout = file_basename(files[ii], '.flat')+'.gain'
     outdir = self.out_dir+'/gaintables/'
     h = head
                                ;
                                ; Output gaintable
     file_mkdir, outdir
     print, 'red::makegains : saving '+outdir+namout
     fzwrite, float(gain), outdir+namout, h
  endfor
                                ;
  return
end
