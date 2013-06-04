pro red::sumpinh, nthreads = nthreads, descatter = descatter, ustat = ustat ;, pref=pref
                                ;
  inam = 'red::sumpinh : '
                                ;
  if(~self.dopinh) then begin
     print, 'red::sumpinh : ERROR : undefined pinh_dir'
     return
  endif
  if(~keyword_set(nthreads)) then nthread = 2 else nthread = nthreads
                                ; if(n_elements(pref) eq 0) then pref = ''
                                ;
                                ;
                                ; create file list for each camera
  spawn, 'find ' + self.pinh_dir + '/' + self.camt + '/ | grep im.ex', tfiles
  spawn, 'find ' + self.pinh_dir + '/' + self.camr + '/ | grep im.ex', rfiles
  spawn, 'find ' + self.pinh_dir + '/' + self.camwb + '/ | grep im.ex', wfiles
  nt = n_elements(tfiles)
  nr = n_elements(rfiles)
  nw = n_elements(wfiles)
  

                                ;
  if((tfiles[0] eq '') OR (rfiles[0] eq '') OR (wfiles[0] eq '')) then begin
     print, 'red::sumpinh : files not found for all cameras:'
     print, '   '+self.camt+' -> '+ nt
     print, '   '+self.camr+' -> '+ nr
     print, '   '+self.camw+' -> '+ nw
     return
  endif
                                ;
                                ; sort files
  tfiles = red_sortfiles(tfiles)
  rfiles = red_sortfiles(rfiles)
  wfiles = red_sortfiles(wfiles)
                                ;
                                ; Get states for the cameras
  print, 'red::sumpinh : extracting states for all cameras' 
  tstat = red_getstates(tfiles)
  rstat = red_getstates(rfiles)
  wstat = red_getstates(wfiles)

                                ;
                                ; Flagging first image after tunning (will use tstat for the states)
                                ;
  
  state = tstat.state
  if(~keyword_set(ustat)) then ustat = state[uniq(state, sort(state))]
                                ;
  rstate = rstat.state
  wstate = wstat.state
  
                                ;
                                ; output sum

  outdir = self.out_dir+ '/' + 'pinh/'
  file_mkdir, outdir

                                ;
                                ; Load darks
                                ;
  tcam = red_camtag(tfiles[0])
  rcam = red_camtag(rfiles[0])
  wcam = red_camtag(wfiles[0])
  dt = f0(self.out_dir + '/darks/'+tcam+'.dark')
  dr = f0(self.out_dir + '/darks/'+rcam+'.dark')
  dw = f0(self.out_dir + '/darks/'+wcam+'.dark')


                                ;
                                ; loop
                                ;
  firsttime = 1B
  ns = n_elements(ustat)
  for ii = 0L, ns-1 do begin
     print, 'red::sumpinh : summing pinh for state -> ' + ustat[ii]
                                ;
     pos = where((state eq ustat[ii]), count)
     if(count GT 0) then tc = red_sumfiles(tstat.files[pos]) - dt else begin
        print, inam + 'WARNING-> No files found for the Transmitted camera -> '+ ustat[ii]
     endelse
     head = 'n_aver=' + red_stri(count)

     pos = where((rstate eq ustat[ii]) , count)
     if(count GT 0) then rc = red_sumfiles(rstat.files[pos]) - dr else begin
        print, inam + 'WARNING-> No files found for the Reflected camera -> '+ ustat[ii]
     endelse

     
     pos = where((wstate eq ustat[ii]) , count)
     if(count GT 0) then wc = red_sumfiles(wstat.files[pos]) - dw else begin
        print, inam + 'WARNING-> No files found for the WB camera -> '+ ustat[ii]
     endelse


     tflatf = self.out_dir + 'flats/' + strjoin([tcam, ustat[ii], 'flat'],'.')
     rflatf = self.out_dir + 'flats/' + strjoin([rcam, ustat[ii], 'flat'],'.')
     pref = (strsplit(ustat[ii], '.',/extract))[0]
     wflatf = self.out_dir + 'flats/' + strjoin([wcam, pref, 'flat'],'.')

     if(file_test(tflatf)) then begin
        tflat = f0(tflatf)
     endif else begin
        print, inam + 'ERROR -> flat not found for '+tcam+'.'+ustat[ii]
        stop
     endelse
     if(file_test(rflatf)) then begin
        rflat = f0(rflatf)
     endif else begin
        print, inam + 'ERROR -> flat not found for '+rcam+'.'+ustat[ii]
        stop
     endelse
     if(file_test(wflatf)) then begin
        wflat = f0(wflatf)
     endif else begin
        print, inam + 'ERROR -> flat not found for '+wcam+'.'+ustat[ii]
        stop
     endelse
     
                                ;
                                ; Descatter data?
                                ;

     if(keyword_set(descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
        if(firsttime) then begin
           ptf = self.descatter_dir+ '/' + tcam + '.psf.f0'
           prf = self.descatter_dir+ '/' + rcam + '.psf.f0'
           pwf = self.descatter_dir+ '/' + wcam + '.psf.f0'
           btf = self.descatter_dir+ '/' + tcam + '.backgain.f0'
           brf = self.descatter_dir+ '/' + rcam + '.backgain.f0'
           bwf = self.descatter_dir+ '/' + wcam + '.backgain.f0'
           
           if(file_test(ptf) and file_test(btf)) then begin
              Psft = f0(ptf)
              bgt = f0(btf)
              et = 1
           endif else etr = 0
           if(file_test(prf) and file_test(brf)) then begin
              Psfr = f0(prf)
              bgr = f0(brf)
              er = 1
           endif else er = 0
           if(file_test(pwf) and file_test(bwf)) then begin
              Psfw = f0(pwf)
              bgw = f0(bwf)
              ew = 1
           endif else ew = 0
           
           firsttime = 0B
        endif

        if(et) then begin
           tflat = red_cdescatter(tflat, bgt, psft, /verbose, nthreads = nthread)
           tc = red_cdescatter(tc, bgt, psft, /verbose, nthreads = nthread)
        endif
        if(er) then begin
           rflat = red_cdescatter(rflat, bgr, psfr, /verbose, nthreads = nthread)
           rc = red_cdescatter(rc, bgr, psfr, /verbose, nthreads = nthread)
        endif
        if(ew) then begin
           wflat = red_cdescatter(wflat, bgw, psfw, /verbose, nthreads = nthread)
           wc = red_cdescatter(wc, bgw, psfw, /verbose, nthreads = nthread)
        endif
     endif
     
     tg = red_flat2gain(tflat)
     rg = red_flat2gain(rflat)
     wg = red_flat2gain(wflat)


     print, inam + 'Applying dark + flat corrections to images'
     tc = red_fillpix(temporary(tc) * tg, mask=tg ne 0, nthreads = nthreads)
     rc = red_fillpix(temporary(rc) * rg, mask=rg ne 0, nthreads = nthreads)
     wc = red_fillpix(temporary(wc) * wg, mask=wg ne 0, nthreads = nthreads)
     

                                ;
                                ; Save
                                ;
     namout0 = tcam+'.' +ustat[ii]+'.pinh'
     namout1 = rcam+'.' +ustat[ii]+'.pinh'
     namout2 = wcam+'.' +ustat[ii]+'.pinh'
     print, 'red::sumpinh : saving ' + outdir + namout0
     fzwrite, fix(round(10. * tc)), outdir+namout0, head
                                ;
     print, 'red::sumpinh : saving ' + outdir + namout1
     fzwrite, fix(round(10. * rc)), outdir+namout1, head   
                                ;
     print, 'red::sumpinh : saving ' + outdir + namout2
     fzwrite, fix(round(10. * wc)), outdir+namout2, head
  endfor
end
