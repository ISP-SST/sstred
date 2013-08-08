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
;    scan  : 
;   
;   
;   
;    cam  : 
;   
;   
;   
;    pref  : 
;   
;   
;   
;    debug  : 
;   
;   
;   
;    lc  : 
;   
;   
;   
;    psf  : 
;   
;   
;   
;    nthreads : 
;   
;   
;   
; 
; 
; :history:
;
;   2013-07-01: JdlCR : Created!
; 
;   2013-07-24 : Use red_show rather than show.
;
;
;-
pro red::sum_data_intdif, cam = cam, t1 = t1, nthreads = nthreads, pref = pref, $
                          verbose = verbose, overwrite=overwrite, $
                          min = min, max = ma, smooth = smooth, bad = bad, $
                          descatter = descatter, show=show

  inam = 'red::sum_data_intdif : '

  outdir = self.out_dir + '/cmap_intdif/'
  ucam = [self.camr, self.camt]


;
; Search files
;
  if self.ndir gt 1 then begin
     for ii = 0, self.ndir-1 do print, stri(ii)+$
                                       ' -> '+file_basename(self.data_list[ii])
     idx = 0L
     read, idx, prom = inam + 'Please select folder : '
     print, inam + 'Using -> '+self.data_list[idx]
     dir = self.data_list[idx]
  endif else dir = self.data_list[0]
  imdir = strsplit(dir,'/',/extract)
  outdir += '/'+imdir[n_elements(imdir)-1]+'/'
  file_mkdir, outdir

  spawn, 'find ' + dir + '/' + self.camr + '/ | grep camX | grep -v ".lcd."', files
  nf = n_elements(files)
  if(nf eq 0) then begin
     print, inam + 'ERROR, data folder is empty'
     if(debug) then stop else return
  endif

;
; Extract tags from file names
; 
  state = red_getstates(files)


;
; build my states
;
  mpref = state.pref
  upref = reform([mpref[uniq(mpref, sort(mpref))]])
  np = n_elements(upref)


  sel_pref = 1B
  if(keyword_set(pref)) then begin
     pos = where(upref eq pref, count)
     if(count eq 0) then begin
        print, inam + 'User supplied prefilter is wrong -> '+ pref
     endif else sel_pref = 0B
     
  endif 

  if(sel_pref) then begin
     
     if(np eq 1) then begin
        ip = 0 
     endif else begin
        print, inam + 'Found prefilters:'
        for ii = 0, np - 1 do print, string(ii, format='(I3)') + ' -> ' + upref[ii]
        ip = 0L
        read, ip, prompt = 'Select state id: '
     endelse

     pref = upref[ip]
  endif
  print, inam + 'Selected prefilter -> ' + pref


;
; Isolate files from prefilter
;
  idx = where(mpref eq pref, count)
  mwav = state.wav[idx]
  mdwav = state.dwav[idx]
  mlc = state.lc[idx]
  mscan = state.scan[idx]
  mnum = state.nums[idx]
  mfiles = state.files[idx]
; 
  uscan = mscan[uniq(mscan, sort(mscan))]
  ulc = mlc[uniq(mlc, sort(mlc))]
  uwav = mwav[uniq(mwav, sort(mdwav))]
  udwav = float(mdwav[uniq(mdwav, sort(mdwav))] - double(pref))
  ns = n_elements(uscan)
  nw = n_elements(uwav)
  nlc = n_elements(ulc)

;
; Get camera tags
;
  self -> getcamtags, dir = self.data_dir
  ctag = [self.camrtag, self.camttag]

;
; Start summing 
;
  ;; if(n_elements(t0) ne 0) then begin
  ;;    lscan = long(uscan)
  ;;    tt0 = where(lscan eq t0, count)
  ;;    if(count eq 0) then begin
  ;;       print, inam + 'Scan',t0,' not found -> starting at t0 = 0.'
  ;;       tt0 = 0L
  ;;    endif
  ;; endif else tt0 = 0L
  if(n_elements(t1) ne 0) then begin
     lscan = long(uscan)
     tt1 = where(lscan eq t1, count)
     if(count eq 0) then begin
        print, inam + 'Scan',t1,' not found -> ending at t1 = ',ns-1
        tt1 = ns-1L
     endif
  endif else tt1 = ns-1

  
  for cc = 0, 1 do begin
     if(keyword_set(cam)) then begin
        if ucam[cc] ne cam then begin
           print, inam + 'skipping cam -> '+ucam[cc]+' != '+cam
           continue
        endif
     end
     
     ;;
     ;; Descatter?
     ;;
     if(keyword_set(descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
        pf = self.descatter_dir+ '/' + ctag[cc] + '.psf.f0'
        bf = self.descatter_dir+ '/' + ctag[cc] + '.backgain.f0'
        if(~file_test(pf) OR ~file_test(bf)) then begin
           print, inam + 'ERROR, descattering data not found, returning.'
           return
        endif 
        pff = f0(pf)
        bff = f0(bf)
     endif


     tempdir=dir + '/' + ucam[cc]+'/'
     longi = strlen(file_dirname(mfiles[0])+'/'+self.camrtag)
     mmfiles = tempdir + ctag[cc]+strmid(mfiles,longi,200)
     cfile = outdir + '/' + ctag[cc] + '.'+ pref + '.intdif.icube'
     dfile = outdir + '/' + ctag[cc] + '.'+ pref + '.intdif.save'
     
     if(file_test(dfile) and ~keyword_set(overwrite)) then begin
        print, inam + 'WARNING! Data files already exist:'
        print, dfile
        print, cfile
        print, inam + 'you can overwrite them with /overwrite, otherwise we will resume the summing!'
        restore, dfile
        idx = where(done eq 0, count)
        if count eq 0 then begin
           print, inam + 'nothing to do for '+ctag[cc]
           continue
        endif else tt0 = idx[0]
        openu, lun, cfile, /get_lun
     endif else begin
        done = bytarr(ns)
        openw, lun, cfile, /get_lun
        tt0 = 0L
     endelse
     dim = size(f0(mmfiles[0]),/dim)
     head = red_pol_lpheader(dim[0], dim[1], n_elements(ulc)*(tt1+1L)*n_elements(uwav))
     writeu,lun, head
     dat = assoc(lun, intarr(dim[0],dim[1],/nozer), 512)
     nx = dim[0]
     ny = dim[1]
                                ;
                                ; load dark file
                                ;
     df = self.out_dir+'/darks/'+ctag[cc]+'.dark'
     if(file_test(df)) then dd = f0(df) else begin
        print, inam + 'ERROR, dark-file not found -> '+df
        stop
     endelse


                                ;
                                ; Load fitgains_ng results
                                ;
     cmf = self.out_dir + '/flats/spectral_flats/' + ctag[cc] + '.' +$
           pref + '.' + 'fit_results.sav'

     
                                ;
                                ; Loop wavelengths within the scan
                                ;
     for ss = tt0[0], tt1[0] do begin
        
                                ;
        for ll = 0L, nlc - 1 do begin

                                ;
           scstate = strarr(nw)
                                ;
           for ww = 0L, nw - 1 do begin


              istate = strjoin([uscan[ss], pref,uwav[ww], ulc[ll]],'.')
              scstate[ww] = istate
              
              idx = where((mwav EQ uwav[ww]) AND (mlc EQ ulc[ll]) AND $
                          (mscan EQ uscan[ss]), count)
              
              
              
              if(count eq 0) then begin
                 print, inam + 'no files found for state -> ' + istate
              endif else begin
                 print, inam + 'Found ' + red_stri(count) + ' images -> '+istate 
              endelse
              
              gf = self.out_dir+'/flats/'+strjoin([ctag[cc], pref, uwav[ww], 'unpol.flat'], '.')
              
                                ;
                                ; Load flat
                                ;
              if(file_test(gf)) then g =  red_flat2gain(f0(gf),/preserve, mi=mi,ma=ma,smooth=smooth,bad=bad) else begin
                 print, inam + 'ERROR, gain file not found -> ', gf
                 stop
              endelse
              
              if(keyword_set(verbose)) then print, transpose(file_basename(mmfiles[idx]))
              tmp = red_sumfiles(mmfiles[idx], /check) - dd

              if(keyword_set(descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
                 tmp = red_cdescatter(temporary(tmp), bff, pff, /verbose, nthreads = nthreads)
              endif
              
              tmp = red_fillpix((temporary(tmp)*g), $
                                nthreads=nthreads)

            
              ele = ss*nlc*nw + ll*nw + ww
              dat[ele] = fix(round(7.0 * tmp))

              if(keyword_set(show)) then begin
                 if n_elements(mydum) eq 0 then begin
                    mydum = 1
                    red_show, histo_opt(tmp)
                 endif
                 red_show, histo_opt(tmp),/now
              endif
           endfor
        endfor
        done[ss] = 1B
        save, file=dfile, done, uwav, ulc, uscan, nw, nlc, ns, pref, nx, ny, udwav
     endfor
     free_lun, lun
  endfor

  save, file=dfile, done, uwav, ulc, uscan, nw, nlc, ns, pref, nx, ny, udwav

end

