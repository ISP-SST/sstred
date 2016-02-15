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
;    no_descatter : in, optional, type=boolean 
;   
;      Don't do back-scatter compensation.
;   
;   
;   
; 
; 
; :History:
;
;   2013-07-01 : JdlCR : Created!
; 
;   2013-07-24 : Use red_show rather than show.
;
;   2014-01-02 : PS take out nremove (selection only done in prepmomfbd)
;                  use linked data, to skip incomplete scans
;
;
;   2015-05-05 : MGL. With a *, run all directories.
; 
;   2016-02-15 : MGL. Use loadbackscatter. Remove keyword descatter,
;                new keyword no_descatter.
; 
;
;-
pro red::sum_data_intdif, cam = cam, t1 = t1, nthreads = nthreads, pref = pref, $
                          verbose = verbose, overwrite=overwrite, $
                          no_descatter = no_descatter, show=show, LINK_DIR = link_dir

  inam = 'red::sum_data_intdif : '
  
  IF NOT keyword_set(link_dir) THEN link_dir = '/data'
  
;  outdir = self.out_dir + '/cmap_intdif/'
  ucam = [self.camr, self.camt]
  
    ;;; Search files.  Use links created by link_data

  dirs = file_search(self.out_dir + link_dir+'/*', COUNT = nd,/test_dir)
  IF nd EQ 0 THEN BEGIN
      print, inam + 'ERROR: No data found - did you run link_data?'
      if(debug) then stop else return
  ENDIF
  IF nd GT 1 THEN BEGIN
      FOR ii = 0, nd-1 DO print, red_stri(ii)+' -> '+ $
                                 file_basename(dirs[ii])

      idx = ''
      read, idx, prom = inam + 'Please select folder (* for all of them): '
      if idx ne '*' then begin
         dirs = dirs[long(idx)]
         
         print, inam + 'Using -> '+dirs
      endif
  ENDIF

  for idir = 0, n_elements(dirs)-1 do begin

     dir = dirs[idir]

                                ; Get camera tags
     self -> getcamtags, dir = dir
     ctag = [self.camrtag, self.camttag]
     
     ;outdir += file_basename(dir)+'/'
     outdir = self.out_dir + '/cmap_intdif/' + file_basename(dir)+'/'
     file_mkdir, outdir
     
     files = file_search(dir+ '/' + self.camr + '/cam*', count = nf)

     if(nf eq 0) then begin
        print, inam + 'ERROR, data folder is empty'
        if(debug) then stop else return
     endif
     
     files = red_sortfiles(temporary(files))

;
; Extract tags from file names
; 
     state = red_getstates(files, /LINKS)


                                ;
                                ; Remove frames
                                ;
                                ;red_flagtuning, state, remove
     

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
     mstar = state.star[idx]
; 
     uscan = mscan[uniq(mscan, sort(mscan))]
     ulc = mlc[uniq(mlc, sort(mlc))]
     uwav = mwav[uniq(mwav, sort(mdwav))]
     udwav = float(mdwav[uniq(mdwav, sort(mdwav))] - double(pref))
     ns = n_elements(uscan)
     nw = n_elements(uwav)
     nlc = n_elements(ulc)


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
        if(~keyword_set(no_descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
           self -> loadbackscatter, ctag[cc], pref, bff, pff
;           pf = self.descatter_dir+ '/' + ctag[cc] + '.psf.f0'
;           bf = self.descatter_dir+ '/' + ctag[cc] + '.backgain.f0'
;           if(~file_test(pf) OR ~file_test(bf)) then begin
;              print, inam + 'ERROR, descattering data not found, returning.'
;              return
;           endif 
;           pff = f0(pf)
;           bff = f0(bf)
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
                                ; Load fitgains results
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
                             (mscan EQ uscan[ss]) AND mstar eq 0, count)
                 
                 
                 
                 if(count eq 0) then begin
                    print, inam + 'no files found for state -> ' + istate
                 endif else begin
                    print, inam + 'Found ' + red_stri(count) + ' images -> '+istate 
                 endelse
                 
                 gf = self.out_dir+'/flats/'+strjoin([ctag[cc], pref, uwav[ww], 'unpol.flat'], '.')
                 
                                ;
                                ; Load flat
                                ;
                 ;;if(file_test(gf)) then g =
                 ;;red_flat2gain(f0(gf),/preserve,
                 ;;mi=mi,ma=ma,smooth=smooth,bad=bad) else begin
                 if(file_test(gf)) then begin
                    g = f0(gf)
                    g = median(g) / g
                    badpixels = where(~finite(g), count)
                    if(count gt 0) then g[badpixels] = 0.0
                 endif else begin
                    print, inam + 'ERROR, gain file not found -> ', gf
                    stop
                 endelse
                 
                 if(keyword_set(verbose)) then print, transpose(file_basename(mmfiles[idx]))
                 tmp = red_sumfiles(mmfiles[idx], /check) - dd

                 if(~keyword_set(no_descatter) AND self.dodescatter AND (pref eq '8542' OR pref eq '7772')) then begin
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

  endfor                        ; idir

end

