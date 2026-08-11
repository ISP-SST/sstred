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
;   2013-06-04 : Split from monolithic version of crispred.pro.
;   2013-06-27 : JdlCR: Fixed bug with cam_t, it used files from cam_r
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
; 
;-
pro red::make_cmap_intdif, scan = scan, cam = cam, pref = pref, debug = debug, lc = lc, psf = psf, nthreads=nthreads, verbose = verbose, overwrite=overwrite
;
  inam = 'red::make_cmap_intdif : '

  outdir = self.out_dir + '/cmap_intdif/'
  file_mkdir, outdir
  ucam = [self.camr, self.camt]

;
; Search files
;

  spawn, 'find ' + self.data_dir + '/' + self.camr + '/ | grep camX', files
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
  self -> getdetectors, dir = self.data_dir
  ctag = [self.camrtag, self.camttag]

;
; Start summing 
;


  for ss = 0L, ns - 1 do begin
                                ;
                                ; Check for user supplied scan number
                                ;
     if(keyword_set(scan)) then begin
        if(uscan[ss] ne string(scan, format='(I05)')) then begin
           print, inam + 'skipping scan -> ' + uscan[ss]
           continue
        endif
     endif

     for cc = 0, 1 do begin
        if(keyword_set(cam)) then begin
           if ucam[cc] ne cam then begin
              print, inam + 'skipping cam -> '+ucam[cc]+' != '+cam
              continue
           endif
        end
        tempdir=self.data_dir + '/' + ucam[cc]+'/'
        longi = strlen(file_dirname(mfiles[0])+'/'+self.camrtag)
        mmfiles = tempdir + ctag[cc]+strmid(mfiles,longi,200)
        
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
        cmf = self.out_dir + '/flats/spectral_flats/' + ctag[cc] + '.' + pref + '.' + 'fit_results.sav'
        if(file_test(cmf)) then begin
           print, inam + 'loading cavity-map from ' + cmf
           restore, cmf
           cmap = float(reform(fit.pars[1,*,*]))
           fit = 0B
        endif else begin
           print, inam + 'ERROR, cavity-map file not found -> '+cmf
           print, '  -> have you executed fitgains_ng?'
           stop
        endelse

        
        

        
                                ;
                                ; Loop wavelengths within the scan
                                ;
        for ll = 0L, nlc - 1 do begin
                                ;
           if(keyword_set(lc)) then begin
              if(ulc[ll] ne lc) then begin
                 print, inam + 'skipping lc = '+ulc[ll]+' != '+lc
                 continue
              endif
           endif
                                ;
           scstate = strarr(nw)
                                ;
           for ww = 0L, nw - 1 do begin
              istate = strjoin([uscan[ss], pref,uwav[ww], ulc[ll]],'.')
                                ;istate = strjoin([uscan[ss], pref,uwav[ww]],'.')

              scstate[ww] = istate
              
              idx = where((mwav EQ uwav[ww]) AND (mlc EQ ulc[ll]) AND $
                          (mscan EQ uscan[ss]), count)
              
                                ; idx = where((mwav EQ uwav[ww]) AND (mscan EQ uscan[ss]), count)
              

              if(count eq 0) then begin
                 print, inam + 'no files found for state -> ' + istate
              endif else begin
                 print, inam + 'Found ' + red_stri(count) + ' images -> '+istate 
              endelse
              
              gf = self.out_dir+'/flats/'+strjoin([ctag[cc], pref, uwav[ww], 'unpol.flat'], '.')
              
                                ;
                                ; Load flat and dark
                                ;
              if(file_test(gf)) then g =  self->flat2gain(f0(gf)) else begin
                 print, inam + 'ERROR, gain file not found -> ', gf
                 stop
              endelse
              
              if(keyword_set(verbose)) then print, transpose(file_basename(mmfiles[idx]))

              if(ww eq 0) then begin
                 tmp = red_fillpix(((red_sumfiles(mmfiles[idx], /check) - dd)*g), nthreads=nthreads)
                 dim = size(tmp, /dim)
                 cub = fltarr(nw,dim[0], dim[1])
                 cub[ww,*,*] = temporary(tmp)
              endif else cub[ww,*,*] = red_fillpix(((red_sumfiles(mmfiles[idx], /check) - dd)*g), nthreads=nthreads)
              
              
           endfor               ; kk
           cub2 = fltarr(nw, dim[0], dim[1])
           tmp = fltarr(nw)
           for jj = 0, dim[1] - 1 do begin
              for ii = 0, dim[0] - 1 do begin
                 f77_hintep, nw, float(udwav), cub[*,ii,jj], nw, float(udwav + cmap[ii,jj]), tmp
                 cub2[*,ii,jj] = tmp
              endfor
              print, string(13B)+inam+'Shifting cube to a constant wavelength grid -> ',100./(dim[1] - 1.0) * jj, '%', format='(A,F5.1,A,$)'
           endfor
           print, ' '
           
           rat = temporary(cub) / temporary(cub2)
           
           for ww =  0, nw-1 do begin
              if(keyword_set(psf)) then begin
                 rat[ww,*,*] = red_convolve(reform(rat[ww,*,*]), psf/total(psf))
              endif

              outname = outdir + '/' + ctag[cc] + '.' + scstate[ww] + '.intdiff'
              print, inam + 'saving -> '+ outname
              fcwrite, round(reform(rat[ww,*,*]*10000.)), outname, ' '
           endfor
           
        endfor                  ; ll
        
     endfor                     ; cc
  endfor                        ; ss
  
end
