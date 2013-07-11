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
;    rot_dir  : 
;   
;   
;   
;    square  : 
;   
;   
;   
;    tiles : 
;   
;   
;   
;    clips : 
;   
;   
;   
;    scans_only  : 
;   
;   
;   
;    overwrite  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;
;   2013-07-05 : Allow to bypass flat-ratio operations
;
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
; 
; 
;-
pro red::make_unpol_crispex, rot_dir = rot_dir, square = square, tiles=tiles, clips=clips, scans_only = scans_only, overwrite = overwrite, noflats=noflats, iscan=iscan
  inam = 'red::make_unpol_crispex : '
  if(n_elements(rot_dir) eq 0) then rot_dir = 0B

 
  ;;
  ;; select folder
  ;;
  search = self.out_dir +'/momfbd/'
  f = file_search(search+'*', count = ct, /test_dir)
  if(ct eq 0) then begin
     print, inam + 'No sub-folders found in: ' + search
     return
  endif

  if(ct gt 1) then begin
     print, inam + 'Found folders(s): '
     for ii = 0L, ct-1 do print, red_stri(ii) + '  -> '+f[ii]
     idx = 0L
     read, idx, prompt = inam + 'Select folder ID: '
     idx = idx>0 < (ct-1)
     f = f[idx]
  endif

  print, inam + 'Selected -> '+ f
  time_stamp = strsplit(f, '/', /extract)
  time_stamp = time_stamp[n_elements(time_stamp)-1]

  ;;
  ;; Search prefilters in folder
  ;;
  search = f
  f = file_search(f+'/*', /test_dir, count = ct)
  if(ct eq 0) then begin
     print, inam + 'No sub-folders found in: ' + search
     return
  endif
  
  if(ct gt 1) then begin
     print, inam + 'Found prefilters(s): '
     for ii = 0L, ct-1 do print, red_stri(ii) + '  -> '+file_basename(f[ii])
     idx = 0L
     read, idx, prompt = inam + 'Select folder ID: '
     idx = idx>0 < (ct-1)
     f = f[idx]
  endif
  print, inam + 'Selected -> '+ f
  pref = strsplit(f, '/', /extract)
  pref = pref[n_elements(pref)-1]


  f += '/cfg/results/'
  
  ;;
  ;; Look for time-series calib file
  ;;
  if(~keyword_set(scans_only)) then begin
     cfile = self.out_dir + '/calib_tseries/tseries.'+pref+'.'+time_stamp+'.calib.sav'
     if(~file_test(cfile)) then begin
        print, inam + 'Could not find calibration file: ' + cfile
        print, inam + 'Try to execute red::polish_tseries on this dataset first!'
        return
     endif else print, inam + 'Loading calibration file -> '+file_basename(cfile)
     restore, cfile
     tmean = mean(tmean) / tmean
  endif else tmean = replicate(1.0, 10000) ; Dummy time correction

  ;;
  ;; Camera tags
  ;;
  self->getcamtags, dir = self.data_dir


  ;;
  ;; Load prefilter
  ;;
  tpfile = self.out_dir + '/prefilter_fits/'+self.camttag+'.'+pref+'.prefilter.f0'
  tpwfile = self.out_dir + '/prefilter_fits/'+self.camttag+'.'+pref+'.prefilter_wav.f0'
  rpfile = self.out_dir + '/prefilter_fits/'+self.camrtag+'.'+pref+'.prefilter.f0'
  rpwfile = self.out_dir + '/prefilter_fits/'+self.camrtag+'.'+pref+'.prefilter_wav.f0'

  if(file_test(tpfile) AND file_test(tpwfile)) then begin
     print, inam + 'Loading:'
     print, '  -> ' + file_basename(tpfile)
     tpref = f0(tpfile)
     print, '  -> ' + file_basename(tpwfile)
     twav = f0(tpwfile)
  endif else begin
     print, inam + 'prefilter files not found!'
     return
  endelse

  if(file_test(rpfile) AND file_test(rpwfile)) then begin
     print, inam + 'Loading:'
     print, '  -> ' + file_basename(rpfile)
     rpref = f0(tpfile)
     print, '  -> ' + file_basename(rpwfile)
     rwav = f0(tpwfile)
  endif else begin
     print, inam + 'prefilter files not found!'
     return
  endelse

  tfiles = file_search(f+'/'+self.camttag+'.?????.'+pref+'.????_*momfbd', count=tf)
  rfiles = file_search(f+'/'+self.camrtag+'.?????.'+pref+'.????_*momfbd', count=rf)
  wfiles = file_search(f+'/'+self.camwbtag+'.?????.'+pref+'.????_*momfbd', count=wf)
  wbfiles = file_search(f+'/'+self.camwbtag+'.?????.'+pref+'.momfbd', count = wbf)

  st = red_get_stkstates(tfiles)

  if((tf NE rf) OR (tf NE wf)) then begin
     print, inam + 'ERROR -> Irregular number of images:'
     print, '  '+self.camttag+' -> '+red_stri(tf)
     print, '  '+self.camrtag+' -> '+red_stri(rf)
     print, '  '+self.camwbtag+' -> '+red_stri(wf)
     return
  endif
  
  nwav = st.nwav
  nscan = st.nscan
  wav = st.uiwav * 1.e-3

  ;;
  ;; Interpolate prefilters to the observed grid 
  ;;
  tpref = 1./red_intepf(twav, tpref, wav)
  rpref = 1./red_intepf(rwav, rpref, wav)

  ;;
  ;; Load clean flats and gains
  ;;
  if(~keyword_set(noflats)) then begin
     for ii = 0, nwav-1 do begin
        
        tff = self.out_dir + 'flats/'+self.camttag + '.'+pref+'.'+st.uwav[ii]+'.unpol.flat'
        rff = self.out_dir + 'flats/'+self.camrtag + '.'+pref+'.'+st.uwav[ii]+'.unpol.flat'
        tgg = self.out_dir + 'gaintables/'+self.camttag + '.'+pref+'.'+st.uwav[ii]+'.lc4.gain'
        rgg = self.out_dir + 'gaintables/'+self.camrtag + '.'+pref+'.'+st.uwav[ii]+'.lc4.gain'
        if(ii eq 0) then print, inam + 'Loading: '
        print,' -> '+tff
        print,' -> '+rff
        print,' -> '+tgg
        print,' -> '+rgg
        
        if(~file_test(tff) OR ~file_test(rff) OR ~file_test(tgg) OR ~file_test(rgg)) then begin
           print, inam + 'ERROR -> Flat/gain files not found'
           return
        endif
        
        if(ii eq 0) then begin
           dim = size(f0(tff),/dimen)
           tratio = fltarr(dim[0], dim[1], nwav)
           rratio = fltarr(dim[0], dim[1], nwav)
        endif 
        
        tmp = f0(tff)
        tmp1 = f0(tgg)
        idx = where(tmp1 gt 0.0001, count, complement = idx1)
        mask = bytarr(dim) + 1B
        tmp[idx] = mean(tmp[idx]) / tmp[idx]
        tmp[idx] = tmp[idx] / tmp1[idx]
        if(n_elements(idx1) gt 0) then begin
           tmp[idx1] = 0.0d0
           mask[idx1] = 0B
        endif
        tmp = red_fillpix(tmp, mask=red_cleanmask(mask),nthreads=6)
        tratio[*,*,ii] = temporary(tmp)
        
        tmp = f0(rff)
        tmp1 = f0(rgg)
        idx = where(tmp1 gt 0.0001, count, complement = idx1)
        mask = bytarr(dim) + 1B
        tmp[idx] = mean(tmp[idx]) / tmp[idx]
        tmp[idx] = tmp[idx] / tmp1[idx]
        if(n_elements(idx1) gt 0) then begin
           tmp[idx1] = 0.0d0
           mask[idx1] = 0B
        endif
        tmp = red_fillpix(tmp, mask = red_cleanmask(mask),nthreads=6)
        rratio[*,*,ii] = temporary(tmp)
     endfor
  endif


  ;;
  ;; Load WB image and define the image border
  ;;
  tmp = red_mozaic(momfbd_read(wbfiles[0]))
  dimim = red_getborder(tmp, x0, x1, y0, y1, square=square)

  ;;
  ;; Create temporary cube and open output file
  ;;
  d = fltarr(dimim[0], dimim[1], nwav)  
  if(~keyword_set(scans_only)) then begin
     head =  red_unpol_lpheader(dimim[0], dimim[1], nwav*nscan)
  endif else begin
     head = red_unpol_lpheader(dimim[0], dimim[1], nwav)
  endelse
  
  if(n_elements(odir) eq 0) then odir = self.out_dir + '/crispex/' + time_stamp + '/'
  file_mkdir, odir

  if(~keyword_set(scans_only)) then begin
     ;; Open assoc file for output of multi-scan data cube.
     ofile = 'crispex.'+pref+'.'+time_stamp+'.time_corrected.icube'

     if file_test(odir + '/' + ofile) then begin
        if keyword_set(overwrite) then begin
           print, 'Overwriting existing data cube:'
           print, odir + '/' + ofile
        endif else begin
           print, 'This data cube exists already:'
           print, odir + '/' + ofile
           return
        endelse
     endif
     
     openw, lun, odir + '/' + ofile, /get_lun
     writeu, lun, head
     point_lun, lun, 0L
     dat = assoc(lun, intarr(dimim[0], dimim[1], nwav,/nozero), 512)
  endif 
  ;;
  ;; start processing data
  ;; 

  if(~keyword_set(tiles) OR (~keyword_set(clips))) then begin
     tiles = [8,16,32,64]
     clips = [8,4,2,1]
  endif

  for ss = 0L, nscan-1 do begin
     if(n_elements(iscan) gt 0) then if(iscan ne st.uscan[ss]) then continue
     print, inam + 'processing scan -> '+st.uscan[ss]

     IF(SS EQ 0) THEN BEGIN
        fzwrite, wav, odir + '/' + 'wav.'+pref+'.f0',' '
     endif

     if(keyword_set(scans_only)) then begin
        ofile = 'crispex.'+pref+'.'+time_stamp+'_scan='+st.uscan[ss]+'.icube'
        if file_test(odir + '/' + ofile) then begin
           if keyword_set(overwrite) then begin
              print, 'Overwriting existing data cube:'
              print, odir + '/' + ofile
           endif else begin
              print, 'Skip to next scan, this one exists already:'
              print, odir + '/' + ofile
              continue          ; Skip to next iteration of "for ss ..." loop.
           endelse
        endif
     endif

     wb = (red_mozaic(momfbd_read(wbfiles[ss])))[x0:x1, y0:y1]
     
     for ww = 0L, nwav - 1 do begin 
        state = strjoin((strsplit(file_basename(st.ofiles[ww,ss]),'.',/extract))[1:*],'.')
        
        ttf = f + '/' + self.camttag+'.'+state
        rrf = f + '/' + self.camrtag+'.'+state
        wwf = f + '/' + self.camwbtag+'.'+state
        print, inam + 'processing state -> '+state 
     

        ;;
        ;; Get destretch to anchor camera (residual seeing)
        ;;
        wwi = (red_mozaic(momfbd_read(wwf)))[x0:x1, y0:y1]
        grid1 = dsgridnest(wb, wwi, tiles, clips)
        print, 'computed grid'
        tmp_raw0 = momfbd_read(ttf)
        tmp_raw1 = momfbd_read(rrf)

        ;;
        ;; Apply flat ratio after convolving with the PSF
        ;; of the patch: red_img2momfbd
        ;;
        tmp0 = (red_mozaic(tmp_raw0))[x0:x1, y0:y1] * tpref[ww]
        tmp1 = (red_mozaic(tmp_raw1))[x0:x1, y0:y1] * rpref[ww]
        if(~keyword_set(noflats)) then begin
           trat = (red_mozaic(red_img2momfbd(tmp_raw0, tratio[*,*,ww])))[x0:x1, y0:y1]
           rrat = (red_mozaic(red_img2momfbd(tmp_raw1, rratio[*,*,ww])))[x0:x1, y0:y1]
           
           tmp0 = temporary(tmp0) * temporary(trat) 
           tmp1 = temporary(tmp1) * temporary(rrat) 
        endif

        ;;
        ;; combine cameras, compute scale factor avoiding borders...
        ;;
        dim = size(tmp0,/dim)
        xx0 = round(dim[0] * 0.15)
        xx1 = round(dim[0] * 0.85)
        yy0 = round(dim[1] * 0.15)
        yy1 = round(dim[1] * 0.85)
        
        me = mean(tmp0[xx0:xx1,yy0:yy1] + tmp1[xx0:xx1,yy0:yy1]) * 0.5
        sclt = me / (mean(tmp0[xx0:xx1,yy0:yy1]))
        sclr = me / (mean(tmp1[xx0:xx1,yy0:yy1]))
        
        tmp = (temporary(tmp0) * sclt + temporary(tmp1) * sclr) * tmean[ss]
        
        ;;
        ;; Apply destretch to anchor camera and prefilter correction
        ;;

        tmp = stretch(temporary(tmp), grid1)
        
        ;;
        ;; Apply derot, align, dewarp
        ;;
        if(~keyword_set(scans_only)) then d[*,*,ww] = rotate(stretch(red_rotation(temporary(tmp), ang[ss], total(shift[0,ss]), total(shift[1,ss])), reform(grid[ss,*,*,*])), rot_dir) else d[*,*,ww] = rotate( temporary(tmp), rot_dir)
        
     endfor
     
     if(ss eq 0) then begin
        imean = fltarr(nwav)
        for ii = 0, nwav-1 do imean[ii] = median(d[*,*,ii])
        cscl = 15000./max(imean)
     endif
     
     if(~keyword_set(scans_only)) then begin
        ;; Write this scan's data cube to assoc file
        dat[ss] = fix(round(d*cscl))
     endif else begin
        ;; Write this scan's data cube as an individual file.
        print, inam + 'saving to '+ odir + '/' + ofile
        openw, lun, odir + '/' + ofile, /get_lun
        writeu, lun, head
        writeu, lun, fix(round(d))
        free_lun, lun
     endelse
  endfor
  
  if(~keyword_set(scans_only)) then begin
     ;; Close assoc file for output of multi-scan data cube.
     free_lun, lun
     print, inam + 'done'
     print, inam + 'result saved to -> '+odir+'/'+ofile 
  endif

end
