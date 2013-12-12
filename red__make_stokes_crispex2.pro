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
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
; 
; 
;-
pro red::make_stokes_crispex2, rot_dir = rot_dir, square = square
  inam = 'red::make_stokes_crispex : '
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


  f += '/cfg/results/stokes/'
  
  ;;
  ;; Look for time-series calib file
  ;;
  cfile = self.out_dir + '/calib_tseries/tseries.'+pref+'.'+time_stamp+'.calib.sav'
  if(~file_test(cfile)) then begin
     print, inam + 'Could not find calibration file: ' + cfile
     print, inam + 'Try to execute red::polish_tseries on this dataset first!'
     return
  endif else print, inam + 'Loading calibration file -> '+file_basename(cfile)
  restore, cfile
  tmean = mean(tmean) / tmean


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

  files = file_search(f+'/'+self.camttag+'.?????.'+pref+'.????_*', count=tf)
  st = red_get_stkstates(files)
  
  nwav = st.nwav
  nscan = st.nscan
  wav = st.uiwav * 1.e-3

  ;;
  ;; Interpolate prefilters to the observed grid 
  ;;
  pref = 2.0 / (red_intepf(twav, tpref, wav) +  1./red_intepf(rwav, rpref, wav))
  stop
  ;;
  ;; Load clean flats and gains
  ;;
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
     tmp[idx] = mean(tmp[idx]) / tmp[idx]
     tmp[idx] = tmp[idx] / tmp1[idx]
     if(n_elements(idx1) gt 0) then tmp[idx1] = 0.0d0
     tratio[*,*,ii] = temporary(tmp)
     
     tmp = f0(rff)
     tmp1 = f0(rgg)
     idx = where(tmp1 gt 0.0001, count, complement = idx1)
     tmp[idx] = mean(tmp[idx]) / tmp[idx]
     tmp[idx] = tmp[idx] / tmp1[idx]
     if(n_elements(idx1) gt 0) then tmp[idx1] = 0.0d0
     rratio[*,*,ii] = temporary(tmp)
  endfor

  ;;
  ;; Load WB image and define the image border
  ;;
  tmp = red_mozaic(momfbd_read(wbfiles[0]))
  dimim = red_getborder(tmp, x0, x1, y0, y1, square=square)

  ;;
  ;; Create temporary cube and open output file
  ;;
  d = fltarr(dimim[0], dimim[1], nwav)
  
  head =  red_unpol_lpheader(dimim[0], dimim[1], nwav*nscan)
  
  if(n_elements(odir) eq 0) then odir = f + '/crispex/'
  file_mkdir, odir
  ofile = 'crispex.'+pref+'.'+time_stamp+'.time_corrected.icube'
  
  openw, lun, odir + '/' + ofile, /get_lun
  writeu, lun, head
  point_lun, lun, 0L
  dat = assoc(lun, intarr(dimim[0], dimim[1], nwav,/nozero), 512)
  
  ;;
  ;; start processing data
  ;; 
  for ss = 0L, nscan-1 do begin
     wb = (red_mozaic(momfbd_read(wbfiles[ss])))[x0:x1, y0:y1]
     
     for ww = 0L, nwav - 1 do begin 
        state = strjoin((strsplit(file_basename(st.ofiles[ww,ss]),'.',/extract))[1:*],'.')
        
        ttf = f + '/' + self.camttag+'.'+state
        rrf = f + '/' + self.camrtag+'.'+state
        wwf = f + '/' + self.camwbtag+'.'+state
        print, inam + 'processing state -> '+state 
        ;;
        ;; Load flats used and compute ration with the clean one
        ;;
        tgain = self.out_dir + '/gaintables/'+self.camttag + '.' + st.ostate[ww,ss]


        ;;
        ;; Get destretch to anchor camera (residual seeing)
        ;;
        
        grid1 = dsgridnest(wb, (red_mozaic(momfbd_read(wwf)))[x0:x1, y0:y1], tiles, clips)
        
        tmp_raw0 = momfbd_read(ttf)
        tmp_raw1 = momfbd_read(rrf)

        ;;
        ;; Apply flat ratio after convolving with the PSF
        ;; of the patch: red_img2momfbd
        ;;
        tmp0 = (red_mozaic(tmp_raw0))[x0:x1, y0:y1]
        tmp1 = (red_mozaic(tmp_raw1))[x0:x1, y0:y1]
        trat = (red_mozaic(red_img2momfbd(tmp_raw0, tratio[*,*,ww])))[x0:x1, y0:y1]
        rrat = (red_mozaic(red_img2momfbd(tmp_raw1, rratio[*,*,ww])))[x0:x1, y0:y1]
        tmp0 = temporary(tmp0) * temporary(trat) * tpref[ww]
        tmp1 = temporary(tmp1) * temporary(rrat) * rpref[ww]

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
        d[*,*,ww] = rotate(stretch(red_rotation(temporary(tmp), ang[ss], total(shift[0,ss]), total(shift[1,ss])), reform(grid[ss,*,*,*])), rot_dir)
        
     endfor
     
     if(ss eq 0) then begin
        imean = fltarr(nwav)
        for ii = 0, nwav-1 do imean[ii] = median(d[*,*,ii])
        cscl = 15000./max(imean)
     endif
     
     dat[ss] = fix(round(d*cscl))
  endfor

  free_lun, lun
  print, inam + 'done'
  print, inam + 'result saved to -> '+odir+'/'+ofile 
  
end
