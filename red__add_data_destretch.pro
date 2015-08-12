; docformat = 'rst'

;+
; :history:
; 
;   2015-04-01 : MGL. Use red_download to get turret log file.
; 
;-
pro red::add_data_destretch, scan = scan, min = min, max = max, smooth = smooth, $
                             bad = bad, nthreads=nthreads, nostretch = nostretch,$
                             scans_only = scans_only, no_cross_talk =no_cross_talk, $
                             mask = mask, overwrite = overwrite, extraclip = extraclip, $
                             t0 = t0, t1 = t1

  ;; Get procedure name
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0]) + ': '


  ;; Defaults
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(smooth) eq 0) then smooth = 3.0
  if(n_elements(nthreads) eq 0) then nthreads = 6L
  ;; Look for calibration data
  files = file_search(self.out_dir+'/calib_tseries/header.add.*idlsave', count = ct)
  
  if(ct eq 0) then begin
     print, 'ERROR, calibration data with destretch info not found'
     return
  endif else begin
     toread = 0L
     if(ct gt 1) then begin
        print, inam+'Found calibration files:'
        for ii = 0, ct-1 do print, ii, files[ii]
        read, toread, prompt='Select file ID: '
     endif

     restored_file = files[toread]

     print, inam + 'loading ', files[toread]
     headfile = files[toread]
     restore, headfile
  endelse
  tag = strjoin((strsplit(file_basename(headfile),'.',/extract))[1:3],'.')
  cfiles = file_search(self.out_dir+'/calib_tseries/tstep.'+tag+'.?????.idlsave', count = calcount)
  cnum = strarr(calcount)
  for ii = 0, calcount -1 do cnum[ii] = (strsplit(file_basename(cfiles[ii]), '.',/extract))[4]
  
  self->getcamtags, dir = dir
  cams = [self.camwbtag,self.camttag,self.camrtag]
  
  
  
  ;; Get darks and flats for all cameras
  for cc = 1, 2 do begin
     nam = self.out_dir + '/darks/'+cams[cc]+'.dark'
     if(~file_test(nam)) then begin
        print,inam+'ERROR, darks not found'
        return
     endif

     print, inam+'loading dark '+file_basename(nam)
     tmp = f0(nam)
     if cc eq 1 then begin
        dim = size(tmp,/dim)
        dd = fltarr(dim[0],dim[1],3)
        ff = fltarr(dim[0],dim[1],n_elements(uwav),3)
     endif
     dd[*,*,cc] = temporary(tmp)
     

     
     ;; load FF data
     for ww = 0, n_elements(uwav)-1 do begin
        if(cc eq 0) then dum = 0 $; ffile = self.out_dir+'/flats/'+cams[cc]+'.'+pref+'.flat' $
        else ffile = self.out_dir+'/flats/'+strjoin([cams[cc],pref,uwav[ww],'unpol.flat'],'.')

        if(~file_test(ffile)) then begin
           print, unam+' ERROR, flat not found -> '+ffile
        endif
        if(cc eq 0 ) then begin
           ;if(ww eq 0) then begin
           ;   print, inam+'loading flat '+file_basename(ffile);
;
;              ff[*,*,0,0] = red_flat2gain(f0(ffile), mi = min, ma = max, $
;                                          smooth = smooth, bad = bad, /preserve)
;              for tt = 1,  n_elements(uwav)-1 do ff[*,*,tt,cc] = ff[*,*,0,0]
;           endif
        endif else begin
           print, inam+'loading flat '+file_basename(ffile)
           ff[*,*,ww,cc] = red_flat2gain(f0(ffile), mi = min, ma = max, $
                                         smooth = smooth, bad = bad, /preserve)
        endelse
     endfor
  endfor ;; cc (cams)



  ;;
  ;; Prepare output file names and folders
  ;;
  time = (strsplit(file_basename(restored_file),'.',/extract))[3]
  odir = self.out_dir +'/summed_cubes/'+time+'/'
  odir1 = odir+'wb/'
  file_mkdir, odir
  file_mkdir, odir1
  
  ofiles_root = odir +strjoin( ['data',pref],'.')
  self->getcamtags
  wblen = strlen(self.camwbtag)

  data_fold = file_dirname(ifiles[0])
  data_root = strsplit(data_fold, '/',/extract) 
  dum = n_elements(data_root)
  data_root = '/'+strjoin(data_root[0:dum-2],'/') +'/'
  
  files_root = strmid(file_basename(ifiles), wblen, 500)

  ;; Per camera files
  tfiles = data_root + '/'+self.camt+'/'+self.camttag+files_root
  rfiles = data_root + '/'+self.camr+'/'+self.camrtag+files_root

  ;;
  ;; Load offset files
  ;;
  offiles = file_search(self.out_dir+'/calib/'+self.camttag+'.*'+pref+'*.xoffs', count = count)
  if(count eq 0) then begin
     print, inam+'ERROR, offset files not found'
     return
  endif
  ;; Use last one, likely to be continuum
  xto = offiles[count-1]
  yto = file_dirname(xto)+'/'+file_basename(xto,'xoffs')+'yoffs'
  xro = file_dirname(xto)+'/'+strjoin([self.camrtag,(strsplit(file_basename(xto),'.',/extract))[1:*]],'.')
  yro = file_dirname(xto)+'/'+strjoin([self.camrtag,(strsplit(file_basename(yto),'.',/extract))[1:*]],'.')
  
  xto = f0(xto)
  yto = f0(yto)
  xro = f0(xro)
  yro = f0(yro)

  if(nlc eq 4) then begin
     tloadfile = self.out_dir+'/polcal/'+self.camttag+'.'+pref+'.polcal.f0'
     rloadfile = self.out_dir+'/polcal/'+self.camrtag+'.'+pref+'.polcal.f0'
     if(file_test(tloadfile) AND file_test(rloadfile)) then begin
        print, inam+'loading polcal: '+file_basename(tloadfile)
        tmat = f0(tloadfile)
        print, inam+'loading polcal: '+file_basename(rloadfile)
        rmat = f0(rloadfile)
        pol = 1

        nx = abs(cl[0,1]-cl[1,1]) + 1
        ny = abs(cl[2,1]-cl[3,1]) + 1

        itmat = fltarr(4,4, nx, ny)
        irmat = fltarr(4,4, nx, ny)
        
        dim = size(tmat, /dim)
        tmat = reform(temporary(tmat), [4,4,dim[1], dim[2]])
        rmat = reform(temporary(rmat), [4,4,dim[1], dim[2]])

        print, inam +'inverting modulation matrix ... ', format='(A, $)'
        for yy = 0, dim[2]-1 do for xx=0,dim[1]-1 do tmat[*,*,xx,yy] = invert(tmat[*,*,xx,yy])
        for yy = 0, dim[2]-1 do for xx=0,dim[1]-1 do rmat[*,*,xx,yy] = invert(rmat[*,*,xx,yy])
        print, 'done'

        ;; Remove CCD tabs
        for jj = 0,3 do for ii = 0, 3 do begin
           tmat[ii,jj,*,*] = red_mask_ccd_tabs(reform(tmat[ii,jj,*,*]))
           rmat[ii,jj,*,*] = red_mask_ccd_tabs(reform(rmat[ii,jj,*,*]))
        end

        ;; Clip modulation matrix and apply offsets
        for jj = 0,3 do for ii = 0, 3 do begin
           itmat[ii,jj,*,*] = red_applyoffsets( red_clipim( reform( tmat[ii,jj,*,*]), cl[*,1]), xto, yto)
           irmat[ii,jj,*,*] = red_applyoffsets( red_clipim( reform( rmat[ii,jj,*,*]), cl[*,2]), xro, yro)
        endfor

        tmat = temporary(itmat)
        rmat = temporary(irmat)
        
     endif else begin
        print, inam + 'ERROR, polcal data not found in '+self.out_dir+'/polcal/'
        return
     endelse
  endif

  ;;
  ;;  telescope pointing log
  ;;
  if(self.isodate eq '') then begin
     isodate = ''
     read, isodate, prompt=inam+'please enter the ISO-date of the observations (YYYY-MM-DD): '
  endif else isodate = self.isodate
  date = strjoin(strsplit(isodate,'-./', /extract),'.')
  
  ;; Download Turret position log file if needed and read azel angles
  ;; from it. MGL 2015-04-01
  red_download, date = isodate, /turret, pathturret = turretfile
  if turretfile then begin
     print, inam + 'Using SST position LOG -> ' + turretfile
     telpos = red_read_azel(turretfile, date)
  endif else begin
     print, inam + 'No Turret log file'
     stop
  endelse 
 
  ;;
  ;; Load prefilter transmission
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
     print, inam + ' : prefilter files not found!'
     return
  endelse
  
  if(file_test(rpfile) AND file_test(rpwfile)) then begin
     print, inam + 'Loading:'
     print, '  -> ' + file_basename(rpfile)
     rpref = f0(tpfile)
     print, '  -> ' + file_basename(rpwfile)
     rwav = f0(tpwfile)
  endif else begin
     print, inam + ' : prefilter files not found!'
     return
  endelse

  tempwav = udwav - double(pref)
  totpref = 2./(red_intepf(twav, tpref, tempwav) + red_intepf(rwav, rpref, tempwav))
  


  
  ;;
  ;; Loop scans
  ;;
  dim  = size(refs,/dim)
  nx = dim[0]
  ny = dim[1]
  nadd = max(numadd)
  
  nb = fltarr(nx, ny, nlc, nwav,2)
  wb = fltarr(nx, ny)

  iscale = total(total(refs,1),1) /  (float(nx)*ny) 
  iscale = mean(iscale) / iscale
  
  if(n_elements(t0) eq 0) then t0 = 0
  if(n_elements(t1) eq 0) then t1 = nscan - 1
  
  for tt = t0, t1 do begin
     
     print, ' '
     ;; Check that file does not exist
     ofile = ofiles_root+'.'+uscan[tt]+'.fcube'
     if(file_test(ofile) AND ~keyword_set(overwrite)) then begin
        print, inam+'Warning, skipping processing of existing file -> '+file_basename(ofile)+''
        continue
     endif
     pos = where(cnum eq uscan[tt], mycount)
     if(mycount ne 1) then continue
     print, inam+'loading '+cfiles[pos]
     restore, cfiles[pos]
     
     
     ;; load data
     nadd = 0L
     wb[*] = 0.0
     nb[*] = 0.0
     print, inam +'loading/correcting data for scan ',tt,' -> t='+times[tt], format='(A,I0,A,$)'

     aver_shift = [0.,0.]
     aver_corr = corrs[*,*,*,0,0,0] * 0.
     
     for ww = 0, nwav-1 do begin
        for ll = 0, nlc-1 do begin
           for nn = 0, numadd[ww,ll,tt]-1 do begin
              ;;
              ;; Read, dark, flat, clip, fillpix
              ;;
              ;fzread, iwb,  ifiles[nn,ww,ll,tt], h
            
              
            ;  iwb = $
            ;     red_fillpix(red_clipim((float(iwb) - dd[*,*,0]) * ff[*,*,ww,0], cl[*,0]), nt=nthreads)
              inbt= $
                 red_fillpix(red_clipim((f0(tfiles[nn,ww,ll,tt]) - dd[*,*,1]) * ff[*,*,ww,1], cl[*,1]), nt=nthreads)
              inbr = $
                 red_fillpix(red_clipim((f0(rfiles[nn,ww,ll,tt]) - dd[*,*,2]) * ff[*,*,ww,2], cl[*,2]), nt=nthreads)
              nadd++

              
              ;;
              ;; Apply offsets to nb images
              ;;
              inbt = red_applyoffsets(temporary(inbt), xto, yto)
              inbr = red_applyoffsets(temporary(inbr), xro, yro)
              

              ;;
              ;; ishift
              ;;
              ishift = shifts[*,nn,ww,ll]
              ;if(n_elements(shifts_refs) gt 0) then ishift += shifts_refs[*,tt]

              
              ;;
              ;; Rotation? and shift
              ;;
              if(rot) then begin
                ; iwb =  red_rotation(temporary(iwb),  ang[tt], float(ishift[0]), float(ishift[1]))
                 inbt = red_rotation(temporary(inbt), ang[tt], float(ishift[0]), float(ishift[1]))
                 inbr = red_rotation(temporary(inbr), ang[tt], float(ishift[0]), float(ishift[1]))
              endif else begin
                ; iwb = shift(temporary(iwb),  ishift)
                 inbt= shift(temporary(inbt), ishift)
                 inbr= shift(temporary(inbr), ishift)
              endelse
              stop
              ;;
              ;; Apply distortion correction?
              ;;
              icorr = corrs[*,*,*,nn,ww,ll]
             ; if(n_elements(tgrid) gt 0) then icorr += reform(tgrid[tt,*,*,*])
              
              if(~keyword_set(nostretch)) then begin
                ; iwb = red_stretch(temporary(iwb),  icorr)
                 inbt= red_stretch(temporary(inbt), icorr)
                 inbr= red_stretch(temporary(inbr), icorr)
              endif

              
              ;;
              ;; add to state
              ;;
           ;   wb += iwb
              nb[*,*,ll,ww,0] += temporary(inbt)
              nb[*,*,ll,ww,1] += temporary(inbr) 


              ;;
              ;; Add shift and corrs
              ;;
              aver_shift +=  ishift
              aver_corr +=  icorr
              
           endfor
           nb[*,*,ll,ww,*] *= (iscale[tt] / float(numadd[ww,ll,tt]))
        endfor
     endfor
     print, ' ... done'
    ; wb /= float(nadd)
     writefits, odir1+'wb.'+pref+'.'+uscan[tt]+'.fits', wb, times[tt]
     
     ;;
     ;; Distort demod matrix with average im shift and average distortion
     ;;
     aver_shift /= float(nadd)
     aver_corr /= float(nadd)


     ;time /= nadd
     itime = times[tt];red_time2double(time, /dir)
     print, inam + 't_aver -> '+itime

     ;;
     ;; Load telescope matrix
     ;;
     mtel = red_telmat(pref, telpos, itime, /no_zero)
     imtel = invert(mtel)
     imtel /= imtel[0]

     
     
     ;;
     ;; prepare the demodulation matrix
     ;;
     itmat = tmat
     irmat = rmat
     
     for jj = 0,3 do for ii = 0, 3 do begin
        ;; rotation
        if(rot) then begin
           itmat[ii,jj,*,*] = red_rotation(reform(itmat[ii,jj,*,*]), ang[tt],  float(aver_shift[0]),  float(aver_shift[1]))
           irmat[ii,jj,*,*] = red_rotation(reform(irmat[ii,jj,*,*]), ang[tt],  float(aver_shift[0]),  float(aver_shift[1]))
        endif else begin
           
           ;; shift
           itmat[ii,jj,*,*] =   shift( reform(itmat[ii,jj,*,*]), aver_shift)
           irmat[ii,jj,*,*] =   shift( reform(irmat[ii,jj,*,*]), aver_shift)
        endelse
        
        ;; destretch?
        if(~keyword_set(nostretch)) then begin
           itmat[ii,jj,*,*] =   red_stretch( reform(itmat[ii,jj,*,*]), aver_corr)
           irmat[ii,jj,*,*] =   red_stretch( reform(irmat[ii,jj,*,*]), aver_corr)
        endif
        
     endfor
     
     ;;
     ;; Demodulate states
     ;;
     for ww = 0, nwav-1 do begin
        nb[*,*,*,ww,0] = red_demodulate_simple(itmat, nb[*,*,0,ww,0], nb[*,*,1,ww,0], nb[*,*,2,ww,0], nb[*,*,3,ww,0])
        nb[*,*,*,ww,1] = red_demodulate_simple(irmat, nb[*,*,0,ww,1], nb[*,*,1,ww,1], nb[*,*,2,ww,1], nb[*,*,3,ww,1])


        tmp = reform(nb[*,*,*,ww,*])
        nb[*,*,*,ww,*] = 0.0
        FOR j=0, 3 DO FOR i=0, 3 DO begin
           nb[*, *, j,ww, 0] += tmp[*, *, i, 0] * imtel[i, j]
           nb[*, *, j,ww, 1] += tmp[*, *, i, 1] * imtel[i, j]
        endfor
        tmp = 0B
     endfor


     ;; Combine beams
     tmean = median(nb[*,*,0,*,0])
     rmean = median(nb[*,*,0,*,1])
     bmean = 0.5*(tmean + rmean)
     
     t_scale = 0.5 * bmean / tmean
     r_scale = 0.5 * bmean / rmean
     
     cub = (nb[*,*,*,*,0] * t_scale + nb[*,*,*,*,1] * r_scale)
     for ww = 0, nwav-1 do cub[*,*,*,ww] *= totpref[ww]


     ;; Mask spectra and FOV for crosstalk?
     if(~keyword_set(no_cross_talk)) then begin
        if(keyword_set(mask) AND tt eq 0) then begin
           ppc = red_select_spoints(udwav, total(total(reform(cub[*,*,0,*]),1)/nx,1)/ny)
        endif else ppc = indgen(nwav)
        crt = red_get_ctalk(cub, idx=ppc)
        for zz=1,3 do for ww = 0, nwav-1 do cub[*,*,zz,ww] -= crt[zz]*cub[*,*,0,ww]
     endif

     nx1 = nx
     ny1 = ny

     if(n_elements(extraclip) eq 4) then begin
        x0 = extraclip[0]
        x1 = nx - extraclip[1] - 1
        y0 = extraclip[2]
        y1 = ny - extraclip[3] - 1
        nx1 = x1 - x0 + 1
        ny1 = y1 - y0 + 1
        cub = cub[x0:x1,y0:y1,*,*]
     endif

     
     print, inam+'saving '+ofile
     head = red_pol_lpheader(nx1, ny1, nwav*4L, /float)
     openw, lun, ofile, /get_lun
     writeu,lun, head
     writeu,lun, float(transpose(cub,[0,1,3,2]))
     free_lun, lun
     if(tt eq 0) then fzwrite, tempwav, file_dirname(ofile)+'/wav.'+pref+'.f0',' '
     
  endfor
  
  
  
  stop
  return
end
