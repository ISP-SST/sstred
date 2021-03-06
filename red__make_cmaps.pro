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
;    wbpsf  : 
;   
;   
;   
;    reflected  : 
;   
;   
;   
;    square : 
;   
;   
;   
;    rot_dir  : 
;   
;   
;   
;    fwhm  : 
;   
;   
;   
;    only_scans : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
pro red::make_cmaps,  wbpsf = wbpsf, reflected = reflected, square=square, rot_dir = rot_dir, $
                      fwhm = fwhm, only_scans=only_scans, remove_smallscale = remove_smallscale, $
                      wavelength_cube = wavelength_cube, float=float

  inam  = 'red::make_cmaps : '
;  if(n_elements(cmap) eq 0) then begin
;     print, inam + 'Error, provide a valid image input'
;     return
                                ; endif
  if(n_elements(rot_dir) eq 0) then rot_dir = 0
  if(~keyword_set(tiles) OR (~keyword_set(clip))) then begin
     itiles = [8,16,32,64]
     iclip = [8,4,4,2]
  endif
  if(n_elements(fwhm) eq 0) then fwhm = 7.0
  if(keyword_set(wavelength_cube)) then dowav=1B else dowav=0B

  ;;
  ;; Search directories
  ;;
  root = self.out_dir + 'momfbd/'
  dir = red_select_sub(root)
  if(~file_test(dir+'/cfg')) then dir = red_select_sub(dir)
  time = (strsplit(dir,'/',/extract))
  pref = time[n_elements(time)-1]
  time = time[n_elements(time)-2]


  ;;
  ;; Search files
  ;; 
  self->getdetectors
  
  cam = self.camttag
  if(keyword_set(reflected)) then cam = self.camrtag

  wbf = file_search(dir+'/cfg/results/'+self.camwbtag+'.?????.????.momfbd', count = cw)
  if(keyword_set( wavelength_cube)) then begin
     nbf = file_search(dir+'/cfg/results/'+cam+'.*.lc1.momfbd', count = cn)
     if(cn eq 0) then nbf = file_search(dir+'/cfg/results/'+cam+'.*.lc4.momfbd', count = cn)
     ;; wbff= file_search(dir+'/cfg/results/'+cam+'.?????.????.momfbd', count = cw1)
     ;;
     ;; States
     ;;
     st = red_get_stkstates(nbf)
     pref = (strsplit(st.state[0],'.',/extract))[1]
  endif


  ;;
  ;; Cavity map
  ;;
  cfile = self.out_dir+'flats/spectral_flats/'+cam+'.'+pref+'.fit_results.sav'
  if(~file_test(cfile)) then begin
     print, inam + 'Error, calibration file not found -> '+cfile
     stop
  endif
  restore, cfile
  cmap = reform(fit.pars[1,*,*])
  fit = 0B
  
  ;;
  ;; Small scale already corrected?
  if(keyword_set(remove_smallscale)) then begin
     npix = 30
     cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
     cpsf /= total(cpsf, /double)
     cmap = red_convolve(temporary(cmap), cpsf)
     ;cmap -= lscale
  endif
  

  
  ;;
  ;; time calibration data
  ;;
  time1 = time
  if(~keyword_set(only_scans)) then begin
     cfile = self.out_dir + 'calib_tseries/tseries.'+pref+'.'+time+'.calib.sav'
     if(~file_test(cfile)) then begin
        print, inam + 'Error: Could not find calibration file -> '+cfile
        stop
     endif
     restore, cfile
     full = 0
     if(n_elements(ff) eq 5) then full = 1
  endif

  ;;
  ;; Load clipping info for simple cmap
  ;; 
  cfile = self.out_dir + '/calib/align_clips.'+pref+'.sav'
  if(~file_test(cfile)) then begin
     print, inam + 'Error: Could not find calibration file -> '+cfile
     return
  endif
  print, inam +'loading -> '+cfile
  restore, cfile
  

  if(keyword_set(wavelength_cube)) then begin
     ;;
     ;; loop case wavelength dep.
     ;;
     nw = st.nwav
     nscan = st.nscan
     ;;
     ;; get image border
     ;;
  endif else nscan = n_elements(ang)
  if(n_elements(x0) eq 0) then dimim = red_getborder(red_mozaic(momfbd_read(wbf[0])), x0, x1, y0, y1, square=square)


  nx = x1 - x0 + 1L
  ny = y1 - y0 + 1L
  if(full) then begin
     nx = nd[0]
     ny = nd[1]
  endif
  
  ;;
  ;; load offsets file
  ;;
  cfile = file_search(self.out_dir + '/calib/'+cam+'.*'+pref+'.*.xoffs', count = ct)
  if(ct eq 0) then begin
     print, inam + 'Could not find offset files, skipping correction may result in small errors at the edges'
     xoffs = cmap *0.
     yoffs = cmap *0.
  endif else begin
     cfile = cfile[ct-1] ;; just pick one close to continuum (last state?)
     print, inam+'loading offset files:'
     print, '  -> '+file_basename(cfile)
     print, '  -> '+file_basename(cfile,'xoffs')+'yoffs'
     xoffs = f0(cfile)
     yoffs = f0(file_dirname(cfile)+'/'+file_basename(cfile,'xoffs')+'yoffs')
  endelse
     
  ;;
  ;; init cmap2 for simple case
  ;;
  idx = 1
  if(keyword_set(reflected)) then idx = 2
  if(~keyword_set(remove_smallscale)) then begin
     npsf = round(fwhm * 7.)
     if((npsf/2)*2 eq npsf) then npsf += 1L

     psf = red_get_psf(npsf, npsf, fwhm, fwhm)
     psf /= total(psf, /double)
     cmap1 = red_convolve(cmap, psf)
  endif else cmap1 = cmap
  cmap1 = (red_applyoffsets(red_clipim(temporary(cmap1), cl[*,idx]), xoffs,yoffs))[x0:x1,y0:y1]


  ;;
  ;; Open output files
  ;;
  odir = self.out_dir + '/crispex/cavity_maps/'
  file_mkdir, odir
  root1 = 'nx='+red_stri(nx)+'_ny='+red_stri(ny)+'_nt='+red_stri(nscan)
  if(keyword_set(wavelength_cube)) then root2 = 'nx='+red_stri(nx)+'_ny='+red_stri(ny)+'_nwav='+red_stri(nw)+'_nt='+red_stri(nscan)
  root3 = 'nx='+red_stri(nx)+'_ny='+red_stri(ny)

  if(keyword_set(float)) then exten='.fcube' else exten='.icube'
  
  file1 = odir + 'cmap.'+pref+'_'+time1+'.'+root1+'.simple'+exten
  if(dowav) then file2 = odir + 'cmap.'+pref+'_'+time1+'.'+root2+'.allwav'+exten


  
  if(keyword_set(only_scans)) then begin
     cmap11 = rotate(temporary(cmap1), rot_dir)
     ofil = odir+'cmap.'+pref+'.'+root3+'.f0'
     fzwrite, cmap11, ofil, ' '
     print, inam+'result saved to -> '+ofil
     return
  endif
  ;;
  ;; Lp format header
  ;;
  openw, lun1, file1, /get_lun
  if(dowav) then openw, lun2, file2, /get_lun
  writeu, lun1, red_unpol_lpheader(nx, ny, nscan,float=float)
  if(dowav) then writeu, lun2, red_unpol_lpheader(nx, ny, nscan*nw, float=float)

  for ss = 0L, nscan-1 do begin
     ;;
     ;; CASE 1
     ;;
     ;; Derotate and shift (case 1)
     if(full) then begin
        cmap11 = red_rotation(cmap1, ang[ss], total(shift[0,ss]), total(shift[1,ss]), full=ff)
     endif else begin
        cmap11 = red_rotation(cmap1, ang[ss], total(shift[0,ss]), total(shift[1,ss]))
     endelse
        
     ;; Time de-warp (case 1)
     cmap11 = red_stretch(temporary(cmap11), reform(grid[ss,*,*,*]))
     
     ;; Flip any of the axes? (case 1)
     cmap11 = rotate(temporary(cmap11), rot_dir)
     
     print, string(13B)+ inam+'processing scan '+string(ss, format='(I0)')+' of '+string(nscan-1, format='(I0)'), format='(A,$)'
     
     ;; write to disk (case 1)
     if(~keyword_set(float)) then writeu, lun1, fix(round(temporary(cmap11)*1000.)) else  writeu, lun1, float(temporary(cmap11)*1000.)
     
     ;;
     ;; CASE 2
     ;;
     if(dowav) then begin
        wb = (red_mozaic(momfbd_read(wbf[ss])))[x0:x1,y0:y1]
        for ww = 0L, nw - 1 do begin
           
           iwbf = strjoin([self.camwbtag, (strsplit(file_basename(st.ofiles[ww,ss]),'.',/extract))[1:*]],'.')
           iwbf = file_dirname(st.ofiles[ww,ss]) + '/'+iwbf
           
           ;; Load images
           iwb = (red_mozaic(momfbd_read(iwbf)))[x0:x1,y0:y1]
           im = momfbd_read(st.ofiles[ww,ss])
           
           ;; get dewarp from WB
           igrid = red_dsgridnest(wb, iwb, itiles, iclip)
           
           ;; Convolve CMAP and apply wavelength dep. de-warp
           cmap2 =  red_stretch((red_mozaic(red_conv_cmap(cmap, im)))[x0:x1, y0:y1], igrid)
           
           ;; Derotate and shift
           if(full) then begin
              cmap2 = red_rotation(temporary(cmap2), ang[ss], total(shift[0,ss]), $
                                   total(shift[1,ss]), full=ff)
           endif else begin
              cmap2 = red_rotation(temporary(cmap2), ang[ss], total(shift[0,ss]), $
                                   total(shift[1,ss]))
           endelse
           
           ;; Time de-warp
           cmap2 = red_stretch(temporary(cmap2), reform(grid[ss,*,*,*]))
           
           ;; Flip any of the axes?
           cmap2 = rotate(temporary(cmap2), rot_dir)
           
           ;; write to disk
;           writeu, lun2, fix(round((cmap2)*1000.))
           if(~keyword_set(float)) then writeu, lun2, fix(round(temporary(cmap2)*1000.)) else  writeu, lun2, float((temporary(cmap2)*1000.))
        endfor
     endif
  endfor
  print, ' '
  ;; 
  ;; Close files
  ;;
  free_lun, lun1
  if(dowav) then free_lun, lun2
  
  print, inam + 'Output file created -> '+ file1
  if(dowav) then print, inam + 'Output file created -> '+ file2

end
