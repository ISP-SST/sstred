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
;    state  : 
;   
;   
;   
;    tiles  : 
;   
;   
;   
;    clip  : 
;   
;   
;   
;    no_destretch  : 
;   
;   
;   
;    no_filter  : 
;   
;   
;   
;    overwrite  : 
;   
;   
;   
;    img  : 
;   
;   
;   
;    sp  : 
;   
;   
;   
;    power_sp  : 
;   
;   
;   
;    cmap  : 
;   
;   
;   
;    nosave : 
;   
;   
;   
;    noclip  : 
;   
;   
;   
;    savecams  : 
;   
;   
;   
;    mdate  : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-12 : MGL. Use red_read_azel, not read_azel.
; 
; 
;-
pro pol::demodulate, state = state, tiles = tiles, clip = clip, no_destretch = no_destretch, $
                      no_filter = no_filter, overwrite = overwrite, img = res, sp = sp, power_sp = power_sp, $
                      cmap = cmap, nosave=nosave, noclip = noclip, savecams = savecams, mdate = mdate, $
                      noflat = noflat, ext_time = ext_time
   
  inam = 'pol::demodulate : '
   
  ;; shall the state be demodulated?
   
  if(keyword_set(state)) then begin
     bla = strjoin((strsplit(self.state, '.',/extract))[1:*], '.')
     if(state NE bla) then begin
        print, inam + self.state + ' != ' + state + ' -> Skipping demodulation!'
        return
     endif     
  endif

  if(~keyword_set(tiles) OR (~keyword_set(clip))) then begin
     tiles = [8,16,32,64]
     clip = [8,4,4,2]
  endif

  ;; Check output file
   
  outdir = file_dirname(self.tfiles[0]) + '/stokes/'
  outname = 'stokesIQUV.'+self.state+'.f0'
   
  if(~keyword_set(overwrite)) then begin
     if file_test(outdir + outname) then begin
        print, inam + 'WARNING, output file exists -> skipping'
        print, inam + ' -> use keyword overwrite = 1 to overwrite the file!'
        return
     endif
  endif
   
  destretch = 1B
  if(keyword_set(no_destretch) OR (self.destretch eq 0)) then destretch = 0B
   
  ;; load images into the object
   
  self -> loadimages
  
  ;; load files
  if(~keyword_set(noflat)) then begin
     utflat = f0(self.utflat)
     urflat = f0(self.urflat)
     tgain = fltarr([size(utflat,/dim), 4])
     rgain = fltarr([size(urflat,/dim), 4])
     for ii = 0,3 do begin
        print, inam + 'loading -> '+self.ftfiles[ii]
        tgain[*,*,ii] = f0(self.ftfiles[ii])
        print, inam + 'loading -> '+self.frfiles[ii]
        rgain[*,*,ii] = f0(self.frfiles[ii])
     endfor
  endif

  ;; Apply flat ratio to IMMT (remember! -> IMMT = Inverse Modulation Matrix!)
   
  immt = *(*self.immt)
  immr = *(*self.immr)
  
  if(~keyword_set(noflat)) then begin
     print, inam + 'applying flat ratios to the inverse modulation matrix (s) ... ', format='(A,$)'
     immt = red_demat(tgain, utflat, temporary(immt))
     immr = red_demat(rgain, urflat, temporary(immr))
     print, 'done'
  endif

  ;; fugly fix for the case when the WB channel is mirrored
  wbclip = self.wbclip
  if wbclip[0] gt wbclip[1] then begin
    print, inam + 'Wideband channel is mirrored: Reversing x'
    immt = reverse( immt, 2, /OVERWRITE )
    immr = reverse( immr, 2, /OVERWRITE )
  endif
  if wbclip[2] gt wbclip[3] then begin
    print, inam + 'Wideband channel is mirrored: Reversing y'
    immt = reverse( immt, 3, /OVERWRITE )
    immr = reverse( immr, 3, /OVERWRITE )
  endif

  ;; Convert the demodulation matrix into patches (like the momfbd
  ;; images) and convolve with the PSFs from the patches.
   
  immt = red_matrix2momfbd((*self.timg[0]), (*self.timg[1]),$
                           (*self.timg[2]),(*self.timg[3]), immt)
  immr = red_matrix2momfbd((*self.rimg[0]), (*self.rimg[1]),$
                           (*self.rimg[2]),(*self.rimg[3]), immr)

  ;; Create arrays for the mozaics

  ;;tmp = red_mozaic(*self.timg[0])
  ;;dim = red_getborder(tmp, x0, x1, y0, y1)
  ;;dim = size(tmp, /dim)
  x0 = self.x0
  x1 = self.x1
  y0 = self.y0
  y1 = self.y1
  if(keyword_set(noclip)) then begin
     tmp = size(red_mozaic((*self.timg[1])), /dim)
     x0 = 0L
     x1 = tmp[0] - 1
     y0 = 0L
     y1 = tmp[1] - 1
  endif
  
  dim = [x1 - x0 + 1, y1 - y0 + 1]
                                ;
  img_t = fltarr(dim[0], dim[1], 4)
  img_r = fltarr(dim[0], dim[1], 4)
  mymt = fltarr(4,4,dim[0], dim[1])
  mymr = fltarr(4,4,dim[0], dim[1])

  ;; Printout

  print, inam + 'Image dimensions after clipping:'
  print, ' -> Nx = ' + red_stri(dim[0])
  print, ' -> Ny = ' + red_stri(dim[1])

  ;; Mozaic images and demodulation matrix
   
  for ii = 0L, 3 do begin
     img_t[*,*,ii] = (red_mozaic(*self.timg[ii]))[x0:x1, y0:y1]
     img_r[*,*,ii] = (red_mozaic(*self.rimg[ii]))[x0:x1, y0:y1]
     for jj = 0L, 3 do mymt[jj,ii,*,*] = (red_mozaic(immt[jj,ii]))[x0:x1, y0:y1]
     for jj = 0L, 3 do mymr[jj,ii,*,*] = (red_mozaic(immr[jj,ii]))[x0:x1, y0:y1]
  endfor
  if(keyword_set(cmap)) then cmap = (red_mozaic(red_conv_cmap(cmap,*self.timg[1])))[x0:x1,y0:y1]
   
  ;; power spectrum (r-camera)
   
  if(keyword_set(sp)) then begin
     tmp = img_r[*,*,3]
     dim = size(tmp, /dim)
     me = median(tmp)
     ;;w = red_taper(dim[0] * [1., 0., 1. / 24.])
     w = red_taper2(dim[0], dim[1], 1./24.)
     power_sp = fft((tmp - me) * w, 1)
     power_sp = shift(conj(power_sp) * power_sp, -dim/2)
     power_sp/= max(power_sp)
     print, inam + 'returning power_spectrum (output keyword)'
     return
  endif
      
  ;; filter ?
   
  if(~keyword_set(no_filter)) then begin
     filterf = file_dirname(self.tfiles[0])+'/filter.'+red_stri(dim[0])+'.f0'
     if(file_test(filterf)) then begin
        print, inam+'Using filter image: '+filterf
        filter = f0(filterf)
        filter = smooth(float(temporary(filter)), 3, /edge_truncate)
        filter = complex(filter, filter)
                                ;
        print, inam+'Filtering ... ', format='(A,$)'
        for ii = 0L, 3 do img_t[*,*,ii] = red_fftfilt(reform(img_t[*,*,ii]), filter)
        for ii = 0L, 3 do img_r[*,*,ii] = red_fftfilt(reform(img_r[*,*,ii]), filter)
        print, 'done'
   
     endif
  endif

  ; Filter states together?
  
  


  ;; Destretch ?
   
  if(destretch) then begin
   
     ;; Load WB images
   
     img_wb = fltarr(dim[0], dim[1], 4)
     for jj = 0L, 3 do img_wb[*,*,jj] = (red_mozaic(momfbd_read(self.wbfiles[jj])))[x0:x1, y0:y1]
     wb = (red_mozaic(momfbd_read(self.wb)))[x0:x1, y0:y1]
      
     ;; measure offsets and apply
      
     grid0 = red_dsgridnest(wb, img_wb[*,*,0], tiles, clip)
     grid1 = red_dsgridnest(wb, img_wb[*,*,1], tiles, clip)
     grid2 = red_dsgridnest(wb, img_wb[*,*,2], tiles, clip)
     grid3 = red_dsgridnest(wb, img_wb[*,*,3], tiles, clip)
      
     if(keyword_set(cmap)) then cmap = red_stretch(temporary(cmap), grid1)

     rest = fltarr(dim[0], dim[1], 4)
     resr = fltarr(dim[0], dim[1], 4)
      
     for ii = 0L, 3 do begin
        rest[*,*,ii] = red_stretch(reform(mymt[0,ii,*,*]) * img_t[*,*,0], grid0)  + $
                       red_stretch(reform(mymt[1,ii,*,*]) * img_t[*,*,1], grid1)  + $
                       red_stretch(reform(mymt[2,ii,*,*]) * img_t[*,*,2], grid2)  + $
                       red_stretch(reform(mymt[3,ii,*,*]) * img_t[*,*,3], grid3) 
 
        resr[*,*,ii] = red_stretch(reform(mymr[0,ii,*,*]) * img_r[*,*,0], grid0)  + $
                       red_stretch(reform(mymr[1,ii,*,*]) * img_r[*,*,1], grid1)  + $
                       red_stretch(reform(mymr[2,ii,*,*]) * img_r[*,*,2], grid2)  + $
                       red_stretch(reform(mymr[3,ii,*,*]) * img_r[*,*,3], grid3) 
     endfor
     img_t = temporary(rest)
     img_r = temporary(resr)
  endif else begin
     print, inam + 'Demodulating data without de-stretching'
                                
     rest = fltarr(dim[0], dim[1], 4)
     resr = fltarr(dim[0], dim[1], 4)
                                
     for ii = 0L, 3 do begin
        rest[*,*,ii] = reform(mymt[0,ii,*,*]) * img_t[*,*,0] + $
                       reform(mymt[1,ii,*,*]) * img_t[*,*,1] + $
                       reform(mymt[2,ii,*,*]) * img_t[*,*,2] + $
                       reform(mymt[3,ii,*,*]) * img_t[*,*,3]
                                
        resr[*,*,ii] = reform(mymr[0,ii,*,*]) * img_r[*,*,0] + $
                       reform(mymr[1,ii,*,*]) * img_r[*,*,1] + $
                       reform(mymr[2,ii,*,*]) * img_r[*,*,2] + $
                       reform(mymr[3,ii,*,*]) * img_r[*,*,3]
     endfor
                                
     img_t = temporary(rest)
     img_r = temporary(resr)
  endelse

  ;;  if(keyword_set(experiment)) then begin
  ;;     for ii = 0, 3 do img_r[*,*,ii] = deconvolve_straylight(img_r[*,*,ii], 0.059,fwhm = 7.59)
  ;; endif 
                                
  ;; Combine cameras scaling to the median
   
  dum = size(img_t,/dim)
  drm = dum / 8.0 
  xx0 = drm[0] - 1
  xx1 = dum[0] - drm[0] - 1
  yy0 = drm[1] - 1
  yy1 = dum[1] - drm[1] - 1
   
  aver = 0.5 * (mean(img_t[xx0:xx1,yy0:yy1,0]) + mean(img_r[xx0:xx1,yy0:yy1,0]))
  sct = aver / mean(img_t[xx0:xx1,yy0:yy1,0])
  scr = aver / mean(img_r[xx0:xx1,yy0:yy1,0])
   
  res = (sct * (img_t) + scr * (img_r)) * 0.5
   
  print, inam + 'Combining data from transmitted and reflected camera'
  print, ' -> Average Intensity = ' + red_stri(aver)
  print, ' -> Tcam scale factor -> '+red_stri(sct)
  print, ' -> Rcam scale factor -> '+red_stri(scr)
   
  ;; Average obs. time
   
  time_obs = 0.d0
  for ii = 0L, 3 do time_obs+= red_time2double((*self.timg[ii]).time)
  time_obs = red_time2double(time_obs * 0.25d0, /dir)
  
  if(n_elements(ext_time) gt 0) then begin
     iscan = long(self.scan)
     if(n_elements(ext_time) gt iscan + 1) then time_obs = ext_time[iscan]
  endif
   
  ;; telescope model
   
  line = (strsplit(self.state,'.',/extract))[1]
  print, inam+'Detected spectral line -> '+line
  if(file_test(self.telog)) then begin

     ;; mtel =
     ;; sst_mueller_all((*self.timg[0]).date, time_obs, self.telog, line)
     if(n_elements(mdate) eq 0) then mdate = strjoin(strsplit((*self.timg[0]).date,'-/.',/extra),'/')
     year = long((strsplit(mdate,'/',/extract))[0])
     telpos = red_read_azel( self.telog,mdate)
     print, inam + 'time_obs = '+time_obs
     if(line eq '6300') then begin
        mtel = red_telmat('6302', telpos, time_obs, /no_zero, year=year)
     endif else mtel = red_telmat(line, telpos, time_obs, /no_zero, year=year)
     
     imtel = invert(mtel) 
     imtel /= imtel[0]

     ;; Apply the matrix

     res1 = fltarr(dim[0], dim[1], 4)
     FOR j=0, 3 DO FOR i=0, 3 DO res1[*, *, j] += res[*, *, i] * imtel[i, j]
     res = temporary(res1)
     
  endif else print, inam + 'WARNING, SST position LOG not found -> telescope polarization not corrected!!!!'
      
  ;; Save result
                                
  if(keyword_set(savecams)) then begin
     odir = file_dirname(self.tfiles[0]) + '/demodulated_cameras/'
     file_mkdir, odir
     
     tcam = (strsplit(file_basename(self.tfiles[0]), '.',/extract))[0]
     rcam = (strsplit(file_basename(self.rfiles[0]), '.',/extract))[0]

     ofil = tcam + '.'+self.state + '.f0'
     res1 = fltarr(dim[0], dim[1],4)
     FOR j=0, 3 DO FOR i=0, 3 DO res1[*, *, j] += img_t[*, *, i] * imtel[i, j]
     fzwrite,  temporary(res1) * sct, odir + ofil,' '
     
     ofil = rcam + '.'+self.state + '.f0'
     res1 = fltarr(dim[0], dim[1],4)
     FOR j=0, 3 DO FOR i=0, 3 DO res1[*, *, j] += img_r[*, *, i] * imtel[i, j]
     fzwrite,  temporary(res1) * scr, odir + ofil,' '

     return
  endif

  if(~keyword_set(nosave)) then begin
     file_mkdir, outdir
     print, inam + 'saving file -> '+ outdir + outname
     head = 'TIME_OBS='+time_obs+' DATE_OBS='+(*self.timg[0]).date
     fzwrite, res, outdir + outname, head+''
                                
     if(keyword_set(cmap)) then begin
        odir = file_dirname(self.tfiles[0]) + '/cavity_map/'
        file_mkdir, odir
        ofile = 'cmap.' + self.state + '.f0'
        print, inam + 'saving cavity map -> '+odir+ofile
        fzwrite, float(cmap), odir+ofile, head + ''
     endif
  endif

  ;; Do not remove! otherwise the IDL session will eat a lot of memory!
  ;; At least the images need to be unloaded

  self -> unloadimages          ; Deallocate pointers

  ;; Is this needed?

  immt = 0B
  immr = 0B
  img_t = 0B
  img_r = 0B
  img_wb = 0B
  mymt = 0B
  mymr = 0B
  wb = 0B
  tmp = 0B
  tiles = 0B
  clip = 0B
  ;;res = 0B
   
  return
end
