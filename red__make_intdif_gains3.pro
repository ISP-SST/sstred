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
;   all : in, optional, type=boolean
;   
;      Process all data sets.
;   
;   
;   
; 
; 
; :history:
;   2013-07-01: JdlCR : Created!
;
;   2013-12-17  PS make smallscale default
;                  if SUMLC is given, don't recompute the gain and use links
;
;   2015-05-07: MGL. With a *, run all directories.
; 
;   2016-02-17 : MGL. New keyword "all".
;
;
;-
pro red::make_intdif_gains3, timeaver = timeaver, sumlc = sumlc, pref = pref, debug = debug $
                             , cam = cam, min=min, max=max, bad=bad, smooth=smooth, psfw = psfw $
                             , preserve = preserve, scan = scan, smallscale = smallscale $
                             , all = all

  inam = 'red::make_intdif_gains3 : '
  if(n_elements(timeaver) eq 0) then timeaver = 1L
  if(n_elements(min) eq 0) then min = 0.1
  if(n_elements(max) eq 0) then max = 4.0
  if(n_elements(smooth) eq 0) then smooth = 3.0
  if(n_elements(bad) eq 0) then bad = 1.0
  if(n_elements(smallscale) EQ 0) THEN smallscale = 1

  ;;
  ;; Search directories
  ;;
  dirs = file_search(self.out_dir+'/cmap_intdif/*', /test_dir, count = count)
  if(count eq 0) then begin
     print, inam + 'You should run red::sum_data_intdif first! -> returning'
     return
  endif
  if(count gt 1) and ~keyword_set(all) then begin
     for ii = 0L, count -1 do print, red_stri(ii)+' -> '+dirs[ii]
     idx = ''
     read, idx, prom = inam+'Select folder (* for all of them): '
     if idx ne '*' then begin
        dirs = dirs[long(idx)]
        
        print, inam + 'Using -> '+dirs
     endif
  ENDIF

  for idir = 0, n_elements(dirs)-1 do begin

     dir = dirs[idir]

     imdir = file_basename(dir)
     print, inam + 'Using folder -> '+imdir
     outdir = self.out_dir + '/gaintables/'+imdir+'/'
     file_mkdir, outdir

     self -> getdetectors, dir = self.data_dir
     cams = [self.camttag, self.camrtag]
     ncam = 2
     
     if(keyword_set(psfw)) then begin
        psf = red_get_psf(2*psfw-1, 2*psfw-1, double(psfw), double(psfw))
        psf /= total(psf,/double)
     endif

     for cc = 0, ncam-1 do begin
        if(n_elements(cam) gt 0) then if(cams[cc] ne cam) then continue
        ;;
        ;; search files
        ;;
        IF (keyword_set(pref)) THEN $
           search = dir + '/' + cams[cc]+'.'+pref+'.intdif.icube' $
        ELSE $
           search = dir + '/' + cams[cc]+'.*.intdif.icube'
        cfile = file_search(search, count = count)
        dfile = dir+'/'+file_basename(cfile,'icube')+'save'
        
        idx = intarr(count)
        k = 0L
        for ii = 0L, count-1 do begin
           if(~file_test(dfile[ii])) then continue
           print, red_stri(k)+' -> '+dfile[k]
           idx[k] = ii
           k+=1
        endfor
        k -= 1
        toread = 0
        if(k gt 0) then read, toread, promp=inam+'select file id : '
        if(k lt 0) then continue
        cfile = cfile[toread]
        dfile = dfile[toread]
        print, inam + 'using -> '+cfile

        ;;
        ;; Open files
        ;;
        restore, dfile
        ;; variables in there:  
        ;;  done, uwav, ulc, uscan, nw, nlc, ns, pref, nx, ny, udwav
        openr, lun, cfile, /get_lun
        dat = assoc(lun, intarr(nx, ny, nw, nlc,/nozero), 512)
        fff = self.out_dir +'/flats/spectral_flats/'+cams[cc]+'.'+pref+'.fit_results.sav'
        if(~file_test(fff)) then begin
           print, inam + 'ERROR, could not find/load -> ', fff
           return
        endif
        print, inam + 'loading -> '+file_basename(fff)
        restore, fff
        cmap = reform((temporary(fit)).pars[1,*,*])
        udwav = double(udwav)

        ;;
        ;; The real cavity-map has quite large shifts, but only the local
        ;; fine structure affects momfbd. With this option, we only
        ;; compensate for the fine local structure, removing the large
        ;; scale features
        ;;
        if(keyword_set(smallscale)) then begin
           print, inam+'correcting only for local line shifts ... ', format='(A,$)'
           npix = 30
           cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
           cpsf /= total(cpsf, /double)
           lscale = red_convolve(cmap, cpsf)
           cmap -= lscale
           print, 'done'
        endif


        ;;
        ;; load cavity-error free flats
        ;;
        flats = fltarr(nx, ny, nw)
        for ww=0, nw-1 do begin
           ffile = self.out_dir + '/flats/'+strjoin([cams[cc],pref,uwav[ww]],'.')+'.unpol.flat'
           if(~file_test(ffile)) then begin
              print, inam + 'ERROR, cannot find flat file '+ffile
              return
           endif
           print, inam + 'loading -> '+file_basename(ffile)
           flats[*,*,ww] = f0(ffile)
        endfor


        ;;
        ;; Existing scans?
        ;;
        idx = where(done ne 0, count)
        if(count eq 0) then begin
           print, inam + 'have you run red::sum_data_intdif? -> returning'
           return
        endif
        t0 = idx[0]
        t1 = idx[n_elements(idx)-1]
        print, inam + 't0 = '+red_stri(t0)+', t1 = '+red_stri(t1)
        print, inam + 'Looping through scans'
        print, ' '
        ;;
        ;; sum images
        ;;
        for ss = t0, t1 do begin 
           if(n_elements(scan) gt 0) then begin
              if(uscan[ss] ne scan) then begin
                 print, inam + 'skipping scan -> '+uscan[ss]
                 continue
              endif
           endif

           ;;
           ;; Get timeaver bounds
           ;;
           dt = timeaver/2
           x0 = (ss-dt)>t0
           x1 = (x0 + timeaver-1)
           if(x1 gt t1) then begin
              x0 = (t1-timeaver+1)>t0
              x1 = t1
           endif
           print, inam + 'x0 = '+red_stri(x0)+', x1 = '+red_stri(x1)

           ;;
           ;; read data
           ;;
           if((ss eq t0) or (n_elements(cub) eq 0)) then begin
              for ii = x0, x1 do begin
                 print, string(13B),inam +'adding t = '+red_stri(ii)+' / '+red_stri(x1), format='(A,A,I0,A,I0,$)' 
                 if(ii eq x0) then cub = float(dat[ii]) else cub += dat[ii]
              endfor
              print, ' '
           endif else begin
              if(ox0 ne x0 OR ox1 NE x1) then print, inam + 'correcting loaded cube:' else $
                 print, inam + 'time-average window has not changed -> not loading data for this scan'
              if(ox0 ne x0) then begin
                 print, '   -> x0(='+red_stri(x0)+') != ox0(='+red_stri(ox0)+') -> removing t='+red_stri( ox0)
                 cub -= dat[ox0]
              endif
              if(ox1 ne x1) then begin
                 print, '   -> x1(='+red_stri(x1)+') != ox1(='+red_stri(ox1)+') -> adding t='+red_stri( x1)
                 cub += dat[x1]
              endif
           endelse
           ox0 = x0
           ox1 = x1

           cub1 = dblarr(nw, nx, ny)
           for ll = 0, nlc-1 do begin
              if(keyword_set(sumlc) and ll gt 0) then begin
                 print, inam + 're-using spectra for '+ulc[ll]
              endif else begin
                 if(~keyword_set(sumlc)) then begin
                    cub2 = float(transpose(cub[*,*,*,ll], [2,0,1])) 
                 endif else begin
                    cub2 = float(transpose(total(cub,4,/double), [2,0,1]))
                 endelse
                 cub1 = double(cub2)

                 ;;
                 ;; Convolve data
                 ;;
                 if(keyword_set(psfw)) then begin
                    print, inam + 'convolving data ... ', format='(A,$)'
                    for ww = 0, nw-1 do begin
                       cub1[ww,*,*] = red_convolve(reform(cub1[ww,*,*]), psf) 
                       cub2[ww,*,*] = red_convolve(reform(cub2[ww,*,*]), psf)
                    endfor
                    print, 'done'
                 endif

                 ;;
                 ;; Shift spectra
                 ;;
                 print, inam+'shifting cube ... ', format='(A,$)'
                 for yy = 0, ny-1 do for xx=0, nx-1 do begin
                    cub1[*,xx,yy] = cbezier3(udwav, cub1[*,xx,yy], udwav+cmap[xx,yy])
                 endfor
                 print, 'done'
              endelse
              
              

              ;;
              ;; Compute new gains
              ;;
              FOR ww = 0, nw-1 DO BEGIN
                 ofile = strjoin([cams[cc],uscan[ss],pref,uwav[ww], ulc[ll]],'.')+'.gain'
                 
                 IF keyword_set(sumlc) AND ll GT 0 THEN BEGIN
                    print, 'creating link '+outdir+ofile
                    file_delete, outdir+ofile, /ALLOW_NONEXISTENT
                    ofile_0 = strjoin([cams[cc],uscan[ss],pref,uwav[ww], ulc[0]],'.')+'.gain'
                    file_link, outdir+ofile_0, outdir+ofile
                 ENDIF ELSE BEGIN 
                    rat = flats[*, *, ww] * reform(cub2[ww, *, *]/cub1[ww, *, *])

                    g = float(self->flat2gain(temporary(rat), min = min, max = max, bad = bad, $
                                            smooth = smooth, preserve = preserve))
                    ;;
                    ;; Save gains
                    ;;
                    print, 'saving '+ outdir+ofile
                    fzwrite, float(g), outdir+ofile, ' '
                 ENDELSE
              endfor            ; ww
              
           endfor               ; ll
           

        endfor                  ; ss
        dum = size(temporary(cub))
        free_lun, lun
     endfor                     ; cc
     cub1 = 0B
     rat = 0B
     flats = 0B
  endfor                        ; idir

end
