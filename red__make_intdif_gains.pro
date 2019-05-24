; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz Rodriguez, ISP
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
;    sumlc  : 
;   
;   
;   
;    psfw  : 
;   
;   
;   
;   
;   all : in, optional, type=boolean
;   
;      Process all data sets.
;   
;   
;   
; 
; 
; :History:
;
;   2013-07-01 : JdlCR. Created!
;
;   2013-12-17 : PS. Make smallscale default. If SUMLC is given, don't
;                recompute the gain and use links
;
;   2015-05-07 : MGL. With a *, run all directories.
; 
;   2016-02-17 : MGL. New keyword "all".
; 
;   2017-04-26 : MGL. Use rdx_cbezier3.
; 
;   2018-09-19 : MGL. Adapt to new codebase.
; 
;   2019-03-28 : MGL. Make gains with /preserve for prefilters with
;                backscatter correction, i.e., 8542 and 7772.
;
;-
pro red::make_intdif_gains, all = all $
                            , bad = bad $
                            , cam = cam $
                            , debug = debug $
                            , max = max $
                            , min = min $
                            , pref = pref $
                            , preserve = preserve $
                            , psfw = psfw $
                            , scan = scan $
                            , smallscale = smallscale $
                            , smooth = smooth $
                            , sumlc = sumlc $
                            , timeaver = timeaver 

  inam = red_subprogram(/low, calling = inam1)

  red_make_prpara, prpara, bad
  red_make_prpara, prpara, max
  red_make_prpara, prpara, min
  red_make_prpara, prpara, preserve
  ;;  red_make_prpara, prpara, n_elements(psfw); what to do with this?
  red_make_prpara, prpara, smallscale
  red_make_prpara, prpara, smooth
  red_make_prpara, prpara, sumlc
  red_make_prpara, prpara, timeaver

  
  if n_elements(timeaver) eq 0 then timeaver = 1L
  if n_elements(min) eq 0 then min = 0.1
  if n_elements(max) eq 0 then max = 4.0
  if n_elements(smooth) eq 0 then smooth = 3.0
  if n_elements(bad) eq 0 then bad = 1.0
  if n_elements(smallscale) EQ 0 THEN smallscale = 1

  ;; Search and select output from sum_data_intdif
  dirs = file_search(self.out_dir+'/cmap_intdif/*', /test_dir, count = count)
  if count eq 0 then begin
    print, inam + ' : Subdirectory cmap_intdif/ is empty. Please run sum_data_intdif!'
    return
  endif
  if count gt 1 and ~keyword_set(all) then begin
    tmp = red_select_subset(dirs, qstring = 'Select folder(s)' $
                            , indx = sindx, count = Nselect)
    if Nselect eq 0 then return
;    for ii = 0L, count -1 do print, red_stri(ii)+' -> '+dirs[ii]
;    idx = ''
;    read, idx, prom = inam+'Select folder (* for all of them): '
;    if idx ne '*' then begin
    dirs = dirs[sindx]
    print, inam + ' : Using -> '+dirs
  endif

  ;; Find the narrowband cameras.
  self -> getdetectors
  cams = *self.cameras
;  detectors = *self.detectors
  is_wb = strmatch(cams,'*-[DW]')
  cams = cams[where(~is_wb, Ncams)]
  detectors = (*self.detectors)[where(~is_wb, Ncams)]
  if Ncams eq 0 then stop

  for idir = 0, n_elements(dirs)-1 do begin

    dir = dirs[idir]

    imdir = file_basename(dir)
    print, inam + ' : Using folder -> '+imdir
    outdir = self.out_dir + '/gaintables/'+imdir+'/'
    file_mkdir, outdir

    if keyword_set(psfw) then begin
      psf = red_get_psf(2*psfw-1, 2*psfw-1, double(psfw), double(psfw))
      psf /= total(psf,/double)
    endif

    for icam = 0, Ncams-1 do begin

      if n_elements(cam) gt 0 && cams[icam] ne cam then continue

      ;; Search files
      if n_elements(pref) gt 0 then $
         search = dir + '/' + detectors[icam]+'.'+pref+'.intdif.icube' $
      else $
         search = dir + '/' + detectors[icam]+'.*.intdif.icube'

      cfile = file_search(search, count = count)
      dfile = dir+'/'+file_basename(cfile,'icube')+'save'
      
      idx = intarr(count)
      k = 0L
      for ii = 0L, count-1 do begin
        if ~file_test(dfile[ii]) then continue
        print, red_stri(k)+' -> '+dfile[k]
        idx[k] = ii
        k += 1
      endfor
      k -= 1
      toread = 0
      if k gt 0 then read, toread, promp=inam+' : select file id : '
      if k lt 0 then continue
      cfile = cfile[toread]
      dfile = dfile[toread]
      print, inam + 'using -> '+cfile

      restore, dfile
      ;; Variables in dfile: done, uwav, ulc, uscan, Ntunings, Nlc,
      ;; Nscans, pref, Nx, Ny, udwav

      ;; Open files
      openr, lun, cfile, /get_lun
      dat = assoc(lun, intarr(Nx, Ny, Ntunings, Nlc,/nozero), 512)
      fff = self.out_dir +'/flats/spectral_flats/'+detectors[icam]+'_'+pref+'_fit_results.sav'
      if ~file_test(fff) then begin
        print, inam + 'ERROR, could not find/load -> ', fff
        return
      endif
      print, inam + ' : loading -> '+file_basename(fff)
      restore, fff
      cmap = reform((fit).pars[1,*,*]) ; Cavity map in [Å]
      udwav = double(udwav)*1d10       ; Wavelength coordinates in [Å]

      ;; The real cavity-map has quite large shifts, but only the
      ;; local fine structure affects momfbd. With this option, we
      ;; only compensate for the fine local structure, removing the
      ;; large scale features.

      if keyword_set(smallscale) then begin
        print, inam+' : Correcting only for local line shifts ... ', format='(A,$)'
        npix = 30
        cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
        cpsf /= total(cpsf, /double)
        lscale = red_convolve(cmap, cpsf)
        cmap -= lscale
        print, 'done'
      endif

      ;; Load cavity-error free flats
      flats = fltarr(Nx, Ny, Ntunings)
      fhdrs = list()
      for iwav=0, Ntunings-1 do begin

        ;; Get this through get_calib method ---------------------------------
;        ffile = self.out_dir + '/flats/'+strjoin([cams[icam],pref,uwav[iwav]],'.')+'.unpol.flat'
        ;;ffile = self.out_dir + '/flats/'+strjoin([detectors[icam], pref, uwav[iwav], 'lc0', 'cavityfree'],'_')+'.flat.fits'
        ffile = fit.oname[iwav]

        if ~file_test(ffile) then begin
          print, inam + ' : ERROR, cannot find flat file '+ffile
          return
        endif
        print, inam + ' : loading -> '+file_basename(ffile)
        flats[*,*,iwav] = red_readdata(ffile, header = fhdr)
        fhdrs.add, fhdr
      endfor                    ; iwav

      ;; Existing scans?
      idx = where(done ne 0, count)
      if count eq 0 then begin
        print, inam + ' : Have you run red::sum_data_intdif? -> returning'
        return
      endif
      t0 = idx[0]
      t1 = idx[n_elements(idx)-1]
      print, inam + ' : t0 = '+red_stri(t0)+', t1 = '+red_stri(t1)
      print, inam + ' : Looping through scans'
      print, ' '


      
      ;; Sum images
      for ss = t0, t1 do begin 
        if n_elements(scan) gt 0 then begin
          if uscan[ss] ne scan then begin
            print, inam + ' : skipping scan -> '+uscan[ss]
            continue
          endif
        endif

        ;; Get timeaver bounds
        dt = timeaver/2
        x0 = (ss-dt) > t0
        x1 = (x0 + timeaver-1)
        if(x1 gt t1) then begin
          x0 = (t1-timeaver+1)>t0
          x1 = t1
        endif
        print, inam + ' : x0 = '+red_stri(x0)+', x1 = '+red_stri(x1)

        ;; Read data
        if (ss eq t0) or (n_elements(cub) eq 0) then begin
          for ii = x0, x1 do begin
            print, string(13B),inam +' : adding t = '+red_stri(ii)+' / '+red_stri(x1), format='(A,A,I0,A,I0,$)' 
            if ii eq x0 then cub = float(dat[ii]) else cub += dat[ii]
          endfor                ; ii
          print, ' '
        endif else begin
          if (ox0 ne x0) or (ox1 ne x1) then print, inam + 'correcting loaded cube:' else $
             print, inam + ' : time-average window has not changed -> not loading data for this scan'
          if ox0 ne x0 then begin
            print, '   -> x0(='+red_stri(x0)+') != ox0(='+red_stri(ox0)+') -> removing t='+red_stri( ox0)
            cub -= dat[ox0]
          endif
          if ox1 ne x1 then begin
            print, '   -> x1(='+red_stri(x1)+') != ox1(='+red_stri(ox1)+') -> adding t='+red_stri( x1)
            cub += dat[x1]
          endif
        endelse
        ox0 = x0
        ox1 = x1

        cub1 = dblarr(Ntunings, Nx, Ny)
        for ilc = 0, Nlc-1 do begin
          
          if keyword_set(sumlc) and (ilc gt 0) then begin

            print, inam + ' : re-using spectra for '+ulc[ilc]

          endif else begin

            if ~keyword_set(sumlc) then begin
              cub2 = float(transpose(cub[*,*,*,ilc], [2,0,1])) 
            endif else begin
              cub2 = float(transpose(total(cub,4,/double), [2,0,1]))
            endelse
            cub1 = double(cub2)

            if keyword_set(psfw) then begin
              ;; Convolve data
              print, inam + ' : convolving data ... ', format='(A,$)'
              for iwav = 0, nw-1 do begin
                cub1[iwav,*,*] = red_convolve(reform(cub1[iwav,*,*]), psf) 
                cub2[iwav,*,*] = red_convolve(reform(cub2[iwav,*,*]), psf)
              endfor            ; iwav
              print, 'done'
            endif
            
            ;; Shift spectra
            print, inam+' : shifting cube ... ', format='(A,$)'
            for iy = 0, Ny-1 do for ix=0, Nx-1 do begin
              cub1[*,ix,iy] = rdx_cbezier3(udwav, cub1[*,ix,iy], udwav+cmap[ix,iy])
            endfor              ; ix,iy
            print, 'done'

          endelse
          
          ;; Compute new gains
          for iwav = 0, Ntunings-1 do begin

            ;; Move creation of filename to crisp::filenames? ------------------------------------
            ;; Should include timestamp of data directory!
            
            ofile = strjoin([detectors[icam] $
                             , string(uscan[ss], format = '(I05)') $
                             , pref $
                             , uwav[iwav] $
                             , 'lc'+strtrim(long(ulc[ilc]), 2) $
                            ], '_') + '.gain.fits'

            if keyword_set(sumlc) and ilc gt 0 then begin

              print, 'creating link '+outdir+ofile
              file_delete, outdir+ofile, /allow_nonexistent
              ;; Move creation of filename to crisp::filenames?  ------------------------------------
              ofile_0 = strjoin([detectors[icam] $
                                 , string(uscan[ss], format = '(I05)') $
                                 , pref $
                                 , uwav[iwav] $
                                 , 'lc'+strtrim(long(ulc[0]), 2) $
                                ], '_')+'.gain.fits'
              file_link, outdir+ofile_0, outdir+ofile

            endif else begin 
              
              rat = flats[*, *, iwav] * reform(cub2[iwav, *, *]/cub1[iwav, *, *])

              g = float(self -> flat2gain(temporary(rat), min = min, max = max, bad = bad $
                                          , smooth = smooth $
                                          , preserve = keyword_set(preserve) $
                                          or pref eq '8542' or pref eq '7772'))

              ;; Save gains
              print, 'saving '+ outdir+ofile
              ;; fzwrite, float(g), outdir+ofile, ' '
              overwrite = 1     ; make keyword?-----------------------------------------------------------

              ;; Make header
              hdr = fhdrs[iwav]
              red_fitsaddkeyword, hdr, 'FILENAME', ofile
              red_fitsaddkeyword, hdr, 'SCANNUM', uscan[ss]
              red_fitsaddkeyword, hdr, 'LC', long(ulc[ilc]) ; Or something like this ----------------------------

              self -> headerinfo_addstep, hdr, prstep = 'Gain making' $
                                          , prproc = inam, prpara = prpara
              
              
              red_writedata, outdir+ofile, g, header = hdr,$
                             filetype='fits', overwrite = overwrite
              
            endelse

          endfor                ; iwav
          
        endfor                  ; ilc

      endfor                    ; ss

      undefine, cub
      free_lun, lun

    endfor                      ; icam

    undefine, cub1 
    undefine, rat 
    undefine, flats

  endfor                        ; idir

end
