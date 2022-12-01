; docformat = 'rst'

;+
; Make time-dependent gainfiles.
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
;   mosaic_tag : in, optional, type=string
; 
;      Filename tag that identifies different tiles in data collected
;      in automatic mosaic mode.
;
;   overwrite : in, optional, type=boolean
;
;      Overwrite existing gainfiles.
;
;   tdirs : in, optional, type=strarr
;   
;      Process datasets with tdirs timestamps
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
;   2021-07-06 : OA. Add tdirs keyword.
; 
;   2019-03-28 : MGL. New keyword overwrite. Default is now not to
;                overwrite existing gainfiles.
; 
;   2022-11-15 : MGL. New keyword mosaic_tag.
;
;-
pro red::make_intdif_gains, all = all $
                            , bad = bad $
                            , cam = cam $
                            , debug = debug $
                            , max = max $
                            , min = min $
                            , mosaic_tag = mosaic_tag $
                            , overwrite = overwrite $
                            , pref = pref $
                            , preserve = preserve $
                            , psfw = psfw $
                            , scan = scan $
                            , smallscale = smallscale $
                            , smooth = smooth $
                            , sumlc = sumlc $
                            , tdirs = tdirs $
                            , timeaver = timeaver 

  inam = red_subprogram(/low, calling = inam1)

  if n_elements(timeaver) eq 0 then timeaver = 1L
  if n_elements(min) eq 0 then min = 0.1
  if n_elements(max) eq 0 then max = 4.0
  if n_elements(smooth) eq 0 then smooth = 3.0
  if n_elements(bad) eq 0 then bad = 1.0
  if n_elements(smallscale) EQ 0 THEN smallscale = 1

  red_make_prpara, prpara, bad
  red_make_prpara, prpara, max
  red_make_prpara, prpara, min
  red_make_prpara, prpara, preserve
  ;;  red_make_prpara, prpara, n_elements(psfw); what to do with this?
  red_make_prpara, prpara, smallscale
  red_make_prpara, prpara, smooth
  red_make_prpara, prpara, sumlc
  red_make_prpara, prpara, timeaver

  if ~keyword_set(pref) then pref = '*'
  ;; Search and select output from sum_data_intdif
  if keyword_set(tdirs) then begin
    ww = file_test(self.out_dir+'/cmap_intdif/'+tdirs+'/*'+pref+'*')
    indx = where(ww eq 1)
    nindx = where(ww eq -1)
    if array_equal(nindx, -1,/not_equal) then $
       for ii=0, n_elements(nindx)-1 do $
          print, inam + ' : Subdirectory cmap_intdif/'+tdirs[nindx[ii]]+ ' is empty. Please run sum_data_intdif!'
    if array_equal(indx, -1) then return
    dirs = self.out_dir+'/cmap_intdif/'+tdirs[indx]
  endif else begin
    dirs = file_search(self.out_dir+'/cmap_intdif/*', /test_dir, count = count)
    if count eq 0 then begin
      print, inam + ' : Subdirectory cmap_intdif/ is empty. Please run sum_data_intdif!'
      return
    endif
    if count gt 1 and ~keyword_set(all) then begin
      tmp = red_select_subset(dirs, qstring = 'Select folder(s)' $
                              , indx = sindx, count = Nselect)
      if Nselect eq 0 then return
      dirs = dirs[sindx]
      print, inam + ' : Using -> '+dirs
    endif
  endelse

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
      if keyword_set(mosaic_tag) then begin
        search = dir + '/' + detectors[icam] + '_' + mosaic_tag + '.' + pref + '.intdif.icube'
      endif else begin
        search = dir + '/' + detectors[icam] + '.' + pref + '.intdif.icube'
      endelse
      
      cfile = file_search(search, count = count)
      
      if count eq 0 then begin
        if ~keyword_set(mosaic_tag) then begin
          tmpfiles = file_search(dir + '/' + detectors[icam] + '_mos??.' + pref + '.intdif.icube', count = tmpcount)
          if tmpcount gt 0 then begin
            ;; So this is a directory where data were collected in
            ;; automatic mosaic mode. Find out how many mosaic tiles
            ;; there are and call ourself individually for each tile.
            pos = strpos(file_basename(tmpfiles[-1]), '_mos')
            Nmos = long(strmid(file_basename(tmpfiles[-1]), pos+4, 2))+1
            for imos = 0, Nmos-1 do begin
              mosaic_tag = 'mos'+string(imos, format = '(i02)')
              self ->  make_intdif_gains, all = all $
                                          , bad = bad $
                                          , cam = cams[icam] $
                                          , debug = debug $
                                          , max = max $
                                          , min = min $
                                          , mosaic_tag = mosaic_tag $
                                          , overwrite = overwrite $
                                          , pref = pref $
                                          , preserve = preserve $
                                          , psfw = psfw $
                                          , scan = scan $
                                          , smallscale = smallscale $
                                          , smooth = smooth $
                                          , sumlc = sumlc $
                                          , tdirs = file_basename(dir) $
                                          , timeaver = timeaver
            endfor              ; imos
            continue
          endif
        endif
        print, inam + ' : No files match search string "'+search+'"'
      endif
      
      if count eq 0 then continue
      
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
      ;; Nscans, pref, Nx, Ny, udwav, usettings, ufullstate

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
      udwavd = udwav-mean(udwav)

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

      testflats = 0
      if testflats then begin
        ;; Load cavity-error free flats
        flats = fltarr(Nx, Ny, Ntunings)
        fhdrs2 = list()
        for iwav=0, Ntunings-1 do begin

          self -> get_calib, { camera:cams[icam] $
                               , detector:detectors[icam] $
                               , scannumber:uscan[0] $
                               , framenumber: 0 $
                               , prefilter: pref $
                               , fpi_state: uwav[iwav] $
                               , tuning: uwav[iwav]  $
                               , cam_settings: usettings[iwav] $                               
                               ;;, fullstate: ufullstate[iwav] $                               
                             } $
                             , cflatdata = cflatdata $
                             , cflatname = cflatname 

          print, inam + ' : loading -> '+file_basename(cflatname)

          flats[*,*,iwav] = cflatdata
          fhdrs2.add, red_readhead(cflatname)

        endfor                  ; iwav
      endif
      
      ;; Load cavity-error free gains
      gains = fltarr(Nx, Ny, Ntunings)
      fhdrs = list()
      for iwav=0, Ntunings-1 do begin

        self -> get_calib, { camera:cams[icam] $
                             , detector:detectors[icam] $
                             , scannumber:uscan[0] $
                             , framenumber: 0 $
                             , prefilter: pref $
                             , fpi_state: uwav[iwav] $
                             , tuning: uwav[iwav]  $
                             , cam_settings: usettings[iwav] $ 
                             ;;, fullstate: ufullstate[iwav] $
                           } $
                           , cgaindata = cgaindata $
                           , cgainname = cgainname 

        print, inam + ' : loading -> '+file_basename(cgainname)

        gains[*,*,iwav] = cgaindata
        fhdrs.add, red_readhead(cgainname)

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

        ;; Do we want this scan?
        if n_elements(scan) gt 0 then begin
          if uscan[ss] ne scan then begin
            print, inam + ' : Not wanted -> ' + detectors[icam] + ', Scan '+strtrim(uscan[ss], 2)
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
            if ii eq x0 then cub = double(dat[ii]) else cub += dat[ii]
          endfor                ; ii
          print, ' '
        endif else begin
          if (ox0 ne x0) or (ox1 ne x1) then print, inam + ' : Correcting loaded cube:' else $
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

        for ilc = 0, Nlc-1 do begin
          
          if keyword_set(sumlc) and (ilc gt 0) then begin

            print, inam + ' : re-using spectra for '+ulc[ilc]

          endif else begin

            if ~keyword_set(sumlc) then begin
              ;; Check if scan is already done      
              ofiles = strarr(Ntunings)
              for iwav = 0, Ntunings-1 do begin
                fullstate = strjoin([(strsplit(ufullstate[iwav],'_',/extract))[0:-2],'lc'+strtrim(long(ulc[ilc]), 2)],'_')
                ofiles[iwav] = self -> filenames('scangain' $
                                                 , {camera:cams[icam] $
                                                    , detector:detectors[icam] $
                                                    , cam_settings: usettings[iwav] $ 
                                                    , fullstate: fullstate $
;                                         , fullstate:strjoin([pref $
;                                         , uwav[iwav] $
;                                         , 'lc'+strtrim(long(ulc[ilc]), 2)], '_') $
                                                    , scannumber:uscan[ss] $
                                                   } $
                                                 , mosaic_tag = mosaic_tag $
                                                 , timestamp = imdir $
                                                 , /wild_framenumber $
                                                 , /wild_prefilter $
                                                 , /wild_tuning $
                                                )
                
                
              endfor            ; iwav
              
              done = file_test(ofiles)
              if min(done) eq 1 && ~keyword_set(overwrite) then begin
                print, inam + ' : Already done -> ' + detectors[icam] + ', Scan '+strtrim(uscan[ss], 2)+', LC'+strtrim(ilc, 2) 
                continue
              endif
            end
            
            if ~keyword_set(sumlc) then begin
              cub2 = transpose(cub[*,*,*,ilc], [2,0,1])
            endif else begin
              cub2 = transpose(total(cub,4,/double), [2,0,1])
            endelse

            if keyword_set(psfw) then begin
              ;; Convolve data
              print, inam + ' : convolving data ... ', format='(A,$)'
              for iwav = 0, nw-1 do begin
                cub2[iwav,*,*] = red_convolve(reform(cub2[iwav,*,*]), psf)
              endfor            ; iwav
              print, 'done'
            endif

            ;; Shift spectra, one pixel at a time
            cub1 = cub2
            print, inam+' : shifting cube ... ', format='(A,$)'
            for iy = 0, Ny-1 do begin
              for ix = 0, Nx-1 do begin
                cub1[*,ix,iy] = rdx_cbezier3(udwavd, cub1[*,ix,iy], udwavd+cmap[ix,iy])
              endfor            ; ix
            endfor              ; iy
            print, 'done'

          endelse
          
          ;; Compute new gains
          for iwav = 0, Ntunings-1 do begin
            fullstate = strjoin([(strsplit(ufullstate[iwav],'_',/extract))[0:-2],'lc'+strtrim(long(ulc[ilc]), 2)],'_')
            ofile = self -> filenames('scangain' $
                                      , {camera:cams[icam] $
                                         , detector:detectors[icam] $
                                         , cam_settings: usettings[iwav] $ 
                                         , fullstate: fullstate $
;                                         , fullstate:strjoin([pref $
;                                         , uwav[iwav] $
;                                         , 'lc'+strtrim(long(ulc[ilc]), 2)], '_') $
                                         , scannumber:uscan[ss] $
                                        } $
                                      , mosaic_tag = mosaic_tag $
                                      , timestamp = imdir $
                                      , /wild_framenumber $
                                      , /wild_prefilter $
                                      , /wild_tuning $
                                     )
            
            
            if keyword_set(sumlc) and ilc gt 0 then begin

              print, 'creating link '+ofile
              file_delete, ofile, /allow_nonexistent

              ofile_0 = self -> filenames('scangain' $
                                          , {camera:cams[icam] $
                                             , detector:detectors[icam] $
                                             , cam_settings: usettings[iwav] $ 
                                             , fullstate: fullstate $
;                                             , fullstate:strjoin([pref $
;                                                                  , uwav[iwav] $
;                                                                  , 'lc'+strtrim(long(ulc[0]), 2)], '_') $
                                             , scannumber:uscan[ss] $
                                            } $
                                          , mosaic_tag = mosaic_tag $
                                          , timestamp = imdir $
                                          , /wild_framenumber $
                                          , /wild_prefilter $
                                          , /wild_tuning $
                                         )
              file_link, ofile_0, ofile

            endif else begin 
              
              g = gains[*, *, iwav] * reform(cub1[iwav, *, *]/cub2[iwav, *, *])

              if min(finite(g)) eq 0 then begin
                ;; We need to take care of zeros in cub2, they can cause
                ;; Infinity values in g.
                zindx = where( ~finite(g), Nz)
                ;; Zero pixels that should have been zero anyway
                tmp_img = gains[*, *, iwav] * reform(cub1[iwav, *, *])
                for iz = 0, Nz-1 do $
                   if (tmp_img)[zindx[iz]] eq 0 then g[zindx[iz]] = 0

                ;; Take care of the rest by filling pixels
                zindx = where( ~finite(g), Nz)
                if Nz gt 0 then begin
                  mask = ~finite(g)
                  g[zindx] = 0
                  g = rdx_fillpix(g, mask = mask)
                end
              endif
              
;;;;;; What if we didn't make new gains from flats here but instead
;;;;;; just modified (multiply with cub1/cub2)) the cavityfree gains
;;;;;; we already have on disk? We don't have states to use with
;;;;;; get_calib but we could get states from the flat names. Or use a
;;;;;; makeshift state like when we generate ofile above.
              

              if testflats then begin
                rat = flats[*, *, iwav] * reform(cub2[iwav, *, *]/cub1[iwav, *, *])
                
                g2 = float(self -> flat2gain(temporary(rat), min = min, max = max, bad = bad $
                                             , smooth = smooth $
                                             , preserve = keyword_set(preserve) $
                                             or pref eq '8542' or pref eq '7772'))

                cgplot, g, g2, psym = 3

                stop
              endif
              
              
              ;; Save gains
              print, 'saving '+ ofile
              overwrite = 1     ; make keyword?-----------------------------------------------------------

              ;; Make header
              hdr = fhdrs[iwav]
              red_fitsaddkeyword, hdr, 'FILENAME', ofile
              red_fitsaddkeyword, hdr, 'SCANNUM', uscan[ss]
              red_fitsaddkeyword, hdr, 'LC', long(ulc[ilc]) ; Or something like this ----------------------------

              self -> headerinfo_addstep, hdr, prstep = 'RECIPROCAL' $
                                          , prproc = inam, prpara = prpara
              if min(finite(g)) eq 0 then stop       
              red_writedata, ofile, g, header = hdr,$
                             filetype='fits', overwrite = overwrite
              
            endelse

          endfor                ; iwav
          
        endfor                  ; ilc

      endfor                    ; ss

      undefine, cub
      undefine, mosaic_tag
      free_lun, lun

    endfor                      ; icam

    undefine, cub1 
    undefine, rat 
    undefine, flats

  endfor                        ; idir

end
