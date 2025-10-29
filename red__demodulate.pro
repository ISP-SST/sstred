; docformat = 'rst'

;+
; Demodulate restored images for a particular scan and tuning state.
; 
; This method does (or at least is meant to do) the same as
; pol__demodulate.pro in the master (old CRISP) branch, except that it
; reads the inverse modulation matrices from disk rather than having
; pointers to them from within the pol class object. (We don't use a
; pol class in the new code base.) We rely on red_mozaik to handle
; both the inverse mod matrices and the image data the same way when
; it comes to clipping.
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    outname : in, type=string
; 
;        The file name in which to write the demodulated Stokes cube.
;
;    immr : in, type=array
; 
;        The inverse demodulation matrix for the Crisp-R camera.
;
;    immt : in, type=array
;
;        The inverse demodulation matrix for the Crisp-T camera.
;
;
;
; :Keywords:
; 
;    clips : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;    cmap : in, optional, type=array
;
;       The cavity map, to be distorted like the image data and stored
;       in the stokes cube.
;
;    nthreads : in, optional, type=integer
;
;       The number of threads to use.
; 
;    nbrfac : in, optional, type=float
;
;       A scaling factor to apply to NBR data to get the correct units.
;
;    nbrstates : in, type=structarr
;   
;       An array of structs with the states of the input NBR images.
; 
;    nbtfac : in, optional, type=float
;
;       A scaling factor to apply to NBT data to get the correct units.
;
;    nbtstates : in, type=structarr
;   
;       An array of structs with the states of the input NBT images.
;
;    noremove_periodic : in, optional, type=boolean
;
;       Do not attempt to use Fourier filtering to remove polcal
;       periodic artifacts.
;
;    tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
; 
;    wbg : in, optional, type=image
;
;      The global WB image to use as destretching reference. If this
;      keyword is not used, no destretching will be done.
; 
;    wbstates : in, optional, type=structarr
;   
;       An array of structs with the states of the input WB images. If
;       provided, will use them to stretch the NB images when
;       combining them.
; 
; 
; :History:
; 
;   2018-10-09 : MGL. First version, based on Jaime's pol::demodulate.
; 
;   2018-12-05 : MGL. New keyword cmap.
; 
;   2018-12-21 : MGL. New keywords smooth_by_kernel and
;                smooth_by_subfield.
; 
;   2019-12-06 : MGL. New keyword noremove_periodic. Do optional
;                Fourier filtering to remove periodic artifacts from
;                polcal. Read header info from filter file and include
;                filtering in prstep.
;
;   2019-04-09 : OA. Added check for NaN values in momfbd files they
;                are filled with zeroes (it happens out of the limb).
; 
;   2020-06-16 : MGL. In FITS header, write SOLARNET recommended
;                PRSTEP and POLCCONV.
;
;   2020-10-01 : JdlCR. Use the new red_stretch_linear to apply the
;                wbcorrection. Also allow for nearest neighbor
;                interpolation.
; 
;   2020-12-09 : MGL. Remove statistics calculations.
; 
;   2021-10-20 : MGL. Save some time by storing smoothed versions of
;                the inverse demodulation matrices.
;
;   2021-12-10 : JdlCR. Make use of the new libgrid routines, now
;                ported to rdx and maintainable by us.
; 
;   2022-07-29 : MGL. Change from a CRISP:: method to a RED:: method. 
; 
;   2025-01-22 : MGL. Set step PRMODE for new demodulation /
;                flat-fielding scheme.
;
;   2025-02-04 : MGL. Added the telescope calibration model parameters
;                and the telescope Mueller matrix as PRREF to the
;                demodulation step information.
;
;   2025-02-20 : MGL. Adapt to new camera alignment model.
; 
;   2025-05-28 : MGL. Adapt to mix of polarimetric and
;                non-polarimetric data. New keyword notelmat.
; 
;-
pro red::demodulate, outname, immr, immt $
                     , clips = clips $
                     , cmap = cmap $
                     , nbrfac = nbrfac $
                     , nbrstates = nbrstates $
                     , nbtfac = nbtfac $
                     , nbtstates = nbtstates $
                     , nearest = nearest $
;                     , newalign = nal $
                     , noremove_periodic = noremove_periodic $
                     , notelmat = notelmat $
                     , nthreads = nthreads $
                     , overwrite = overwrite $
                     , smooth_by_kernel = smooth_by_kernel $
                     , smooth_by_subfield = smooth_by_subfield $
                     , tiles = tiles $
                     , units = units $
                     , wbg = wbg $
                     , wbstates = wbstates $
                     , wcs = wcs
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)                                  

  if file_test(outname) and ~keyword_set(overwrite) then begin
    red_message, 'Already exists: '+outname
    return
  endif

  ;; Camera/detector identification
  self -> getdetectors
  wbindx      = where(strmatch(*self.cameras,'*-W'))
  wbcamera    = (*self.cameras)[wbindx[0]]
  wbdetector  = (*self.detectors)[wbindx[0]]
  nbtindx     = where(strmatch(*self.cameras,'*-T')) 
  nbtcamera   = (*self.cameras)[nbtindx[0]]
  nbtdetector = (*self.detectors)[nbtindx[0]]
  nbrindx     = where(strmatch(*self.cameras,'*-R')) 
  nbrcamera   = (*self.cameras)[nbrindx[0]]
  nbrdetector = (*self.detectors)[nbrindx[0]]

  instrument = (strsplit(wbcamera, '-', /extract))[0]

  red_make_prpara, prpara, clips         
  red_make_prpara, prpara, nbrfac 
  red_make_prpara, prpara, nbtfac
  red_make_prpara, prpara, smooth_by_kernel
  red_make_prpara, prpara, smooth_by_subfield
  red_make_prpara, prpara, tiles
  red_make_prpara, prpara, notelmat
  
  ;; Defaults of keywords
  if n_elements(Nthreads) eq 0 then Nthreads = 1
  if n_elements(units) eq 0 then units = 'D.N.'
  
  outdir = file_dirname(outname)
  
  ;; Should be four images from each camera (four LC states)
  Nstokes = 4
;  Nlc = 4
;  if n_elements(nbrstates) ne Nlc then stop
;  if n_elements(nbtstates) ne Nlc then stop
;  if n_elements(wbstates)  ne Nlc then stop
  Nlc = n_elements(wbstates)
  if Nlc ne n_elements(nbtstates) || Nlc ne n_elements(nbrstates) then stop

  ;; Nlc should be 4 (polarimetry) or 1 (no polarimetry, e.g., 3999 Å cont)
  if Nlc ne 4 && Nlc ne 1 then stop
  
  
;  Nelements = 16
  
  rfiles = nbrstates.filename
  tfiles = nbtstates.filename
  wfiles = wbstates.filename

  ;; Get the alignment mode, projective or polywarp.
  alignment_model = red_align_model(wfiles[0])
  red_make_prpara, prpara, alignment_model
  
  if rdx_filetype(rfiles[0]) eq 'MOMFBD' then begin
    ;; MOMFBD output
    for ilc = 0, Nlc-1 do begin
      red_append, rimg, momfbd_read(rfiles[ilc])
      red_append, timg, momfbd_read(tfiles[ilc])
;      red_append, wimg, momfbd_read(wfiles[ilc])  ; Never used
    endfor

    ;; FOV is given in the momfbd output, use it to crop the modulation matrices.
;  mr = momfbd_read(tfiles[0], /names) ; Use /names to avoid reading
;  the data parts
    mr = rimg[0]
    roi = mr.roi
    if n_elements(mr.patch) eq 1 then begin
      ;; Single subfield, no mosaic, no cropping by red_readdata.
      margin = 0
    endif else begin
      ;; This is how much red_readdata crops a mosaic
      margin = mr.margin
    endelse
    x0 = roi[0] + margin
    x1 = roi[1] - margin
    y0 = roi[2] + margin
    y1 = roi[3] - margin
  endif else begin
    ;; FITS files from red::bypass_momfbd
    hdr = red_readhead(wfiles[0])
    x0 = fxpar(hdr, 'CROP_X0')
    x1 = fxpar(hdr, 'CROP_X1')
    y0 = fxpar(hdr, 'CROP_Y0')
    y1 = fxpar(hdr, 'CROP_Y1')
    roi = [x0, x1, y0, y1]
    margin = 0
  endelse
  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1

  prefilter = nbstates[0].prefilter
;  camw = wbstates[0].camera
;  camt = nbtstates[0].camera
;  camr = nbrstates[0].camera
  
  ;; Get some header info
  
  for ilc = 0, Nlc-1 do begin
    wbhead = red_readhead(wfiles[ilc])
    red_fitspar_getdates, wbhead $
                          , date_beg = date_beg $
                          , date_avg = date_avg $
                          , date_end = date_end 
    red_append, tbegs, red_time2double((strsplit(date_beg,'T',/extract))[1])
    red_append, tavgs, red_time2double((strsplit(date_avg,'T',/extract))[1])
    red_append, tends, red_time2double((strsplit(date_end,'T',/extract))[1])
    red_append, xps,   fxpar(wbhead, 'XPOSURE')
    red_append, texps, fxpar(wbhead, 'TEXPOSUR')
    red_append, nsums, fxpar(wbhead, 'NSUMEXP')
  endfor                        ; ilc
  tbeg = min(tbegs)
  tavg = mean(tavgs)
  tend = max(tends)
  
  if n_elements(wcs) gt 0 then wcs.time = tavg
  
  ;; Create arrays for the image mozaics
  img_t = fltarr(Nx, Ny, Nlc)
  img_r = fltarr(Nx, Ny, Nlc)

  dims = size(immt, /dim)
  Nxx = dims[1]                 ; Detector size
  Nyy = dims[2]

  if Nlc gt 1 then begin
    if 0 then begin
      
      ;; Adapt this block to new code base!
      
      utflat = fltarr(Nxx, Nyy)
      urflat = fltarr(Nxx, Nyy)
      tgain = fltarr([Nxx,Nyy,Nlc])
      rgain = fltarr([Nxx,Nyy,Nlc])
      
      ;; utflat and urflat are the "unpolarized" ("cavityfree") flats
      ;; for the current state (ignoring the LC state).
      self -> get_calib, nbtstates[0], cflatdata = utflat
      self -> get_calib, nbrstates[0], cflatdata = urflat
      
      for ilc = 0, Nlc-1 do begin
        
        ;; tgain and rgain are the gaintables for the current state,
        ;; including the LC states.
        
        ;; Should be the same gains that were used when momfbding. That
        ;; is, the scan-specific flats. Make sure to switch to that if
        ;; this part of the code is to be activated!
        
        self -> get_calib, nbtstates[ilc], gaindata = gaindata ; Read the gain
        tgain[0, 0, ilc] = gaindata
        
        self -> get_calib, nbrstates[ilc], gaindata = gaindata ; Read the gain
        rgain[0, 0, ilc] = gaindata

      endfor                    ; ilc

      red_message, 'Applying flat ratios to the inverse modulation matrix.'
      immt_dm = red_demat(tgain, utflat, immt)
      immr_dm = red_demat(rgain, urflat, immr)
      print, 'done'
      
    endif else begin
      immt_dm = immt
      immr_dm = immr
    endelse
  end
  
  ;; Mozaic images
  nmask = bytarr(Nx, Ny)+1
  for ilc = 0L, Nlc-1 do begin
    ;;stop
    im = red_readdata(rfiles[ilc]) * nbrfac
    ;;im = red_mozaic(rimg[ilc], /crop) * nbrfac
    nmask = nmask and finite(im)
    img_r[*,*,ilc] = im
    im = red_readdata(tfiles[ilc]) * nbtfac
    ;;im = red_mozaic(timg[ilc], /crop) * nbtfac
    nmask = nmask and finite(im)
    img_t[*,*,ilc] = im
  endfor                        ; ilc
  
  for ilc = 0L, Nlc-1 do begin
    img_r[*,*,ilc] *= nmask
    img_t[*,*,ilc] *= nmask
  endfor                        ; ilc
  
  if Nlc gt 1 then begin

    if keyword_set(smooth_by_subfield) then begin
      
      red_message, ['The smooth_by_subfield code does not support the new polywarp' $
                    , 'camera alignment model and is anyway not regularly used.']

      ;; Divide modulation matrices into momfbd-subfields, correct for
      ;; image projection, smooth with PSF from each subfield, make
      ;; mosaics of the result.
      
      immts = red_matrix2momfbd(timg, immt_dm, amap = amapt)
      immrs = red_matrix2momfbd(rimg, immr_dm, amap = amapr)

      ;; Mozaic inverse modulation matrix
      mymt = fltarr(Nlc, Nstokes, Nx, Ny)
      mymr = fltarr(Nlc, Nstokes, Nx, Ny)
      for ilc = 0L, Nlc-1 do begin
        for istokes = 0L, Nstokes-1 do begin
          mymr[ilc,istokes,*,*] = red_mozaic(immrs[ilc,istokes], /crop) 
          mymt[ilc,istokes,*,*] = red_mozaic(immts[ilc,istokes], /crop)
        endfor                  ; istokes
      endfor                    ; ilc

      ;; if keyword_set(cmap) then cmap = (red_mozaic(red_conv_cmap(cmap,*self.timg[1])))[x0:x1,y0:y1]

    endif else begin
      
      mymrname = outdir + '/mymr.fits'      
      mymtname = outdir + '/mymt.fits'      
      
      if file_test(mymrname) and file_test(mymtname) then begin
        mymr = readfits(mymrname)    
        mymt = readfits(mymtname)    
        ;; Check dimensions?
      endif else begin
        
        mymt = fltarr(Nlc, Nstokes, Nx, Ny)
        mymr = fltarr(Nlc, Nstokes, Nx, Ny)
        

        ;; Do we need to reform the inverse modulation matrices? At least
        ;; we need to do the projection matrix thing.
        
        immr_dm = reform(immr_dm, [Nlc, Nstokes, Nxx, Nyy])
        immt_dm = reform(immt_dm, [Nlc, Nstokes, Nxx, Nyy])
        
        for ilc = 0L, Nlc-1 do begin
          red_progressbar, ilc, Nlc, 'Apply camera alignment to inverse demodulation matrices for LC'+strtrim(ilc, 2)
          for istokes = 0L, Nstokes-1 do begin
            ;; Apply the geometrical mapping and clip to the FOV of the
            ;; momfbd output.
            tmp = red_apply_camera_alignment(reform(immr_dm[ilc, istokes, *, *]) $
                                             , alignment_model, instrument+'-R' $
                                             , pref = prefilter $
                                             , amap = amapr $
                                             , /preserve_size)
            mymr[ilc,istokes,*,*] = tmp[roi[0]+margin:roi[1]-margin $
                                        , roi[2]+margin:roi[3]-margin]
            tmp = red_apply_camera_alignment(reform(immt_dm[ilc,istokes,*,*]) $
                                             , alignment_model, instrument+'-T' $
                                             , pref = prefilter $
                                             , amap = amapt $
                                             , /preserve_size)
            mymt[ilc,istokes,*,*] = tmp[roi[0]+margin:roi[1]-margin $
                                        , roi[2]+margin:roi[3]-margin]
          endfor                ; istokes
        endfor                  ; ilc

        if keyword_set(smooth_by_kernel) then begin

          ;; Smooth the inverse modulation matrices by a Gaussian kernel

          dpix = round(smooth_by_kernel)*3
          if (dpix/2)*2 eq dpix then dpix -= 1
          dpsf = double(smooth_by_kernel)
          psf = red_get_psf(dpix, dpix, dpsf, dpsf)
          psf /= total(psf, /double)

          ;;for ii=0, Nelements-1 do mm[ii,*,*] = red_convolve(reform(mm[ii,*,*]), psf)
          ;; Smooth the inverse modulation matrix
          for ilc = 0L, Nlc-1 do begin
            red_progressbar, ilc, Nlc, 'Smooth demodulation matrices for LC'+strtrim(ilc, 2)
            for istokes = 0L, Nstokes-1 do begin
              ;;mymr[ilc,istokes,*,*] = red_convolve(reform(immr_dm[ilc, istokes, *, *]), psf)
              ;;mymt[ilc,istokes,*,*] = red_convolve(reform(immt_dm[ilc, istokes, *, *]), psf)
              tmp = red_fillnan(reform(mymr[ilc, istokes, *, *]))
              mymr[ilc,istokes,*,*] = red_convolve(tmp, psf)
              tmp = red_fillnan(reform(mymt[ilc, istokes, *, *]))
              mymt[ilc,istokes,*,*] = red_convolve(tmp, psf)
            endfor              ; istokes
          endfor                ; ilc
          
        endif

        writefits, mymrname, mymr    
        writefits, mymtname, mymt 
        
      endelse                   ; read stored mymr and mymt

    endelse                     ; smooth_by_subfield
  endif

  ;; Load the simultaneous WB images.  
  ;;img_wb = fltarr(dim[0], dim[1], Nlc)
  img_wb = fltarr(Nx, Ny, Nlc)
  for ilc = 0L, Nlc-1 do img_wb[*,*,ilc] = red_readdata(wfiles[ilc]) 

  
  if Nlc gt 1 then begin

    rest = fltarr(Nx, Ny, 4)
    resr = fltarr(Nx, Ny, 4)  

    if n_elements(wbg) gt 0 then begin

      ;; Demodulate with destretching
      
      if n_elements(cmap) ne 0 then cmapp = 0.
      for ilc = 0L, Nlc-1 do begin
        grid = rdx_cdsgridnest(wbg, img_wb[*,*,ilc], tiles, clips, nthreads=nthreads)
        
        for istokes = 0L, Nstokes-1 do begin
          rest[*,*,istokes] += rdx_cstretch(reform(mymt[ilc,istokes,*,*]) * img_t[*,*,ilc] $
                                            , grid, nthreads=nthreads)
          resr[*,*,istokes] += rdx_cstretch(reform(mymr[ilc,istokes,*,*]) * img_r[*,*,ilc] $
                                            , grid, nthreads=nthreads)
        endfor                  ; istokes
        if n_elements(cmap) ne 0 then cmapp += rdx_cstretch(cmap, grid, nthreads=nthreads)
      endfor                    ; ilc
      if n_elements(cmap) ne 0 then cmapp /= Nlc
      
    endif else begin

      ;; No global WB image, demodulate without destretching.
      
      for ilc = 0L, Nlc-1 do begin
        for istokes = 0L, 3 do begin
          rest[*,*,istokes] += reform(mymt[ilc,istokes,*,*]) * img_t[*,*,ilc]
          resr[*,*,istokes] += reform(mymr[ilc,istokes,*,*]) * img_r[*,*,ilc]
        endfor                  ; istokes
      endfor                    ; ilc
      if n_elements(cmap) ne 0 then cmapp = cmap
      
    endelse
    print, 7777
    ;; These are now Stokes cubes!
    img_t = temporary(rest)
    img_r = temporary(resr)

    dum = size(img_t,/dim)
    drm = dum / 8.0 
    xx0 = round(drm[0] - 1)
    xx1 = round(dum[0] - drm[0] - 1)
    yy0 = round(drm[1] - 1)
    yy1 = round(dum[1] - drm[1] - 1)
    
    mediant = median(img_t[xx0:xx1,yy0:yy1,0])
    medianr = median(img_r[xx0:xx1,yy0:yy1,0])
    aver = (mediant + medianr) / 2.
    sct = aver / mediant
    scr = aver / medianr
    
    res = (sct * (img_t) + scr * (img_r)) / 2.
    
    red_message, 'Combining data from transmitted and reflected camera'
    print, '   -> Average Intensity = '  + red_stri(aver)
    print, '   -> Tcam scale factor -> ' + red_stri(sct) + ' (after '+red_stri(nbtfac)+')'
    print, '   -> Rcam scale factor -> ' + red_stri(scr) + ' (after '+red_stri(nbrfac)+')'

    
    ;; Telescope model 
    line = (strsplit(wbstates[0].fpi_state,'_',/extract))[0]
    red_message, 'Detected spectral line -> '+line
    
    red_logdata, self.isodate, tavg, azel = azel,  tilt = tilt
    print, 8888
    if min(finite([azel, tilt])) eq 0 then begin
      red_message, 'red_logdata returned non-finite azimuth, elevation, or tilt for t = ' + strtrim(tavg, 2)
      print, '  azel : ', azel
      print, '  tilt : ', tilt
      lfile = 'downloads/sstlogs/positionLog_'+red_strreplace(self.isodate,'-','.',n=2)+'_final'
      print, '  Will try to fix this by re-downloading the turret logfile,'
      print, lfile
      file_delete, lfile, /allow_nonexistent
      red_logdata, self.isodate, tavg, azel = azel,  tilt = tilt
      if min(finite([azel, tilt])) eq 0 then begin
        red_message, 'Re-downloading the turret log did not help. Giving up.'
        retall
      endif
    endif
    
    year = (strsplit(self.isodate, '-', /extract))[0]
    if ~keyword_set(notelmat) then begin
      mtel = red_telmat(line, {TIME:tavg, AZ:azel[0], ELEV:azel[1], TILT:tilt} $
                        , /no_zero, year=year, model_parameters = telcal_parameters)
      imtel = invert(mtel) 
      imtel /= imtel[0]
      print, 9999
      ;; Apply the telescope Muller matrix
      res1 = fltarr(Nx, Ny, Nstokes)
      for j=0, Nstokes-1 do for i=0, Nstokes-1 do $
         res1[*, *, j] += res[*, *, i] * imtel[i, j]
      if keyword_set(testing) then begin
        window, xs = 4*(xx1-xx0+1), ys = yy1-yy0+1, 12
        for iiii = 0, 3 do tvscl, res[xx0:xx1,yy0:yy1,iiii], iiii
        window, xs = 4*(xx1-xx0+1), ys = yy1-yy0+1, 13
        for iiii = 0, 3 do tvscl, res1[xx0:xx1,yy0:yy1,iiii], iiii
        stop
      endif
      res = temporary(res1)
    endif
    
  endif else begin

    ;; Non-polarimetric: I = img_t + img_r, Q = U = V = 0.
    
    res = fltarr(Nx, Ny, NStokes)

    if n_elements(wbg) gt 0 then begin

      ;; Use destretching
      grid = rdx_cdsgridnest(wbg, img_wb[*,*,0], tiles, clips, nthreads=nthreads)
      img_t = rdx_cstretch(img_t[*,*,0], grid, nthreads=nthreads)
      img_r = rdx_cstretch(img_r[*,*,0], grid, nthreads=nthreads)
      
      if n_elements(cmap) ne 0 then cmap = rdx_cstretch(cmap, grid, nthreads=nthreads)

    endif 

    res[*, *, 0] = (img_t + img_r) / 2.
    
  endelse

  if keyword_set(nosave) then return
  
  ;; Save result as a fitscube with all the usual metadata.
  
  file_mkdir, outdir

  dims = [Nx, Ny, 1, Nstokes, 1]
  res = reform(res, dims)

  ;; Make header
  hdr = red_readhead(nbrstates[0].filename)
  check_fits, res, hdr, /update, /silent
  red_fitsaddkeyword, anchor = anchor, hdr, 'FILENAME', file_basename(outname)
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1


  if Nlc eq 1 then begin

    ;; No polarimetry
    
    ;; Add info about this step
    self -> headerinfo_addstep, hdr $
       , prstep = 'DEMODULATION' $
       , prpara = prpara $
       , prmode = 'Non-polarimetric tuning point' $
       , prproc = inam

  endif else begin
    ;; Set prmode to indicate whether polcal was done with the new LC
    ;; flat fielding or not. When no such mode is set, it should be
    ;; understood to mean that the old method was used. If cavity free
    ;; flats for each LC state exist, then the new method was used If
    ;; cavity free flats without LC state exist, then the old method was
    ;; used. If both exist, then we stop and ask or tell the user to
    ;; start over.
    self -> get_calib, nbtstates[0], cgainname = cgainname, clcgainname = clcgainname
    c_exists   = file_test(cgainname[0])
    clc_exists = file_test(clcgainname[0])
    case 1 of

      ~c_exists && ~clc_exists : begin
        red_message, 'Cannot happen! Investigate!'
        stop
      end
      
      c_exists && clc_exists : begin
        red_message, 'We have both LC cavity free gains and non-LC cavity free gains,' $
                     + 'the latter probably left from earlier processing.' $
                     + 'Please delete momfbd output, the gaintables/, flats/spectral_flats/,' $
                     + 'cmap_intdif/, polcal_subs/, polcal/, and prefilter_fits/ subdirectories,' $
                     + 'and then start over.'
        retall
      endcase

      c_exists : begin
        prmode = 'Non-LC gain polcal'
      endcase

      clc_exists : begin
        prmode = 'LC gain polcal'
      endcase

    endcase

    ;; Add info about this step
    if keyword_set(notelmat) then begin
      self -> headerinfo_addstep, hdr $
         , prstep = 'DEMODULATION' $
         , prpara = prpara $
         , prmode = prmode $
         , prproc = inam
    endif else begin
      self -> headerinfo_addstep, hdr $
         , prstep = 'DEMODULATION' $
         , prpara = prpara $
         , prmode = prmode $
         , prproc = inam $
         , prref = ['Telescope calibration model parameters: '+json_serialize(telcal_parameters) $
                    , 'Telescope Mueller matrix: '+json_serialize(mtel)]
    endelse
    
    if ~keyword_set(noremove_periodic) and $
       (file_test('polcal/periodic_filter_'+prefilter+'.fits') or $
        file_test('polcal/periodic_filter_remapped_'+prefilter+'.fits')) then begin
      ;; Fourier-filter images to get rid of periodic fringes      
      if file_test('polcal/periodic_filter_'+prefilter+'.fits') then begin
        ;; CRISP class
        filt = shift(readfits('polcal/periodic_filter_'+prefilter+'.fits', fhdr), Nxx/2, Nyy/2)
        ;; This reversing should actually depend on the align_map of the
        ;; cameras, i.e., the signs of the diagonal elements [0,0] and
        ;; [1,1]. 
        filt[1:*, 1:*] = reverse(filt[1:*, 1:*], 1)
      endif else begin
        ;; RED class
        filt = shift(readfits('polcal/periodic_filter_remapped_'+prefilter+'.fits', fhdr), Nxx/2, Nyy/2)
        ;; This filter is oriented as WB data. So should the data be at
        ;; this point.
      endelse
      for istokes = 0L, Nstokes-1 do begin
        red_missing, res[*,*,0, istokes, 0], nmissing = Nmissing $
                     , indx_missing = indx_missing, indx_data = indx_data
        if Nmissing gt 0 then begin
          bg = (res[*,*,0, istokes, 0])[indx_missing[0]]
          ;;  if istokes eq 0 then bg = (res[*,*,0, istokes, 0])[indx_missing[0]] else bg = 0
        endif else begin
          bg = median(res[*,*,0, istokes, 0])
        endelse
        frame = red_centerpic(res[*,*,0, istokes, 0], sz = Nxx, z = bg)
        indx_nan = where(~finite(frame), Nnan)
        if Nnan gt 0 then frame[indx_nan] = bg
        filtframe = float(fft(fft(frame)*filt, /inv))
        bg = median(filtframe)
        if Nnan gt 0 then filtframe[indx_nan] = bg
        filtframe = red_centerpic(filtframe, xs = Nx, ys = Ny)
        if Nmissing gt 0 then filtframe[indx_missing] = bg
        ;;   if istokes eq 2 then stop ; <------------------------------------------------------------- ***
        res[*,*,0, istokes, 0] = filtframe
      endfor                    ; ilc
      ;; Add info about this step
      taper_width = fxpar(fhdr, 'HOLEWDTH', count = Nkey)
      if Nkey gt 0 then red_make_prpara, filt_prpara, taper_width        
      hole_x = fxpar(fhdr, 'HOLE_X', count = Nkey)
      if Nkey gt 0 then red_make_prpara, filt_prpara, hole_x
      hole_y = fxpar(fhdr, 'HOLE_Y', count = Nkey)
      if Nkey gt 0 then red_make_prpara, filt_prpara, hole_y
      self -> headerinfo_addstep, hdr $
         , prstep = 'FIXED-PATTERN-REMOVAL' $
         , prpara = filt_prpara $
         , prproc = inam
    endif
  
  endelse
  
  ;; New date, at better position
  red_fitsaddkeyword, anchor = anchor, after = 'SOLARNET', /force, hdr $
                      , 'DATE', red_timestamp(/iso), 'Creation UTC date of FITS header '

  ;; Remove some irrelevant keywords
  red_fitsdelkeyword, hdr, 'TAB_HDUS' ; No tabulated headers in summed file
  red_fitsdelkeyword, hdr, 'FRAMENUM' ; No particular frame number
  red_fitsdelkeyword, hdr, 'FRAMEINC' ; Makes no sense to keep
  red_fitsdelkeyword, hdr, 'CADENCE'  ; Makes no sense to keep
  red_fitsdelkeyword, hdr, 'DADCMODE' ; Too technical?
  red_fitsdelkeyword, hdr, 'DADCTHR'  ; Too technical?
  red_fitsdelkeyword, hdr, 'DADCGR'   ; Too technical?
  
  ;; Set various keywords in hdr. DATE* are useful.
  ;; Also make sure we implement the Stokes dimension properly!
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'POLCCONV', '(+HPLT,-HPLN,+HPRZ)', '1st axis toward Solar N, 2nd E, 3rd observer'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DATE-BEG', self.isodate + 'T' + red_timestring(tbeg) $
                      , 'Start time of combined observation'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DATE-AVG', self.isodate + 'T' + red_timestring(tavg) $
                      , 'Average time of combined observation'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DATE-END', self.isodate + 'T' + red_timestring(tend) $
                      , 'End time of combined observation'
  red_fitsaddkeyword, anchor = anchor, hdr, 'XPOSURE', total(xps)
  red_fitsaddkeyword, anchor = anchor, hdr, 'TEXPOSUR', mean(texps)
  red_fitsaddkeyword, anchor = anchor, hdr, 'NSUMEXP', round(total(nsums))

  ;; Get the frame numbers involved
  for ilc = 0, Nlc-1 do begin
    nhdr = red_readhead(nbrstates[ilc].filename)
    red_append, fnumsum, red_expandrange(fxpar(nhdr, 'FNUMSUM'))
  endfor                        ; ilc
  fnumsum = fnumsum[sort(fnumsum)]
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'FNUMSUM', red_collapserange(fnumsum,ld='',rd='') $
                      , 'List of frame numbers used'
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'CAMERA', nbrcamera+','+nbtcamera
  red_fitsaddkeyword, anchor = anchor, hdr $
                      , 'DETECTOR', nbrdetector+','+nbtdetector

  nbrhead = red_readhead(rfiles[0])
  nbthead = red_readhead(tfiles[0])
  red_fitscopykeyword, anchor = anchor, hdr, 'DETGAIN',  nbrhead, nbthead
  red_fitscopykeyword, anchor = anchor, hdr, 'DETOFFS',  nbrhead, nbthead
  red_fitscopykeyword, anchor = anchor, hdr, 'DETMODEL', nbrhead, nbthead
  red_fitscopykeyword, anchor = anchor, hdr, 'DETFIRM',  nbrhead, nbthead
  red_fitscopykeyword, anchor = anchor, hdr, 'DETTEMP',  nbrhead, nbthead ; Should average over all frames?
;  red_fitsaddkeyword, hdr, 'DETGAIN',  strtrim(fxpar(nbrhead,'DETGAIN'), 2) + ',' + strtrim(fxpar(nbthead,'DETGAIN'), 2)
;  red_fitsaddkeyword, hdr, 'DETOFFS',  strtrim(fxpar(nbrhead,'DETOFFS'), 2)  + ',' + strtrim(fxpar(nbthead,'DETOFFS'), 2)
;  red_fitsaddkeyword, hdr, 'DETMODEL', strtrim(fxpar(nbrhead,'DETMODEL'), 2) + ',' + strtrim(fxpar(nbthead,'DETMODEL'), 2)
;  red_fitsaddkeyword, hdr, 'DETFIRM',  strtrim(fxpar(nbrhead,'DETFIRM'), 2)  + ',' + strtrim(fxpar(nbthead,'DETFIRM'), 2)
;  red_fitsaddkeyword, hdr, 'DETTEMP',  strtrim(fxpar(nbrhead,'DETTEMP'), 2)  + ',' + strtrim(fxpar(nbthead,'DETTEMP'), 2)

  ;; Add "global" metadata
  red_metadata_restore, anchor = anchor, hdr
  
  self -> fitscube_initialize, outname, hdr, lun, fileassoc, dims 
  for istokes = 0, Nstokes-1 do begin
    red_fitscube_addframe, fileassoc, res[*, *, 0, istokes, 0] $
                           , istokes = istokes    
  endfor                        ; istokes
  

  ;; Close the file
  self -> fitscube_finish, lun, wcs = wcs

  
  if n_elements(cmapp) ne 0 then red_fitscube_addcmap, outname, reform(cmapp, Nx, Ny, 1, 1, 1)
  
end
