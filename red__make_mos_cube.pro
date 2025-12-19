; docformat = 'rst'

;+
; Make a fitscube based on a mosaic observation.
;
; This is enteded mostly for showing, not for science. So we are not
; super concerned about getting the different tunings and Stokes
; components perfectly aligned.
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
; :Params:
; 
;    dir : in, type="string or strarr"
; 
;      This could be the path to a cfg directory with results_mos??
;      directories in it, an array of cfg/results_mos?? directory
;      paths, or a regular expression matching the cfg/results_mos??
;      direcories. This method will then make scan cubes for the
;      output in the results_mos?? directories and then stitch those
;      cubes together as a large mosaic scan cube. The dir could also
;      be a path to the directory with the already made scan cubes for
;      the individual mosaic tiles.
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2025-08-20 : MGL. First version. Based in part on
;                 red__quicklook_mosaic.pro.
;
;    2025-10-09 : MGL. Call make_scan_cube with
;                 /no_intensitycorr_timecheck. Use r0 for tile
;                 weighting.
; 
;-
pro red::make_mos_cube, dir $
                        , individual_destretch = individual_destretch $
                        , overwrite = overwrite $
                        , new_scan_cubes = new_scan_cubes $
                        , nthreads = nthreads $
                        , no_cmap = no_cmap $
                        , no_destretch = no_destretch $
                        , _ref_extra = extra
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if keyword_set(new_scan_cubes) || keyword_set(redemodulate) then overwrite = !true
  
  ;; Make prpara
  red_make_prpara, prpara, dir
  red_make_prpara, prpara, no_cmap
  red_make_prpara, prpara, no_destretch

  
  ;; Expand dir, in case it is a regular expression.
  edir = file_search(dir, count = Ndirs)

  if Ndirs eq 0 then begin
    red_message, 'No such directory: ' + dir
    return
  endif 

  dirsplt = strsplit(edir[0], '/', /extract)
  
  
  
  if Ndirs eq 1 then begin

    ;; This could be a cfgdir or a dir with scan cubes in it.
    
    if dirsplt[-1] eq 'cfg' then begin
      ;; If a cfgdir, the results directories should be in it
      cdir = edir               ; The cfg directory
      rdirs = file_search(cdir+'/results_mos??', count = Nrdirs)
      sdir = strjoin(['cubes_mos/', dirsplt[-3:-2]], '/')+'/' ; The scan cubes dir
    endif

  endif else begin

    rdirs = edir                ; The results directories
    Nrdirs = n_elements(rdirs)
    
    sdir = strjoin(['cubes_mos/', dirsplt[-4:-3]], '/') + '/' ; The scan cubes dir
    
  endelse

  if n_elements(sdir) eq 0 then begin

    ;; If sdir is undefined at this point, the input dir should be the
    ;; directory with the scancubes in it.

    sdir = edir                 ; The scan cubes dir
    
  endif else begin

    ;; We need to make the scan cubes, or at least check that they are
    ;; made.
    
    Nfiles = lonarr(Nrdirs)
;    ;; First check if there is momfbd output in all directories 
;    for idir = 0, Nrdirs-1 do begin
;      ;; Should allow for fits output also
;      files = file_search(rdirs[idir]+'/*.momfbd', count = N)
;      Nfiles[idir] = N
;    endfor                      ; idir
    
    red_message, 'Will make scan cubes from data in these directories:'
    print, '   ' + rdirs
    file_mkdir, sdir

    for idir = 0, Nrdirs-1 do begin

      files = file_search(rdirs[idir]+'/*.momfbd', count = N)
      Nfiles[idir] = N
      
      if Nfiles[idir] gt 0 then begin

        tilenum = long(strmid(rdirs[idir],1, /rev))
        
        self -> red::make_scan_cube, rdirs[idir] $
           , filename = sname $
           , mosaic = tilenum $
           , /no_intensitycorr_timecheck $
           , overwrite = new_scan_cubes $
           , odir = sdir $
           , scanno = '0' $
           , _strict_extra = extra
        
        if ~keyword_set(no_cmap) then begin
          undefine, outname
          self -> fitscube_cmapcorr, sname $
             , outname = outname $
             , overwrite = new_scan_cubes
          red_append, sfiles, outname
        endif else begin
          red_append, sfiles, sname
        endelse

      endif else begin

        red_message, 'No momfbd output in '+rdirs[idir]
        
      endelse
      
    endfor                      ; idir 

    if min(Nfiles) eq 0 then begin
      indx = where(Nfiles eq 0)
      red_message, 'There was momfbd output in these directories:'
      print, rdirs[indx]
      red_message, 'Run momfbd on them (or wait for momfbd to finish) and then run make_mos_cube again.'
      return
    endif

    
  endelse 
  
;  ;; At this point, the scan cubes should exist in sdir.
;  if ~keyword_set(no_cmap) then begin
;    ;; Do we want cubes corrected for cavity maps?
;    sfiles = file_search(sdir+'/*_mos[0-9][0-9]*corrected*cmapcorr*fits', count = Ntiles)
;  endif else begin
;    ;; Or do we want whithout cmap correction?
;    sfiles = file_search(sdir+'/*_mos[0-9][0-9]*corrected.fits', count = Ntiles)
;  endelse

  sfiles = sfiles[sort(sfiles)] ; Just to make sure...

  Ntiles = n_elements(sfiles)

  fnsplt = strsplit(file_basename(sfiles[0]), '_', /extract)
  filename = strjoin([fnsplt[0:2], 'mos', fnsplt[4:*]], '_')

  if keyword_set(individual_destretch) then filename = red_strreplace(filename, '_mos_', '_mos_indstr_')

  filename = sdir + filename

  if file_test(filename) && ~keyword_set(overwrite) then begin

    red_message, 'Output mosaic file exists already: ' + filename
    red_message, 'Call with /overwrite to remaking it.'
    
  endif
  
  ;; The plan is now to use the coordinates in the scancube files and
  ;; stitch together mosaics, much like we do in quicklook_mosaic, but
  ;; for all tunings and Stokes parameters, and put the mosaic in a
  ;; scancube file of its own.

  ;; We should use fov_mask and possibly allow the user to define
  ;; another (smaller) mask since crisp2 is bad along the perimeter.


  ;; FOV dimensions could maybe vary slightly between fitscubes? Use
  ;; the max.
  Nx = max(red_fitsgetkeyword_multifile(sfiles, 'NAXIS1'))
  Ny = max(red_fitsgetkeyword_multifile(sfiles, 'NAXIS2'))

  r0 = red_fitsgetkeyword_multifile(sfiles, 'ATMOS_R0') 

  ;; Get date_avg from file headers
  date_avg = red_fitsgetkeyword_multifile(sfiles, 'DATE-AVG', count = Ndates)
  time_avg = red_time2double(strmid(date_avg, 11))

  pref = red_fitsgetkeyword(sfiles[0], 'FILTER1')
  image_scale = self->imagescale(pref)

  Ntuning = red_fitsgetkeyword(sfiles[0], 'NAXIS3')
  Nstokes = red_fitsgetkeyword(sfiles[0], 'NAXIS4')
  
  ;; Start making the FITS header
  hdr = red_sumheaders(sfiles)

  ;; This is probably unnecessary, will get set when we initialize the
  ;; output fitscube.
  red_fitsaddkeyword, hdr, 'NAXIS', 5
  red_fitsaddkeyword, hdr, 'NAXIS1', 0         ; Set later
  red_fitsaddkeyword, hdr, 'NAXIS2', 0         ; Set later
  red_fitsaddkeyword, hdr, 'NAXIS3', Ntuning
  red_fitsaddkeyword, hdr, 'NAXIS4', Nstokes
  red_fitsaddkeyword, hdr, 'NAXIS5', 1
  
  
  ;; Read WB tiles
  wbims = fltarr(Nx, Ny, Ntiles)
  mask = bytarr(Nx, Ny) + 1B    ; We could have individual masks but use a common one for now.

  t = dblarr(Ntiles)
  wave = fltarr(Ntuning)
  hpln = fltarr(2, 2, Ntiles)
  hplt = fltarr(2, 2, Ntiles)

  diskpos = fltarr(2, Ntiles) 
  
  xpos = fltarr(Ntiles)
  ypos = fltarr(Ntiles)

  ;; WCS info
  wcs_alltiles = replicate({ wave:dblarr(2,2) $
                             , hplt:dblarr(2,2) $
                             , hpln:dblarr(2,2) $
                             , time:dblarr(2,2) $
                           }, Ntuning, Ntiles)

  
  for itile = 0, Ntiles-1 do begin

    fxread, sfiles[itile], wbim, extension = 'WBIMAGE'

    ;; The mask, pad with NaN
    red_missing, wbim, /inplace, missing_type_wanted = 'NaN'
    wbim = red_centerpic(wbim, xs = Nx, ys = Ny, z = !Values.F_NaN)
    mask and= finite(wbim)

    ;; The image, pad with median
    red_missing, wbim, /inplace, missing_type_wanted = 'median', missing_value = bg
    wbims[*, *, itile] = red_centerpic(wbim, xs = Nx, ys = Ny, z = bg)

    ;; WCS coordinates
    red_fitscube_getwcs, sfiles[itile], coordinates=coordinates

    wcs_alltiles[*,itile] = coordinates
    
    t[itile] = mean(coordinates.time)
    hpln[*, *, itile] = coordinates[0].hpln
    hplt[*, *, itile] = coordinates[0].hplt

    diskpos[0, itile] = mean(coordinates[0].hpln)
    diskpos[1, itile] = mean(coordinates[0].hplt)
    
    xpos[itile] = mean(coordinates[0].hpln) / image_scale
    ypos[itile] = mean(coordinates[0].hplt) / image_scale

    if itile eq 0 then begin
      ;; Same for all tiles
      wave = coordinates.wave[0, 0]
    endif 
    
  endfor                        ; itile

  ;; When calculating the mosaicking weights, we use the mask but we
  ;; want to ignore any zeros within the main area.
  mmask = mask*0 + 1
  if mask[0, 0] eq 0 then begin
    indx = search2d(mask gt 0, 0,0,0,0)
    if n_elements(indx) gt 0 then mmask[indx] = 0
  endif
  if mask[0,-1] eq 0 then begin
    indx = search2d(mask gt 0, 0,Ny-1,0,0)
    if n_elements(indx) gt 0 then mmask[indx] = 0
  endif
  if mask[-1,0] eq 0 then begin
    indx = search2d(mask gt 0, Nx-1,0,0,0)
    if n_elements(indx) gt 0 then mmask[indx] = 0
  endif
  if mask[-1,-1] eq 0 then begin
    indx = search2d(mask gt 0, Nx-1,Ny-1,0,0)
    if n_elements(indx) gt 0 then mmask[indx] = 0
  endif
  ww = morph_distance(mmask,neigh=3)
  ww /= max(ww)
  weight = sin(ww*!pi/2.)^2

  
  xpos = xpos - min(xpos) + Nx/2 + 50
  ypos = ypos - min(ypos) + Ny/2 + 50

  ;; Size of array
  Sx = round(max(xpos) + Nx/2 + 50)
  Sy = round(max(ypos) + Ny/2 + 50)

  if Sx gt Ntiles*Nx || Sy gt Ntiles*Ny then begin
    red_message, 'Failed to get good tile positions'
    stop
  endif

  fac = 1
  imap = fltarr(Sx/fac, Sy/fac, Ntiles) ; image map
  wmap = fltarr(Sx/fac, Sy/fac, Ntiles) ; weight map
  mmap = fltarr(Sx/fac, Sy/fac, Ntiles) ; mask map

  llx = lonarr(Ntiles)
  urx = lonarr(Ntiles)
  lly = lonarr(Ntiles)
  ury = lonarr(Ntiles)


  gridtiles = [8, 16, 32, 64, 84]
  gridclips = [12, 8,  4,  2 , 1]

  for itile = 0, Ntiles-1 do begin

    red_progressbar, itile, Ntiles, 'Mosaic WB tile '+red_stri(itile)

    llx[itile] = round(xpos[itile]/fac - Nx/2/fac)
    lly[itile] = round(ypos[itile]/fac - Ny/2/fac)
    urx[itile] = llx[itile] + round(Nx/fac - 1)
    ury[itile] = lly[itile] + round(Ny/fac - 1)


    
;    llx2 = 2*(xpos[itile]/fac - Nx/2/fac)
;    lly2 = 2*(ypos[itile]/fac - Ny/2/fac)
;    urx2 = llx2 + Nx/fac - 1
;    ury2 = lly2 + Ny/fac - 1


    ;; Scale the weights with average r0 squared for each tile cube,
    ;; in order to give more weight to better tiles.
    if fac eq 1 then begin
      imap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] = mask * wbims[*, *, itile]
      wmap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] = mask * weight * r0[itile]^2
      mmap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] = mask
    endif else begin
      imap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] $
         = congrid( mask * wbims[*, *, itile], Nx/fac, Ny/fac, cubic = -0.5 )
      wmap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] $
         = congrid( mask * weight * r0[itile]^2, Nx/fac, Ny/fac, cubic = -0.5 ) >0
      mmap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] $
         = congrid( mask, Nx/fac, Ny/fac, cubic = -0.5 ) $
         * (wmap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] gt 0)
    endelse
    
;        red_show, tiles, /re
    

  endfor                        ; itile
  
  qq = Ntiles <9
  colors = red_distinct_colors(qq, /num)
  while n_elements(colors) lt Ntiles do colors = [colors, colors]
  colors = colors[0:Ntiles-1]
  
;  if Ntiles gt 9 then colors=[colors,colors[0:Ntiles-10]]

  xcp = median(diskpos[0, *]) - Sx*image_scale/2.
  ycp = median(diskpos[1, *]) - Sy*image_scale/2.

  debug = 1
  if keyword_set(debug) then begin
 
    cgwindow
    cgplot, /add, /nodata, [0], [0] $
            , xrange = [0, float(Sx)] $
            , yrange = [0, float(Sy)] $
            ;;, xrange = [0, Sx*image_scale], yrange = [0, Sy*image_scale] $
            , aspect = Sy/float(Sx) $
            , xtitle = 'X / 1 pix', ytitle = 'Y / 1 pix' $
            , title = 'Tile positions '+fxpar(hdr,'DATE-OBS')+' '+pref
    cgplot, /add, /over, xpos, ypos, psym=16, color='red'

    cgwindow
    cgplot, /add, /nodata, [0], [0] $
            , xrange = [min(diskpos[0, *]), max(diskpos[0, *])] + image_scale*float(Nx)/2*[-1, 1] $
            , yrange = [min(diskpos[1, *]), max(diskpos[1, *])] + image_scale*float(Ny)/2*[-1, 1] $
            ;;, xrange = [0, Sx*image_scale], yrange = [0, Sy*image_scale] $
            , aspect = Sy/float(Sx) $
            , xtitle = 'HPLN / 1"', ytitle = 'HPLT / 1"' $
            , title = 'Tile coordinates '+fxpar(hdr,'DATE-OBS')+' '+pref
    cgplot, /add, /over, diskpos[0,*], diskpos[1, *], psym=16, color='red'

  endif


  ;; Alignment
  shifts = fltarr(2, Ntiles)

  ;; Start with the center tile, then base the order of
  ;; stitching on the largest overlapping area.
  itile = 0

  xc = median(xpos)
  yc = median(ypos)
  dr = sqrt((xpos-xc)^2 + (ypos-yc)^2)
  ord = sort(dr)      
  red_progressbar, 0, Ntiles, 'Add tile '+red_stri(ord[itile])+' to the mosaic'
  
  edge = sobel(mmap[*, *, ord[itile]])
  indx = array_indices(edge, where(edge gt 0.5*max(edge)))

  if keyword_set(debug) then begin
    cgplot, /add, /over $
            , diskpos[0,ord[itile]] $
            , diskpos[1,ord[itile]] $
;            , xcp+indx[0,ord[itile]]*fac*image_scale $
;            , ycp+indx[1,ord[itile]]*fac*image_scale $
            , psym=16, color=colors[ord[itile]]
    cgtext, /add, align = 1, color = colors[ord[itile]], charsize = 1.5 $
            , diskpos[0,ord[itile]] $
            , diskpos[1,ord[itile]] $
;            , xcp+indx[0,ord[itile]]*fac*image_scale $
;            , ycp+indx[1,ord[itile]]*fac*image_scale $
            , red_stri(ord[itile]) + ' ('+red_stri(itile)+')' 
  endif
  
  ;; Add center tile
  totim = imap[*, *, ord[itile]] * wmap[*, *, ord[itile]]
  totw  = wmap[*, *, ord[itile]]
  totm  = mmap[*, *, ord[itile]]

  commonarea = fltarr(Ntiles)

  ;; Loop over remaining tiles
  for itile = 1, Ntiles-1 do begin

    if itile ne Ntiles-1 then begin
      ;; Find the tile with the most overlap with the mosaic so
      ;; far assembled
      if keyword_set(debug) then print, itile, ord[itile]
      for jmos = itile, Ntiles-1 do begin
        commonmask = mmap[*, *, ord[jmos]] * totm
        commonarea[jmos] = total(commonmask)
        if keyword_set(debug) then print, jmos, ord[jmos], commonarea[jmos]
      endfor                    ; jmos
      if keyword_set(debug) then print, ord[itile:Ntiles-1]
      if keyword_set(debug) then print, commonarea[itile:Ntiles-1]
      ;; Find max
      mxx = max(commonarea[itile:Ntiles-1], maxloc)
      
      ;; Swap to get it first of remaining
      tmp = ord[itile]
      ord[itile] = ord[maxloc+itile]
      ord[maxloc+itile] = tmp
      
    endif 

    red_progressbar, itile, Ntiles, 'Add tile '+red_stri(ord[itile])+' to the mosaic'
    
    ;; We should now add tile ord[itile] to the mosaic
    edge = sobel(mmap[*, *, ord[itile]])
    indx = array_indices(edge, where(edge gt 0.5*max(edge)))
    if keyword_set(debug) then begin
      cgplot, /add, /over $
              , diskpos[0,ord[itile]] $
              , diskpos[1,ord[itile]] $
;              , xcp+indx[0,ord[itile]]*fac*image_scale $
;              , ycp+indx[1,ord[itile]]*fac*image_scale $
              , psym=16, color = colors[ord[itile]]
      cgtext, /add, align = 1, color = colors[ord[itile]], charsize = 1.5 $
              , diskpos[0,ord[itile]] $
              , diskpos[1,ord[itile]] $
;              , xcp+indx[0,ord[itile]]*fac*image_scale $
;              , ycp+indx[1,ord[itile]]*fac*image_scale $
;              , xcp+median(indx[0,*]*fac*image_scale) $
;              , ycp+median(indx[1,*]*fac*image_scale) $
              , red_stri(ord[itile])  + ' ('+red_stri(itile)+')' 
    endif
    

    

    ;; Measure the alignment shifts of the current tile vs the mosaic
    ;; so far within the common, overlapping FOV.
    red_progressbar, itile, Ntiles, 'Add tile '+red_stri(ord[itile])+' to the mosaic : measure alignment'
    commonmask = mmap[*, *, ord[itile]] * totm
    bb = red_boundingbox(commonmask)
    if keyword_set(debug) then print, 'boundingbox : ', bb
    shifts[*, ord[itile]] = red_shc_mask((totim / (1e-5 + totw))[bb[0]:bb[2],bb[1]:bb[3]] $
                                         , imap[bb[0]:bb[2],bb[1]:bb[3], ord[itile]] $
                                         , commonmask[bb[0]:bb[2],bb[1]:bb[3]] $
                                         , range = 25, poly = 2)
    
    if keyword_set(debug) then print, shifts[*, itile]
    
    ;; Accept shifts only if they are less than a few arcsec.
    dr = sqrt(total(shifts[*, ord[itile]]^2)) * fac * image_scale
    if dr lt 10 then begin
      red_progressbar, itile, Ntiles, 'Add tile '+red_stri(ord[itile])+' to the mosaic : apply alignment'
      
      ;; No sub-pixel shifts at this point
      shifts[*, ord[itile]] = round(shifts[*, ord[itile]])

      ;; This shifts the current image within its full-size map
      imap[*, *, ord[itile]] = red_shift_im(imap[*, *, ord[itile]] $
                                            , shifts[0, ord[itile]], shifts[1, ord[itile]], missing = 0.0)
      wmap[*, *, ord[itile]] = red_shift_im(wmap[*, *, ord[itile]] $
                                            , shifts[0, ord[itile]], shifts[1, ord[itile]], missing = 0.0) >0
      mmap[*, *, ord[itile]] = red_shift_im(mmap[*, *, ord[itile]] $
                                            , shifts[0, ord[itile]], shifts[1, ord[itile]], missing = 0.0) $
                               * (wmap[*, *, ord[itile]] gt 0)

      ;; Shift also the corner pixel coordinates. This means we don't
      ;; have to shift the NB images, we just have to put them into
      ;; their maps according to the updated corner cordinates.
      llx[ord[itile]] += shifts[0, ord[itile]]
      urx[ord[itile]] += shifts[0, ord[itile]]
      lly[ord[itile]] += shifts[1, ord[itile]]
      ury[ord[itile]] += shifts[1, ord[itile]]

    endif else stop


    if ~keyword_set(no_destretch) then begin
      red_progressbar, itile, Ntiles, 'Add tile '+red_stri(ord[itile])+' to the mosaic : measure stretch'
      
      ref = totim[llx[ord[itile]] : urx[ord[itile]], lly[ord[itile]] : ury[ord[itile]]] $
            / totw[llx[ord[itile]] : urx[ord[itile]], lly[ord[itile]] : ury[ord[itile]]]
      
      im = imap[llx[ord[itile]] : urx[ord[itile]], lly[ord[itile]] : ury[ord[itile]], ord[itile]]

      red_missing, ref, /inplace, missing_type_wanted = 'median'
      red_missing, im,  /inplace, missing_type_wanted = 'median'

      ref *= commonmask[llx[ord[itile]] : urx[ord[itile]], lly[ord[itile]] : ury[ord[itile]]]
      im  *= commonmask[llx[ord[itile]] : urx[ord[itile]], lly[ord[itile]] : ury[ord[itile]]]

      red_missing, ref, /inplace, missing_type_wanted = 'median'
      red_missing, im,  /inplace, missing_type_wanted = 'median'

      if n_elements(grids) eq 0 then begin
        thisgrid = rdx_cdsgridnest(ref/median(ref), im/median(im), gridtiles, gridclips)
        dims = size(thisgrid, /dim)
        grids = fltarr([dims, Ntiles])
        grids[*, *, *, ord[itile]] = thisgrid
      endif else begin
        grids[*, *, *, ord[itile]] = rdx_cdsgridnest(ref/median(ref), im/median(im), gridtiles, gridclips)
      endelse
        
      red_progressbar, itile, Ntiles, 'Add tile '+red_stri(ord[itile])+' to the mosaic : apply stretch'
      im = rdx_cstretch(imap[llx[ord[itile]] : urx[ord[itile]], lly[ord[itile]] : ury[ord[itile]], ord[itile]] $
                        , grids[*, *, *, ord[itile]], nthreads=nthreads)

      imap[llx[ord[itile]] : urx[ord[itile]], lly[ord[itile]] : ury[ord[itile]], ord[itile]] = im
      
    endif

    
    totim += imap[*, *, ord[itile]] * wmap[*, *, ord[itile]]
    totw  += wmap[*, *, ord[itile]]
    totm  OR= mmap[*, *, ord[itile]]

  endfor                        ; itile


  red_message, 'Tiles were added in this order: '+strjoin(ord, ', ')

  
  indx = where(totm gt 0)
  wbmosaic = fltarr(Sx/fac, Sy/fac)
  wbmosaic[indx] = totim[indx]/totw[indx]
;  red_missing, wbmosaic, /inplace, missing_type_wanted = 'nan'
  
  bb = red_boundingbox(wbmosaic gt 0)

  wbmosaic = wbmosaic[bb[0]:bb[2],bb[1]:bb[3]] >0
;  wbmosaic = red_fillpix(wbmosaic, nthreads=nthreads) $
;             * totm[bb[0]:bb[2],bb[1]:bb[3]]
  red_missing, wbmosaic, /inplace, missing_type_wanted = 'nan'

;  indx = where(wbmosaic gt median(wbmosaic)/100)
;  tmp = red_histo_opt(wbmosaic[indx], cmin = cmin, cmax = cmax)

;  red_show, bytscl(wbmosaic >cmin <cmax) ;, /scroll
  red_show, wbmosaic, /opt, /scroll

  dims = [size(wbmosaic, /dim), Ntuning, Nstokes, 1]

  ;; WCS coordinates
  wcs = replicate({ wave:dblarr(2,2) $
                    , hplt:dblarr(2,2) $
                    , hpln:dblarr(2,2) $
                    , time:dblarr(2,2) $
                  }, Ntuning)
  wcs.wave = wcs_alltiles[*, 0].wave
  for ituning = 0, Ntuning-1 do begin
    wcs[ituning].time = mean(wcs_alltiles[ituning, *].time)
    wcs[ituning].hpln = reform(mean(hpln) + image_scale*float(dims[0])*0.5*[-1,1,-1,1],2,2)
    wcs[ituning].hplt = reform(mean(hplt) + image_scale*float(dims[1])*0.5*[-1,-1,1,1],2,2)
  endfor                        ; ituning


  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
     , prstep = "MOSAICKING,SPATIAL-ALIGNMENT,DESTRETCHING" $
     , prpara = prpara $
     , prproc = inam

  ;; Add some info to the hdr
  date_keyword = red_timestamp(/utc, /iso)
  red_fitsaddkeyword, hdr, 'DATE', date_keyword
  anchor = 'DATE'

  red_fitsaddkeyword, hdr, 'DATE-BEG', self.isodate + 'T' + red_timestring(min(wcs.time))
  red_fitsaddkeyword, hdr, 'DATE-AVG', self.isodate + 'T' + red_timestring(mean(wcs.time))
  red_fitsaddkeyword, hdr, 'DATE-END', self.isodate + 'T' + red_timestring(max(wcs.time))

  
  ;; Initialize fits file, set up for writing the data part.
  self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims, wcs = wcs

  ;; Now make the NB mosaics and add them to the fitscube. We don't
  ;; have to put the tiles into their maps in the same order as the WB
  ;; tiles.  

  for ituning = 0, Ntuning-1 do begin
    for istokes = 0, Nstokes-1 do begin

      print
      print, '(ituning,istokes) = ('+red_stri(ituning)+','+red_stri(istokes)+')'
      
      ;; Initialize the image map
      imap *= 0.0
      
      for itile = 0, Ntiles-1 do begin
        red_progressbar, itile, Ntiles, 'NB tile '+red_stri(itile)
        
        red_fitscube_getframe, sfiles[itile], nbim, ituning = ituning, istokes = istokes

        fxread, sfiles[itile], wbim, extension = 'WBIMAGE'

        ;; Now do similar operations to the nb tile as done to the wb
        ;; tile above

        ;; Some fixing with the areas where data is missing
        red_missing, nbim, /inplace, missing_type_wanted = 'NaN'
        nbim = red_centerpic(nbim, xs = Nx, ys = Ny, z = !Values.F_NaN)
        red_missing, nbim, /inplace, missing_type_wanted = 'median'

        ;; We adjusted the corner coordinates [llx,lly,urx,ury] with
        ;; respect to the alignment shifts meassured in the WB above.
        ;; So we don't have to shift the NB image maps, we just have
        ;; to put the NB tiles into their maps according to the corner
        ;; coordinates.
        if fac eq 1 then begin
          imap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] = mask * nbim
        endif else begin
          imap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] = congrid( masks[*, *, itile] * nbim $
                                                                                   , Nx/fac, Ny/fac, cubic = -0.5 )
        endelse

        
;        dr = sqrt(total(shifts[*, itile]^2)) * fac * image_scale
;        if dr lt 10 then begin
;          ;; Shift image map
;          imap[*, *, itile] = red_shift_im(imap[*, *, itile] $
;                                           , shifts[0, itile], shifts[1, itile], missing = 0.0)
;        endif

        ;; We could have an option to do individual destretch to the
        ;; NB tiles rather than using the stretch grids from WB.   <-------------------- !!!!!!
        ;; But the first test dataset looks good with the WB grids.
        if ~keyword_set(no_destretch) then begin
    
          case !true of 

            keyword_set(individual_destretch) : begin
              ;; Calculate new stretch grids for each tuning
              stop
            end

            else : begin
              ;; Use stretch grids from WB
              im = rdx_cstretch(imap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] $
                                , grids[*, *, *, itile], nthreads=nthreads)
              imap[llx[itile] : urx[itile], lly[itile] : ury[itile], itile] = im
            end

          endcase

        endif
        
        ;;      totim += imap[*, *, itile] * wmap[*, *, itile]
        
      endfor                    ; itile    

      totim = totm * total(imap*wmap, 3) / totw
      
      indx = where(totm gt 0)
      
      nbmosaic = fltarr(Sx/fac, Sy/fac)
      nbmosaic[indx] = totim[indx]
      red_missing, nbmosaic, /inplace, missing_type_wanted = 'nan'
      nbmosaic = nbmosaic[bb[0]:bb[2],bb[1]:bb[3]] 
;;      nbmosaic = red_fillpix(nbmosaic, nthreads=nthreads) $
;;                 * totm[bb[0]:bb[2],bb[1]:bb[3]]
      
;      red_show, nbmosaic, /opt  ;, /scroll

      ;; Add frame to file
      red_fitscube_addframe, fileassoc, nbmosaic, ituning = ituning, istokes = istokes

      
    endfor                      ; istokes
  endfor                        ; ituning

  ;; Close fits file 
  self -> fitscube_finish, lun, wcs = wcs, direction = direction

  ;; Add the WB mosaic as an extension as in other scan cubes.
  ehdr=hdr
  fxaddpar, ehdr, 'XTENSION', 'IMAGE'
  sxdelpar, ehdr, 'SIMPLE'
  check_fits, wbmosaic, ehdr, /update
  fxaddpar, ehdr, 'DATE', date_keyword
  anchor = 'DATE'
  red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'WBIMAGE', 'Wideband image'
  red_fitsaddkeyword, anchor = anchor, ehdr, 'PCOUNT', 0
  red_fitsaddkeyword, anchor = anchor, ehdr, 'GCOUNT', 1
  writefits, filename, wbmosaic, ehdr, /append


  ;; Done with this scan.
  print, inam + ' : Narrowband scan cube stored in:'
  print, filename



end
