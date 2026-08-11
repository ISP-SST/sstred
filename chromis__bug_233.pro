;; Fix chromis NB cubes that were hit by issue 233, see https://dubshen.astro.su.se/redmine/issues/233

;; For this to work, the WB cube used in the call to make_nb_cube still has to exist.

;; Usage:
;; IDL> a = chromisred()
;; IDL> a -> bug_233, filename
;; where filename is the name of a NB cube that needs to be fixed.

pro chromis::bug_233, filename

  inam = red_subprogram(/low, calling = inam1)

  if ~file_test(filename) then begin
    print, 'The file does not exist!'
    print, filename
    retall
  endif

  red_make_prpara, prpara, filename

  ;; We do currently not correct for the small scale cavity map in
  ;; CHROMIS data. (We should get this from earlier meta data!)
  remove_smallscale = 0      


  ;; Collect info from the file header
  hdr = headfits(filename)

  ;; Dimensions
  Naxis = fxpar(hdr, 'NAXIS*')
  Nx      = Naxis[0]
  Ny      = Naxis[1]
  Nwav    = Naxis[2]
  Nstokes = Naxis[3]
  Nscans  = Naxis[4]
  
  ;; Information about processing steps in the formation of the input
  ;; file. 
  prprocs = fxpar(hdr, 'PRPROC*')
  prparas = fxpar(hdr, 'PRPARA*')

  
  ;; Check that the file has the bug. (And the extension!) And
  ;; that it's a CHROMIS file!

  red_fitscube_getframe, filename, im, iframe = 0
  red_fitscube_getwcs, filename, coordinates = coordinates, distortions=distortions

  red_show, im, /scroll, title = 'Image'
  red_show, distortions.wave[*,*,0,0,0], /scroll, title = 'Cavity map'

  s = ''
  print
  print, "Do the position and orientation of the cavity map in its array match"
  read,  "that of the image (it should be obvious if they don't) [Y/n]? ", s
  if strupcase(strmid(s, 0, 1)) ne 'N' then begin
    print
    print, 'Ok, then no correction is needed.'
    print
    return
  endif
  
  ;; Get the name of the WB cube file from the prpara info
  pos_makenb   = where(strmatch(prprocs, '*make_nb_cube'  ), Nmakenb  )
  pos_makescan = where(strmatch(prprocs, '*make_scan_cube'), Nmakescan)

  case 1 of
    Nmakenb gt 0 : begin        ; This is a NB cube
      cube_paras = prparas[pos_makenb[0]]
      cube_paras_struct = json_parse(cube_paras, /tostruct)
      wcfile = cube_paras_struct.wcfile
      wbhdr = headfits(wcfile)
;      wbpref = strtrim(fxpar(wbhdr, 'FILTER1'), 2)
      fxbopen, bunit, cube_paras_struct.wcfile, 'MWCINFO', bbhdr
      fxbreadm, bunit, row = 1 $
                , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
                ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01,   direction
      fxbclose, bunit
      if (json_parse(prparas[pos_makenb[0]])).haskey('CMAP_FWHM') then begin
        fwhm = (json_parse(prparas[pos_makenb[0]]))['CMAP_FWHM']
      endif else begin
        fwhm = 7.0
      endelse

    end
    Nmakescan gt 0 : begin      ; This is a SCAN cube

      stop                      ; Not implemented yet!
      
      cube_paras = prparas[pos_makescan[0]]
      cube_paras_struct = json_parse(cube_paras, /tostruct)
;      wbpref = (stregex(cube_paras_struct.dir $
;                        , '/([0-9][0-9][0-9][0-9])/', /extract,/subex))[1]
      wbimage = mrdfits(filename, 'WBIMAGE', ehdr, STATUS=status, /silent)
      wcTMEAN = median(wbimage)
    end
    else : begin
      print, inam+' : This type of cube is not implemented yet.'
      stop
    end
  endcase

  

  ;; Read the header from the corrected WB cube. Variables begin with
  ;; WC for Wideband Cube. 
  if ~file_test(wcfile) then begin
    print, 'WB cube missing, please run make_wb_cube.'
    print, 'If you deleted your momfbd output you will have to rerun the momfbd processing.'
    print, wcfile
    retall
  endif
  wchead = red_readhead(wcfile)
;  ;; Read parameters from the WB cube
;  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
;  fxbreadm, bunit, row = 1 $
;            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
;            ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01,   direction
;  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
;  ;; wbgfiles (WideBand Global).
;  fxbread, bunit, wbgfiles, 'WFILES', 1
;  fxbclose, bunit

  ;; Default for wb cubes without direction parameter
  if n_elements(direction) eq 0 then direction = 0
  
  x0 = wcX01Y01[0]
  x1 = wcX01Y01[1]
  y0 = wcX01Y01[2]
  y1 = wcX01Y01[3]
  origNx = x1 - x0 + 1
  origNy = y1 - y0 + 1

  ;; Create cubes for science data and scan-adapted cavity maps.
  cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)


  ;; What NB prefilters were used? That is actually not in the NB cube
  ;; header! Reconstruct this info from the prefilters defined in
  ;; linedef.py and the wavelength coordinates. 
  linedefdir = 'downloads/'
  ldfiles = file_search(linedefdir + '/linedef.*', count = nLD)
  ;; Find latest file that is before date, for now use first file
  ild = 0
  linedef = red_readlinedef(ldfiles[ild])
  all_prefs = linedef.nb_wl
  
  lambda = coordinates[*,0].wave[0,0]

  prefs = lonarr(Nwav)
  for iwav = 0, Nwav-1 do begin
    tmp = min(abs(lambda[iwav] - all_prefs/10.), ml)
    prefs[iwav] = all_prefs[ml]
  endfor                        ; iwav

  ;; Read the original cavity map
  pindx = where(prefs ne '3999') ; No cavity map for the Ca II H continuum
  pindx = pindx[uniq(prefs[pindx], sort(prefs[pindx]))]
  cprefs = prefs[pindx]
  Ncprefs = n_elements(cprefs)

  texp = string(fxpar(hdr, 'TEXPOSUR')*1000., format = '(f5.2)')+'ms'

  for icprefs = 0, Ncprefs-1 do begin
    cfile = file_search('flats/spectral_flats/*' $
                        + '_' + texp + '_G*_' $
                        + string(cprefs[icprefs], format = '(i04)') $
                        + '_fit_results.sav' $
                        , count = Nfiles)

    if Nfiles ne 1 then stop
    
    if ~file_test(cfile) then begin
      print, inam + ' : Error, calibration file not found -> '+cfile
      print, 'Please run the fitprefilter for '+cprefs[icprefs]+' or continue without'
      print, 'cavity map for '+cprefs[icprefs]
      stop
    endif
    restore, cfile                   ; The cavity map is in a struct called "fit". 
    cmap = reform(fit.pars[1,*,*])   ; Unit is [Angstrom]
    ;;  cmap = rotate(temporary(cmap), direction)
    cmap /= 10.                 ; Make it [nm]
    cmap = -cmap                ; Change sign so lambda_correct = lambda + cmap
    fit = 0B                    ; Don't need the fit struct anymore.
    
    if keyword_set(remove_smallscale) then begin
      ;; If the small scale is already corrected, then include only the
      ;; low-resolution component in the metadata. The blurring kernel
      ;; should match how the low resolution component was removed when
      ;; making flats.
      npix = 30                 ; Can we get this parameter from earlier headers?
      cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
      cpsf /= total(cpsf, /double)
      cmap = red_convolve(temporary(cmap), cpsf)
      cmap1 = cmap
    endif else begin
      ;; If the small scale is not already corrected, then we still want
      ;; to blur the cavity map slightly.
      npsf = round(fwhm * 7.)
      if((npsf/2)*2 eq npsf) then npsf += 1L
      psf = red_get_psf(npsf, npsf, fwhm, fwhm)
      psf /= total(psf, /double)
      ;; Leave the orignal cmap alone, we might need it later.
      cmap1 = red_convolve(cmap, psf)
    endelse

    ;; Read the output of the pinhole calibrations so we can do the same
    ;; to the cavity maps as was done to the raw data in the momfbd
    ;; step. This output is in a struct "alignments" in the save file
    ;; 'calib/alignments.sav'
    restore,'calib/alignments.sav'
    ;; Should be based on state1 or state2 in the struct? make_cmaps
    ;; says "just pick one close to continuum (last state?)".
    indx = where(string(cprefs[icprefs], format = '(i04)') eq alignments.state2.prefilter, Nalign)
    case Nalign of
      0    : stop               ; Should not happen!
      1    : amap = invert(      alignments[indx].map           )
      else : amap = invert( mean(alignments[indx].map, dim = 3) )
    endcase
    cmap1 = rdx_img_project(amap, cmap1, /preserve) ; Apply the geometrical mapping
    cmap1 = red_rotate(cmap1, direction)            ; Correct orientation
    cmap1 = cmap1[x0:x1,y0:y1]                      ; Clip to the selected FOV

    
    ;; Now make rotated copies of the cavity map
    for iscan = 0L, Nscans-1 do begin

      ;; Apply the same derot, align, dewarp as for the science data
      cmap11 = red_rotation(cmap1, ang[iscan], $
                            wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF, $
                            stretch_grid = reform(wcGRID[iscan,*,*,*]))
      ;;cmap11 = red_stretch_linear(temporary(cmap11), reform(wcGRID[iscan,*,*,*]))
      
      cavitymaps[0, 0, 0, 0, iscan] = cmap11

    endfor                      ; iscan

    if icprefs eq 0 or max(abs(cavitymaps)) gt fxpar(hdr, 'CWERR3') then begin
      fxaddpar, hdr, 'CWERR3', max(abs(cavitymaps)), '[nm] Max total distortion'
    endif
    
    tindx = where(prefs eq cprefs[icprefs], Nt)

    fits_open, filename, fcb, /update
    exten_no = where(fcb.extname eq 'WCSDVARR', Nexten)
    if Nexten ne Ncprefs then stop
    modfits, fcb, cavitymaps, exten_no = exten_no[icprefs], errmsg=errmsg
    if n_elements(errmsg) gt 0 then stop
    fits_close,fbc
    
  endfor                        ; icprefs

  self -> headerinfo_addstep, hdr $
                              , prstep = 'BUGFIX' $
                              , prpara = prpara $
                              , prproc = inam

  red_fitscube_newheader, filename, hdr

end

;; Usage:

;; Obviously this should be the name of your nb cube file:
filename = 'cubes_nb_mats/nb_3950_2020-07-14T08:40:49_scans=100-300_corrected_im.fits'

;; Do the correction:
a=chromisred(/dev)
a -> bug_233, filename

end
