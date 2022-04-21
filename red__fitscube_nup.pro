; docformat = 'rst'

;+
; Fix the orientation and rotation of old fitscubes.
;
; This is for fitscubes made before we implemented the rotation and
; orientation setup parameters and with the average derotation angle
; subtracted. The goal is to produce a cube that is oriented along
; HPLN/HPLT coordinate axes with N up as they normally are done.
; 
; The direction and orientation parameters, whether from config.txt or
; given as keywords, are used to determine how the cube *should* be
; rotated. The mirrorx and mirrory keywords on the other hand, or if
; they are missing possibly info from a momfbd config file, are used
; to determine how it is *currently* oriented.
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
; :Params:
; 
;    inname : in, type=string
; 
;       The path to the file to be converted.
; 
; :Keywords:
;
;    direction : in, optional, type=integer, default="from config file"
;
;      The relative orientation of reference cameras of different
;      instruments. Note that only if direction is set in the config
;      file will it be assumed to be correct when setting the CSYERR
;      FITS header keyword.
;
;    mirrorx : in, optional, type=boolean
;
;      The cube wb data are mirrored in the X direction compared to
;      the raw data.
;
;    mirrory : in, optional, type=boolean
;
;      The cube wb data are mirrored in the Y direction compared to
;      the raw data.
;
;    outname : in, out, optional, type=string, default = inname+'.fits'
;
;      Where to write the output. A spectral cube might also be
;      written, the file name for this will be generated based on
;      outname. If a variable, will be returned with the actual file
;      name.
;
;    overwrite : in, optional, type=boolean
;
;       Don't care if cube is already on disk, overwrite it
;       with a new version.
;
;    point_id : in, optional, type=string, default="YYYY-MM-DDThh:mm:ss"
;
;      Value for the POINT_ID header keyword. 
;
;    rotation : in, optional, type=float
;
;      Offset angle to be added to the field rotation angles.
;
;    do_wb_cube : in, optional, type=boolean
;
;      Set to correct WB cube as well.
;
;    wcs_improve_spatial  : in, optional, type=boolean
;
;      Set to align FOV of a fitscube with HMI intensity image.
;      Use when there are prominent features in FOV.
;
;    no_checksum : in, optional, type=boolean
;
;      Do not calculate DATASUM and CHECKSUM.
;
;    release_date : in, optional, type=string
;
;      The value of the RELEASE keyword, a date after which the data
;      are not proprietary.
;
;    release_comment : in, optional, type=string
;
;      The value of the RELEASEC keyword, a comment with details about
;      the proprietary state of the data, and whom to contact.
;
;    feature : in, optional, type=string
;
;      The value of FEATURE keyword, names of solar features in the FOV.
;
;    observer  : in, optional, type=string
;
;      The value of the OBSERVER keyword, names of observers.
; 
; 
; :History:
; 
;    2022-03-29 : MGL. First version.
;
;    2022-04-11 : OA. Added recursive call for wb cube, keywords for
;                 finalizing the header, time-dependent solar
;                 coordinates in WCS.
; 
;-
pro red::fitscube_nup, inname  $
                       , direction = direction $
                       , wcdirection = wcdirection $
                       , rotation = rotation $
                       , flip = flip $
                       , do_wb_cube = do_wb_cube $
                       , wcs_improve_spatial = wcs_improve_spatial $
                       , mirrorx = mirrorx $                             
                       , mirrory = mirrory $                             
                       , outdir = outdir $
                       , outname = outname $                                    
                       , overwrite = overwrite $
                       , no_checksum = no_checksum $
                       , release_date = release_date $
                       , release_comment = release_comment $
                       , feature = feature $
                       , observer = observer $
                       , point_id = point_id $
                       , nthreads = nthreads $
                       , status = status

  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  
  
  if n_elements(direction) eq 0 then direction = self.direction
  if n_elements(rotation)  eq 0 then rotation  = self.rotation

  ;; May want to check the ALIGN_CLIP keyword if the cfg files are
  ;; still around. In case there was some mirroring done to the WB
  ;; object (the first one). Then the (default) direction might have
  ;; to be changed.
  
  ;; Alternatively, read a raw WB image and compare its orientation to
  ;; the input cube!
  
  red_make_prpara, prpara, inname
  red_make_prpara, prpara, direction
  red_make_prpara, prpara, flip
  red_make_prpara, prpara, mirrorx
  red_make_prpara, prpara, mirrory
  red_make_prpara, prpara, outdir
  red_make_prpara, prpara, outname
  red_make_prpara, prpara, point_id
  red_make_prpara, prpara, rotation
  red_make_prpara, prpara, no_checksum
  red_make_prpara, prpara, do_wb_cube
  red_make_prpara, prpara, wcs_improve_spatial  
  
  if n_elements(outdir) eq 0 then outdir = 'cubes_converted/'
  file_mkdir, outdir
 
  ;; Cameras and detectors
  self->getdetectors
  wbindx = where(strmatch(*self.cameras,'*-W')) ; Matches CRISP and CHROMIS WB cameras
  wbcamera = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx = where(strmatch(*self.cameras,'*-[NRT]')) ; Matches CRISP and CHROMIS NB cameras
  nbcamera = (*self.cameras)[nbindx[0]]
  nbdetector = (*self.detectors)[nbindx[0]]

  ;; Read header and fix PRSTEP info 
  red_fitscube_correct_prstep, inname, header = inhdr, /nowrite
  
  dims = fxpar(inhdr, 'NAXIS*')
  Nx       = dims[0]  
  Ny       = dims[1]  
  Ntunings = dims[2]  
  Nstokes  = dims[3]  
  Nscans   = dims[4]
  bitpix = fxpar(inhdr,'BITPIX')
  
  ;; NB cube of WB cube?
  filetype = (strsplit(file_basename(inname),'_',/extract))[0]

  ans=''
  no_wcfile = 0B
  if filetype eq 'nb' then begin
    indx = where(strmatch(inhdr, '*WCFILE*'), Nwhere)
    if Nwhere eq 0 then begin
      print, 'There is no information about wb cube in the header.'
      read, 'Do you want to continue (y/N) : ',ans
      if strupcase(ans) ne 'Y' then $
        return $
      else $
        no_wcfile = 1B
    endif else begin
      par = strmid(inhdr[indx],0,7)
      par_key = red_fitsgetkeyword(inhdr, par[0])
      wcfile = (json_parse(par_key))['WCFILE']
      if ~file_test(wcfile) then begin
        print, inam + ' : Could not find the WB cube -- ', wcfile
        read, 'Do you want to continue (y/N) : ',ans
        if strupcase(ans) ne 'Y' then $
          return $
        else $
          no_wcfile = 1B
      endif
    endelse
    if ~no_wcfile then begin
      wb_hdr = headfits(wcfile)
      wb_type = fxpar(wb_hdr,'BITPIX')
      if wb_type ne 16 then begin
        ;; With old nb cubes we need to rely on integer wb cubes for padding
        ;; Because old padding might be not really same [median] value
        ;; after rotation and red__missing fails.
        print,'Corresponding WB cube is not integer (i.e. not very old).'
        print,'This method might be not optimal for curation.'
        read,'Do you want to continue (y/N) : ',ans
        if strupcase(ans) ne 'Y' then return
      endif
    endif
  endif else wcfile = inname
  
  ;; Prefilter
  pref = strtrim(red_fitsgetkeyword(inhdr, 'FILTER1'), 2)
  
  ;; Scan numbers
  tmp = red_fitsgetkeyword(inname, 'SCANNUM', variable_values = variable_values)
  s_array = reform(variable_values.values)

  ;; POINT_ID
  pid = strtrim(red_fitsgetkeyword(inhdr, 'POINT_ID', count=cnt), 2)
  if cnt eq 0 then pid = strtrim(red_fitsgetkeyword(inhdr, 'DATE-OBS', count=cnt), 2)
  if cnt eq 0 then pid = strtrim(red_fitsgetkeyword(inhdr, 'STARTOBS', count=cnt), 2)
  if cnt ne 0 then point_id = pid else begin
    print, inam + ' : No POINT_ID info in file header, please call again with point_id keyword.'
    return
  endelse

  ;; Time stamp
  timestamp = (strsplit(point_id, 'T', /extract))[1]

  if n_elements(point_id) eq 0 then point_id = self.isodate + 'T' + timestamp

  if n_elements(outname) eq 0 then begin
    
    ;; Find old datatags in infile. They are all segments after the
    ;; last segment that ends with a digit.
    fsplt = strsplit(file_basename(inname, '.fits'), '_', /extract)
    indx = where(strmatch(fsplt,'*[0-9]'), Nwhere)
    if Nwhere lt n_elements(fsplt) then begin
      ;; We want all those datatags in the output filename.
      datatags = fsplt[indx[-1]+1:-1]
      ;; We also want to add 'nup'.
      if n_elements(datatags) eq 1 then datatags = ['nup', datatags] $
      else datatags = [datatags[0:-2], 'nup', datatags[-1]]
    endif

    oname = outdir + '/' $
            + red_fitscube_filename(filetype $
                                    , pref $
                                    , timestamp $
                                    , red_collapserange(s_array, /nobrack) $
                                    , point_id $ 
                                    , datatags = datatags)
  endif else begin
    oname = outname
  endelse

  ;; Already done?
  if file_test(oname) then begin
    if keyword_set(overwrite) then begin
      print, 'Overwriting existing data cube:'
      print, oname
    endif else begin
      print, 'This data cube exists already:'
      print, oname
      status = 0
      return
    endelse
  endif

  
  ;; Establish the orientation of the file
  if n_elements(mirrorx) eq 0 or n_elements(mirrory) eq 0 then begin

    cfgfiles = file_search('*mfbd*/'+timestamp+'/'+pref+'/cfg/*cfg', count = Ncfg)  

    if Ncfg eq 0 then begin
      print, inam + ' : Cannot find a matching momfbd config file.'
    endif else begin

      print, inam + ' : Found a matching momfbd config file: '
      print, cfgfiles[0]
      
      ;; Check ALIGN_CLIP for mirroring of the anchor object       
      align_clip = redux_cfggetkeyword(cfgfiles[0],'object0.channel0.align_clip')
      if n_elements(align_clip) ne 0 then begin
        align_clip = long(strsplit(align_clip, ',', /extract))
        xmirrorx = align_clip[0] gt align_clip[1]    
        xmirrory = align_clip[2] gt align_clip[3]    
      endif else print, inam + ' : No ALIGN_CLIP in momfbd cfg file'
      ;; Check ALIGN_MAP for mirroring of the anchor object 
      align_map = redux_cfggetkeyword(cfgfiles[0],'object0.channel0.align_map')
      if n_elements(align_map) ne 0 then begin
        align_map = float(reform(strsplit(align_map,',',/extract),3,3))
        xmirrorx = align_map[0, 0] lt 0          
        xmirrory = align_map[1, 1] lt 0          
      endif else print, inam + ' : No ALIGN_MAP in momfbd cfg file'
      ;; There shouldn't be any cfg files with both ALIGN_CLIP and ALIGN_MAP!

      if n_elements(mirrorx) gt 0 then begin
        print, inam + ' : keyword mirrorx on input : ', mirrorx
      endif else begin
        if n_elements(xmirrorx) gt 0 then begin
          mirrorx = xmirrorx
          print, inam + ' : keyword mirrorx from momfbd cfg : ', mirrorx
        endif
      endelse 

      if n_elements(mirrory) gt 0 then begin
        print, inam + ' : keyword mirrory on input : ', mirrory
      endif else begin
        if n_elements(xmirrory) gt 0 then begin
          mirrory = xmirrory
          print, inam + ' : keyword mirrory from momfbd cfg : ', mirrory
        endif
      endelse 

    endelse

  endif  

  if n_elements(wcdirection) eq 0 and $
    n_elements(mirrorx) eq 0 and $
    n_elements(mirrory) eq 0 then begin

    print,'There is no information about previous 90 degrees rotations or mirroring of the cube.'
    read, 'Do you want to continue (Y/n): ',ans
    if strupcase(ans) eq 'N' then return
  endif
  if  n_elements(wcdirection) eq 0 then wcdirection = 0

  ;; Make header
  outhdr = inhdr
  red_fitsaddkeyword, outhdr, 'BITPIX', -32
  red_fitsaddkeyword, outhdr, 'NAXIS1', Nx
  red_fitsaddkeyword, outhdr, 'NAXIS2', Ny
  
  if keyword_set(do_wb_cube) or keyword_set(wcs_improve_spatial) then begin
    if no_wcfile then begin
      print, "There is no wb cube. Please rerun without 'do_wb_cube' or 'wcs_improve_spatial' options."
      return
    endif
    print,'==============================='
    print,'Processing WB cube first.'
    print,'==============================='
    if ~keyword_set(outname) then $
      wb_oname = outdir + '/' $
            + red_fitscube_filename('wb' $
                                    , pref $
                                    , timestamp $
                                    , red_collapserange(s_array, /nobrack) $
                                    , point_id $ 
                                    , datatags = datatags) $
    else begin
      prfx = strmid(outname,0,2)
      if prfx eq 'nb' then $
        wb_oname = outdir + '/wb' + strmid(outname,2,strlen(outname)) $
      else $
        wb_oname = outdir + '/wb_' + outname
    endelse
    
    red_fitsaddkeyword, outhdr, key $
                         , 'Align reference: '+wb_oname $
                         , 'WB cube file name'
    

    if keyword_set(wcs_improve_spatial) then no_check_sum = 1b else no_check_sum = keyword_set(no_checksum)

    self -> fitscube_nup, wcfile $
                          , direction = direction $
                          , rotation = rotation $
                          , mirrorx = mirrorx $                             
                          , mirrory = mirrory $                             
                          , outdir = outdir $
                          , outname = wb_oname $                                    
                          , overwrite = overwrite $   
                          , no_checksum = no_check_sum $
                          , release_date = release_date $
                          , release_comment = releasec $
                          , feature = feature $
                          , observer = observer $
                          , point_id = point_id $
                          , status = status
    if ~status then return
    if keyword_set(wcs_improve_spatial) then begin
      self -> fitscube_wcs_improve_spatial, wb_oname
      if ~keyword_set(no_checksum) then red_fitscube_checksums,wb_oname
    endif
  endif

  if filetype eq 'wb' then $
    red_fitscube_getwcs, inname, coordinates = wcs $
  else $
    red_fitscube_getwcs, inname, distortions = distortions, coordinates = wcs
     
  ;; Original derotation angles
  t_array = reform(wcs[0,*].time[0,0])
  ang = red_lp_angles(t_array, self.isodate, /from_log, offset_angle = defrotation) ; [radians]

  if ~no_wcfile then begin
    ;; Read parameters from the WB cube
    fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
    fxbreadm, bunit, row = 1 $
              , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
              , wcANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01, wcdirection  
    fxbclose, bunit
  endif else begin
    ;; assume that mean angle was subtracted
    mang = median(ang)
    wcang = ang - mang
  endelse
  
  ;; Angles to use are original derotation angles minus rotation
  ;; already done. This should be a scalar.
  ang = mean(ang - wcang) - (rotation-defrotation)*!dtor
  maxangle = abs(ang)
  ff = [maxangle, 0., 0., 0., 0., ang]

  red_fitscube_open, inname, in_fileassoc, in_fitscube_info
  ;; Read one frame
  red_fitscube_getframe, in_fileassoc, frame, iframe = 0, fitscube_info = in_fitscube_info
   
  ;; Get the frame with size implied by ff
  if wcdirection ne 0 then frame = rotate(temporary(frame), -wcdirection) ; Undo old direction
  if keyword_set(mirrorx) then frame = reverse(frame, 1, /over)        
  if keyword_set(mirrory) then frame = reverse(frame, 2, /over)
                                   
  frame = rotate(temporary(frame), direction)    ; Do the new direction
  frame = red_rotation(frame, full = ff, 0, 0, 0, nthreads = nthreads)

  ;; Set the spatial dimensions of the output file.
  dims[0:1] = size(frame, /dim)
  Nx = dims[0]  
  Ny = dims[1]    

  if keyword_set(wcs_improve_spatial) then begin
    red_fitscube_getwcs,  wb_oname, coordinates = wb_wcs
    for iscan = 0L, Nscans-1 do begin
      wcs[*, iscan].hpln = wb_wcs[0,iscan].hpln
      wcs[*, iscan].hplt = wb_wcs[0,iscan].hplt
    endfor                      ; iscan
  endif else begin
    
    ;; We need to restore coordinates from log files
    red_logdata, self.isodate, t_array, diskpos = pointing

    if no_wcfile then begin
      hpln = pointing[0,*]
      hplt = pointing[1,*]
    endif else begin
      hpln = pointing[0,*] - double(self.image_scale) * wcSHIFT[0,*]
      hplt = pointing[1,*] - double(self.image_scale) * wcSHIFT[1,*]
    endelse

    ;; Let's smooth coordinates.
    dt = (t_array[-1] - t_array[0]) / 60. ; minutes
    if dt le 15. or Nscans le 3 then fit_expr = 'P[0] + X*P[1]'
    if dt gt 15. and Nscans gt 3 then fit_expr = 'P[0] + X*P[1] + X*X*P[2]'
    ;; We have to exclude NaNs before running mpfitexpr.
    indx = where(finite(hpln))
    lln = hpln[indx]
    tt = t_array[indx]
    pp = mpfitexpr(fit_expr, tt, lln)
    hpln = red_evalexpr(fit_expr, t_array, pp)
    indx = where(finite(hplt))
    llt = hplt[indx]
    tt = t_array[indx]
    pp = mpfitexpr(fit_expr, tt, llt)
    hplt = red_evalexpr(fit_expr, t_array, pp)
  
    ;; But what we want to tabulate is the pointing in the corners of
    ;; the FOV. Assume hpln and hplt are the coordinates of the center
    ;; of the FOV.
    for iscan = 0L, Nscans-1 do begin
      wcs[*, iscan].hpln[0, 0] = hpln[iscan] - double(self.image_scale) * (Nx-1)/2.d
      wcs[*, iscan].hpln[1, 0] = hpln[iscan] + double(self.image_scale) * (Nx-1)/2.d
      wcs[*, iscan].hpln[0, 1] = hpln[iscan] - double(self.image_scale) * (Nx-1)/2.d
      wcs[*, iscan].hpln[1, 1] = hpln[iscan] + double(self.image_scale) * (Nx-1)/2.d
      
      wcs[*, iscan].hplt[0, 0] = hplt[iscan] - double(self.image_scale) * (Ny-1)/2.d
      wcs[*, iscan].hplt[1, 0] = hplt[iscan] - double(self.image_scale) * (Ny-1)/2.d
      wcs[*, iscan].hplt[0, 1] = hplt[iscan] + double(self.image_scale) * (Ny-1)/2.d
      wcs[*, iscan].hplt[1, 1] = hplt[iscan] + double(self.image_scale) * (Ny-1)/2.d        
    endfor
  endelse

  self -> fitscube_header_finalize, outhdr $
                       , no_checksum = no_checksum $
                       , coordinates = wcs $
                       , release_date = release_date $
                       , release_comment = release_comment $
                       , feature = feature $
                       , observer = observer $
                       , point_id = point_id $
                       , status = status
  if ~(status) then return
  
  ;; Add info about this step
  self -> headerinfo_addstep, outhdr $
                              , prstep = 'DATA-CURATION' $
                              , prpara = prpara $
                              , prproc = inam
     

  ;; Initialize the output file
  self -> fitscube_initialize, oname, outhdr, olun, fileassoc, dims, wcs = wcs  
  

  ;; Copy data frames
  Nframes = round(product(dims[2:*]))
  iframe = 0
  for iscan = 0, Nscans-1 do begin
    if filetype eq 'nb' and ~no_wcfile then begin
      ;; With old nb cubes we need to rely on integer wb cubes for padding
      ;; Because old padding might be not really same [median] value
      ;; after rotation and red__missing fails.
      red_fitscube_getframe,wcfile,wb_frame,iframe=iscan
      red_missing, wb_frame, /inplace, missing_value = 0
      pad_indx = where(wb_frame eq 0)
    endif
    for ituning = 0, Ntunings-1 do begin
      for istokes = 0, Nstokes-1 do begin

        red_progressbar, iframe, Nframes, /predict, 'Copying frames'
        
        red_fitscube_getframe, in_fileassoc, frame $
                               , iscan = iscan, ituning = ituning, istokes = istokes $
                               , fitscube_info = in_fitscube_info

        if filetype eq 'wb' or no_wcfile then begin
          red_missing, frame, image_out = frm, missing_value = 0
          pad_indx = where(frm eq 0)
        endif

        ;; Old wb cubes are integer and can't be padded with NaNs
        if bitpix eq 16 then frame = float(frame)
;        red_missing, frame, /inplace, missing_value = !Values.F_NaN
        frame[pad_indx] = !Values.F_NaN
        
        ;; Undo any mirroring done in the momfbd processing (already
        ;; included in -wcdirection?)
        if keyword_set(mirrorx) then frame = reverse(frame, 1, /over)        
        if keyword_set(mirrory) then frame = reverse(frame, 2, /over)        

        if wcdirection ne 0 then frame = rotate(temporary(frame), -wcdirection) ; Undo old
        frame = rotate(temporary(frame), direction)    ; Do new

        ;; Frames from the input cube are already rotated to a common
        ;; angle, we just need to adjust the angle to Solar N.
        frame = red_rotation(frame, full = ff, ang, 0, 0, background = !Values.F_NaN, nthreads = nthreads)
        
        red_fitscube_addframe, fileassoc, frame  $
                               , iscan = iscan, ituning = ituning, istokes = istokes $
                               , fitscube_info = fitscube_info
        
 
        iframe++

      endfor                    ; istokes
    endfor                      ; ituning   
  endfor                        ; iscan   

;;  free_lun, ilun
  red_fitscube_close, fileassoc, fitscube_info
  free_lun, olun
  if keyword_set(wcs_improve_spatial) then $
    red_fitscube_addwcs, oname, wcs $
                         , csyer_spatial_value = 5. $ ; 5 arcsec, minor rotation error may remain
                         , csyer_spatial_comment = 'Aligned with HMI images' $
                         , dimensions = fxpar(outhdr, 'NAXIS*') $
  else $
    red_fitscube_addwcs, oname, wcs $
                         , dimensions = fxpar(outhdr, 'NAXIS*') $
                         , csyer_spatial_value = 60. $ ; 1 arc minute
                         , csyer_spatial_comment = '[arcsec] Orientation known'

  if filetype eq 'nb' then begin
    for icmap = 0, n_elements(distortions)-1 do begin
      ;; Add cavity maps as WAVE distortions     
      cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)
      ;; cavity maps in older cubes were stored as (Nx, Ny, Nscans)
      ;; let's make reform in case we need to 'curate' newer cubes
      dd = reform(distortions[icmap].wave)
      for iscan = 0L, Nscans-1 do begin
        ;; Same operations as on data frames:        
        frame = dd[*, *, iscan]
        red_missing, frame, /inplace, missing_value = !Values.F_NaN
        if keyword_set(mirrorx) then frame = reverse(frame, 1, /over)        
        if keyword_set(mirrory) then frame = reverse(frame, 2, /over)        
        if wcdirection ne 0 then frame = rotate(temporary(frame), -wcdirection) ; Undo old
        frame = rotate(temporary(frame), direction)    ; Do new
        cavitymaps[*, *, 0, 0, iscan] = red_rotation(frame, full = ff, ang $
                                                     , 0, 0, background = !Values.F_NaN, nthreads = nthreads)
      endfor                      ; iscan
      red_fitscube_addcmap, oname, cavitymaps, cmap_number = icmap+1 $
                            , indx = red_expandrange(distortions[icmap].tun_index)
    endfor                       ; icmap
    ;; Copy all extensions, including varaible keywords. But exclude
    ;; WCS-TAB and WCSDVARR, they are WCS info already added in modified
    ;; form. Exclude also statistics keywords. They should be recalculated.
    red_fitscube_copyextensions, inname, oname $
                               , ext_list = ['WCS-TAB', 'WCSDVARR'], /ext_statistics, /ignore
  endif else $
    red_fitscube_copyextensions, inname, oname $
                               , ext_list = ['WCS-TAB'], /ext_statistics, /ignore
  

  red_fitscube_statistics, oname, /write
  
  if ~keyword_set(no_checksum) then $
    red_fitscube_checksums, oname
  
  if keyword_set(flip) then begin
    ;; Make a flipped version
    print, 'Flip it!'
    red_fitscube_flip, oname $
                       , flipfile = flipfile $
                       , overwrite = overwrite
  endif

  outname = oname
  print
  print, inam + ' : Input: ' + inname
  print, inam + ' : Output: ' + oname
  if keyword_set(flip) then print, inam + ' :    and  ' + flipfile
  
end




wdir = '/scratch/mats/convert_rotation/2018-09-30/CRISP/'
file = 'cubes_nb/nb_6173_2018-09-30T09:33:41_scans=0-19_stokes_corrected_im.fits'


case 2 of

  1 : begin
    wdir = '/scratch/olexa/2021-08-16/CRISP/'
    reffile = 'cubes_nb/nb_8542_2021-08-16T14:04:21_14:04:21=0-37_stokes_corrected_im.fits'
    file = 'cubes_nb/nb_8542_2021-08-16T14:04:21_14:04:21=0-4_corrected_im.fits'
    direction = 6
    rotation = -41.
  end

  2 : begin
    wdir = '/scratch/olexa/2021-08-16/CHROMIS/'
    reffile = 'cubes_nb/nb_3950_2021-08-16T14:02:51_14:02:51=0-19_corrected_im.fits'
    file = 'cubes_nb/nb_3950_2021-08-16T14:02:51_14:02:51=0-4_corrected_im.fits'
    direction = 1
    rotation = -42.
  end
  
endcase



cd, wdir

a =  crispred(/dev, /no)
a -> fitscube_nup, file, outname = outname, /over, direction = direction, rotation = rotation

red_fitscube_getwcs, outname, coordinates=wcs, distortions = distortions

fits_open, outname, fcb

red_fitscube_getframe, outname, outframe, iframe = 0
red_fitscube_getframe, reffile, refframe, iframe = 0

tighttv, refframe, 0
tighttv, outframe, 1


end
