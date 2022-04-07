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
; 
; :History:
; 
;    2022-03-29 : MGL. First version.
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
                       , release_comment = releasec $
                       , feature = feature $
                       , observer = observer $
                       , point_id = point_id $
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

  dir = self.out_dir
  
  if n_elements(outdir) eq 0 then outdir = self.out_dir+'cubes_converted/'
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
  
  ;; NB cube of WB cube?
  filetype = (strsplit(file_basename(inname),'_',/extract))[0]
  
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
      return
    endelse
  endif

  
  ;; Establish the orientation of the file
  if n_elements(mirrorx) eq 0 or n_elements(mirrory) eq 0 then begin

    cfgfiles = file_search(dir+'*mfbd*/'+timestamp+'/'+pref+'/cfg/*cfg', count = Ncfg)  

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

  if filetype eq 'nb' then begin
    indx = where(strmatch(inhdr, '*WCFILE*'), Nwhere)
    if Nwhere eq 0 then begin
      print, 'There is no information about wb cube in the header.'
      return
    endif else begin
      par = strmid(inhdr[indx],0,7)
      par_key = red_fitsgetkeyword(inhdr, par[0])
      wcfile = (json_parse(par_key))['WCFILE']
    endelse
    if ~file_test(wcfile) then begin
      print, inam + ' : Could not find the WB cube -- ', wcfile
      print,'Enter wb cube name manually.'
      wcfile = ''
      read,'WB cube name : ',wcfile
      if strtrim(wcfile,2) eq '' then return
      if file_dirname(wcfile) eq '.' then wcfile = 'cubes_wb/'+wcfile
      if ~file_test(wcfile) then begin
        print,'There is no such wb cube.'
        return
      endif
    endif
  endif else wcfile = inname

  ;; Read parameters from the WB cube
  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
  fxbreadm, bunit, row = 1 $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
            , wcANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01, wcdirection  
  fxbclose, bunit

  if n_elements(wcdirection) eq 0 and $
    n_elements(mirrorx) eq 0 and $
    n_elements(mirrory) eq 0 then begin

    print,'There is no information about previous 90 degrees rotations or mirroring of the cube.'
    ans=''
    read, 'Do you want to continue (Y/n): ',ans
    if strupcase(ans) eq 'N' then return
  endif
  if  n_elements(wcdirection) eq 0 then wcdirection = 0

  if keyword_set(do_wb_cube) or keyword_set(wcs_improve_spatial) then begin
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
      nm = file_basename(inname)
      prfx = strmid(nm,0,2)
      if prfx eq 'nb' then $
        wb_oname = file_dirname(inname) + '/wb' + strmid(inname,2,strlen(nm)) $
      else $
        wb_oname = file_dirname(inname) + '/wb_' + nm
    endelse

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
  
  ;; Angles to use are original derotation angles minus rotation
  ;; already done. This should be a scalar.
  ang = mean(ang - wcang) - (rotation-defrotation)*!dtor
  maxangle = abs(ang)
  ff = [maxangle, 0., 0., 0., 0., ang]    

  ;; Read one frame
  red_fitscube_getframe, inname, frame, iframe = 0
   
  ;; Get the frame with size implied by ff
  if wcdirection ne 0 then frame = rotate(temporary(frame), -wcdirection) ; Undo old direction
  if keyword_set(mirrorx) then frame = reverse(frame, 1, /over)        
  if keyword_set(mirrory) then frame = reverse(frame, 2, /over)
                                   
  frame = rotate(temporary(frame), direction)    ; Do the new direction
  frame = red_rotation(frame, full = ff, 0, 0, 0)

  ;; Set the spatial dimensions of the output file.
  dims[0:1] = size(frame, /dim)
  Nx_old = Nx
  Ny_old = Ny
  Nx = dims[0]  
  Ny = dims[1]  

  ;; Make header
  outhdr = inhdr

  if keyword_set(wcs_improve_spatial) then begin
    red_get_wcs,  wb_oname, coordinates = wb_wcs
    wcs.hpln = wb_wcs.hpln
    wcs.hplt = wb_wcs.hplt
  endif else begin
    hpln = median(wcs.hpln)
    hplt = median(wcs.hplt)
    ;; But what we want to tabulate is the pointing in the corners of
    ;; the FOV. Assume hpln and hplt are the coordinates of the center
    ;; of the FOV.
    wcs.hpln[0, 0, *, *] = hpln - double(self.image_scale) * (Nx-1)/2.d
    wcs.hpln[1, 0, *, *] = hpln + double(self.image_scale) * (Nx-1)/2.d
    wcs.hpln[0, 1, *, *] = hpln - double(self.image_scale) * (Nx-1)/2.d
    wcs.hpln[1, 1, *, *] = hpln + double(self.image_scale) * (Nx-1)/2.d

    wcs.hplt[0, 0, *, *] = hplt - double(self.image_scale) * (Ny-1)/2.d
    wcs.hplt[1, 0, *, *] = hplt - double(self.image_scale) * (Ny-1)/2.d
    wcs.hplt[0, 1, *, *] = hplt + double(self.image_scale) * (Ny-1)/2.d
    wcs.hplt[1, 1, *, *] = hplt + double(self.image_scale) * (Ny-1)/2.d
  endelse

  self -> fitscube_header_finalize, outhdr $
                       , no_checksum = no_checksum $
                       , coordinates = wcs $
                       , release_date = release_date $
                       , release_comment = releasec $
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
    for ituning = 0, Ntunings-1 do begin
      for istokes = 0, Nstokes-1 do begin

        red_progressbar, iframe, Nframes, /predict, 'Copying frames'
        
        red_fitscube_getframe, inname, frame $
                               , iscan = iscan, ituning = ituning, istokes = istokes

        ;;old wb cubes are integer and can't be padded with NaNs
        if size(frame,/type) eq 2 then frame = float(frame)
        red_missing, frame, /inplace, missing_value = !Values.F_NaN
        
        ;; Undo any mirroring done in the momfbd processing (already
        ;; included in -wcdirection?)
        if keyword_set(mirrorx) then frame = reverse(frame, 1, /over)        
        if keyword_set(mirrory) then frame = reverse(frame, 2, /over)        

        if wcdirection ne 0 then frame = rotate(temporary(frame), -wcdirection) ; Undo old
        frame = rotate(temporary(frame), direction)    ; Do new

        ;; Frames from the input cube are already rotated to a common
        ;; angle, we just need to adjust the angle to Solar N.
        frame = red_rotation(frame, full = ff, ang, 0, 0, background = !Values.F_NaN)
        
        self -> fitscube_addframe, fileassoc, frame  $
                                   , iscan = iscan, ituning = ituning, istokes = istokes
        
 
        iframe++

      endfor                    ; istokes
    endfor                      ; ituning   
  endfor                        ; iscan   
  

  self -> fitscube_finish, olun, wcs = wcs  

  if filetype eq 'nb' then begin
    for icmap = 0, n_elements(distortions)-1 do begin
      ;; Add cavity maps as WAVE distortions     
      cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)
      for iscan = 0L, Nscans-1 do begin
        ;; Same operations as on data frames:
        frame = rotate(distortions[icmap].wave[*, *, 0, 0, iscan], direction)
        red_missing, frame, /inplace, missing_value = !Values.F_NaN
        if keyword_set(mirrorx) then frame = reverse(frame, 1, /over)        
        if keyword_set(mirrory) then frame = reverse(frame, 2, /over)        
        frame = rotate(temporary(frame), -wcdirection) ; Undo old
        frame = rotate(temporary(frame), direction)    ; Do new
        cavitymaps[*, *, 0, 0, iscan] = red_rotation(frame, full = ff, ang $
                                                     , 0, 0, background = !Values.F_NaN)
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
