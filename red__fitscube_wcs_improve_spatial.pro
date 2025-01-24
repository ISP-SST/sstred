; docformat = 'rst'

;+
; Check orientation and position of the FOV of a fitscube file by
; comparison with an SDO/HMI intensity image collected close in time.
; 
; Use the iframe, ituning, istokes, and/or iscan keywords to select
; the fitscube frame to compare. For scan cubes, if the frame is not
; specified, will use the WB image in the extension.
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
;    filename : in, type=string
; 
;       The name of a WB, NB, or scan fitscube file.
; 
; 
; :Keywords:
; 
;    iframe : in, optional, type=integer, default="based on ituning, istokes, and iscan"
;
;       The frame index in the data cube seen as a 3D cube. 
;
;    iscan : in, optional, type=integer, default=0
;
;       The scan index, used to calculate iframe.   
; 
;    istokes  : in, optional, type=integer, default=0
;
;       The stokes index, used to calculate iframe.
;
;    ituning  : in, optional, type=integer, default=0
;
;       The tuning index, used to calculate iframe.
;
; :History:
; 
;   2020-03-23 : MGL. First version.
; 
;   2020-03-25 : MGL. New keyword log.
; 
;   2020-03-28 : MGL. Align with respect to two HMI images and
;                interpolate to get shifts at the right time.
; 
;   2020-03-28 : MGL. Rename method.
; 
;   2020-04-24 : MGL. Rewrite.
; 
;   2020-04-27 : MGL. Use WBIMAGE extension in scan cubes.
; 
;-
pro red::fitscube_wcs_improve_spatial, filename $
                                       , iframe = iframe $
                                       , iscan = iscan $
                                       , istokes = istokes $
                                       , ituning = ituning $
                                       , log = log

 
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  red_make_prpara, prpara, iframe
  red_make_prpara, prpara, iscan
  red_make_prpara, prpara, istokes
  red_make_prpara, prpara, ituning

  
  show_offset = 100

  wid_sst = 10
  wid_subim_sst = 11
  
  framespecified = max([n_elements(iframe)  gt 0 $
                        , n_elements(iscan)   gt 0 $
                        , n_elements(istokes) gt 0 $
                        , n_elements(ituning) gt 0 ])
  
  ;; Read the SST data
  
  h = headfits(filename)
  ;; The "direction" and "rotation" parameters are available in the
  ;; class object but you can override it when making the WB cube, so
  ;; better get them from the file header if possible.
  prprocs = fxpar(h,'PRPROC*')
  pos = where(strmatch(prprocs,'*make_*_cube*'), Nmatch)
  if Nmatch eq 0 then stop
  prpara_old = json_parse((fxpar(h,'PRPARA*'))[pos[0]])
  if prpara_old.haskey('DIRECTION') then begin
    ;; Direction in the header
    direction = prpara_old['DIRECTION']
  endif else begin
    ;; Direction from the class object, i.e., from the config file.
    direction = self.direction
    ;; Or maybe this should be 0?
  endelse
  if prpara_old.haskey('ROTATION') then begin
    ;; Rotation in the header
    rotation = prpara_old['ROTATION']
  endif else begin
    ;; Rotation from the class object, i.e., from the config file.
    rotation = self.rotation
    ;; Or maybe this should be 48?
  endelse
  ;; What would it mean if there aren't any ROTATION and
  ;; DIRECTION parameters in the PRPARA keyword? This should only
  ;; happen for old cubes, from before the parameters were introduced.
  ;; So we should stop and say that new cubes have to be made for this
  ;; method to work. 

  ;; Prefilter
  pref = strtrim(red_fitsgetkeyword(h, 'FILTER1'), 2)
  
  ;; Read the data frame. For scan cubes, use the WB image in the
  ;; extension unless a frame is specified.
  if ~framespecified && max(strmatch(prprocs,'*make_scan_cube')) then begin ; Scan cube. 
    im_sst = MRDFITS(filename, 'WBIMAGE', STATUS=status)
    if status lt 0 then return
    iframe = 0
  endif else status = -100
  if status lt 0 then begin     ; Read the specified frame
    if ~framespecified then begin
      ;; If frame is not specified then take a middle frame from a cube.
      scn = red_fitsgetkeyword(filename, 'SCANNUM', comment = comment, variable_values = scannum)
      scans = reform(scannum.values)
      in = n_elements(scans)/2
      iscan = in ;scans[in]-scans[0]
      istokes = 0
      ituning = 0
    endif
    red_fitscube_getframe, filename, im_sst $ 
                           , iframe = iframe $
                           , iscan = iscan $
                           , istokes = istokes $
                           , ituning = ituning
  endif
  
  im_sst = im_sst/median(im_sst)
  dims_sst = size(im_sst, /dim)
  date_beg = fxpar(h, 'DATE-BEG')
  time_beg = red_time2double((strsplit(date_beg, 'T', /extract))[1])

  sst_image_scale = float(self -> imagescale(pref, /use_config))
  
  ;; SST FOV WCS info, array coordinates in arc sec
  red_fitscube_getwcs, filename, coordinates=wcs

  if finite(wcs[iframe].hpln[0]) and finite(wcs[iframe].hplt[0]) then $
     sst_spatial_ok = 1 else sst_spatial_ok = 0
  
  time_sst = wcs[iframe].time[0]
  red_show, im_sst, w = wid_sst, title = 'SST @ '+red_timestring(time_sst, n = 0) $
            , /scroll, offs = show_offset
  

  ;; Due to solar rotation, the features we are aligning with can move
  ;; significantly between HMI exposures. Read the two HMI images
  ;; collected before and after the SST image, resp. Align with both
  ;; and interpolate to get the shifts at the right time.

  ;; Download HMI images - one available every 15 min.
  timestamps_hmi = red_timestring([floor(time_sst/(60*60/4)) $
                                   , ceil(time_sst/(60*60/4))] * (60*60/4) $
                                  , Nsecdecimals=0)
  htimes = red_time2double(timestamps_hmi)
  timestamps_hmi = [red_strreplace(timestamps_hmi, ':', '', n = 2)]
  red_download, timestamps_hmi = timestamps_hmi

  ;; The HMI images closest in time
  hnames = 'downloads/HMI/' $
           + red_strreplace(self.isodate, '-', '', n = 2) $
           + '_' + timestamps_hmi + '_Ic_flat_4k.jpg'

  if file_test(hnames[0]) and file_test(hnames[1]) then begin
    ;; Read the HMI image before the SST image was collected  
    read_jpeg, hnames[0], im_hmi_before
    im_hmi_before = reform(im_hmi_before[1, *, *]) ; use green channel
    im_hmi_before = im_hmi_before/median(red_centerpic(im_hmi_before,sz=1000))

    ;; Read the HMI image after the SST image was collected
    read_jpeg, hnames[1], im_hmi_after
    im_hmi_after = reform(im_hmi_after[1, *, *]) ; use green channel
    im_hmi_after = im_hmi_after/median(red_centerpic(im_hmi_after,sz=1000))

    dims_hmi = size(im_hmi_before, /dim)
    hmi_image_scale = 0.5         ; arc sec / pixel
  endif else begin
    print,'Failed to get HMI images.'
    return
  endelse
  
  ;; Display the "before" HMI image
  red_show, im_hmi_before, /scroll, offs = show_offset $
            , w = wid_hmi_before $
            , title = 'HMI @ '+red_timestring(htimes[0], n = 0)

  ;; Draw axes
  red_draw_lines, dims_hmi[0]/2, dims_hmi[1]/2, color='yellow', /device, /axes

  ;; Draw FOV box
  if sst_spatial_ok then $
     red_draw_lines, color='yellow', /device, /minmax $
                     , wcs[iframe].hpln/hmi_image_scale + dims_hmi[0]/2 $
                     , wcs[iframe].hplt/hmi_image_scale + dims_hmi[1]/2
  
  print
  print, inam + ' : Please use scroll bars (if present) to'
  print, inam + '   1. center the SST window on a structure with recognicable orientation.'
  print, inam + '   2. center the HMI window on the SST FOV. (If there is a yellow'
  print, inam + '      rectangle, it should be near it.)'
  print, inam + '   We will worry about the alignment later but first:'
  print
  s = ''
  read, 'Do the orientations (approximate 90 deg rotations and mirroring) match in the two images [yN]? ', s
  if strupcase(strmid(s, 0, 1)) ne 'Y' then begin

    ;; First invert the current "direction" of the SST data, then
    ;; determine the "direction" needed to get the matching
    ;; orientaiton.
    dims_sst_small = round(dims_sst*sst_image_scale/hmi_image_scale)
    im_sst_orig = congrid(red_rotate(im_sst, -direction) $
                          , dims_sst_small[0], dims_sst_small[1])
    maxdim_small = max(dims_sst_small)
    im_sst_orig = red_centerpic(im_sst_orig, sz = maxdim_small, z = median(im_sst_orig))

    window,9, xs = maxdim_small*8, ys = maxdim_small
    for i = 0, 7 do begin
      tvscl, red_rotate(im_sst_orig, i), i
      cgtext, i*maxdim_small+10, 10, strtrim(i, 2), /device, color = 'yellow'
    endfor                      ; i
    for i = 1, 7 do begin
      red_draw_lines, i*maxdim_small*[1, 1], maxdim_small*[0, 1], /device, color = 'yellow'
    endfor                      ; i

    ;; Select
    print
    print, 'Please click on the panel that is oriented like the HMI image.'
    cursor, x_select, y_select, /device
    direction_select = x_select/maxdim_small
    
    if direction_select eq direction then begin
      print, inam + ' : This is the orientation of the cube, no change.'
    endif else begin
      print, inam + ' : Please edit config.txt, change or add the following line:'
      print, inam + '   direction='+strtrim(direction_select, 2)
      print, inam + '   Then make a new WB cube in a new IDL session and continue from that.'
      wdelete, 9, wid_sst
      return
    endelse

    wdelete, 9    
  endif 

  print, inam + ' : Will now try to measure the misalignment with a subfield.'

  ;; Use XROI GUI to select a rectangular area. 
  
  ;; Initialize the FOV
  X_in = [0, dims_sst[0], dims_sst[0], 0] + [1,-1,-1, 1]*10
  Y_in = [0, 0, dims_sst[1], dims_sst[1]] + [1, 1,-1,-1]*10
  roi_in = OBJ_NEW('IDLgrROI', X_in, Y_in)
  roi_in -> setproperty, name = 'Default'

  print
  print, inam + ' : Use the XROI GUI to either modify the initial ROI/FOV or define a'
  print, inam + '   new one from scratch. Select a FOV with at least one spot or pore,'
  print, inam + '   suitable for cross-correlation.'
  print
  print, inam + '   Then select Quit in the File menu. The last ROI is used.'
  print
  
  ;; Fire up the XROI GUI.
  xroi, bytscl(im_sst), regions_in = [roi_in], regions_out = roi, /block $
        , tools = ['Translate-Scale', 'Rectangle'] $
        , title = 'Modify or define FOV for alignment'

  roi[-1] -> getproperty, roi_xrange = roi_x
  roi[-1] -> getproperty, roi_yrange = roi_y

  obj_destroy, roi_in
  obj_destroy, roi

  ;; We don't want out-of-bounds coordinates here
  x0 = round(roi_x[0]) >0
  y0 = round(roi_y[0]) >0
  x1 = round(roi_x[1]) <(dims_sst[0]-1)
  y1 = round(roi_y[1]) <(dims_sst[1]-1)

  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1

  ;; Draw the selected box in the SST image window
  wshow, wid_sst
  wset, wid_sst
  red_draw_lines, [x0,x1], [y0,y1], color='blue', /device, /minmax

  ;; Show the subimage in its own window
  subim_sst = im_sst[x0:x1, y0:y1]  
  red_show,subim_sst, w = wid_subim_sst $
           , title = 'SST subimage @ '+red_timestring(time_sst, n = 0)

  ;; We define the point (xc,yc) as the center of the SST image,
  ;; (xs,ys) as the center of the subimage selected with the XROI
  ;; interface. The point we mouse-click on later is denoted as
  ;; (xp,yp). Various versions of these coordinates are distinguished
  ;; by tags following an underscore. 

  if sst_spatial_ok then begin
    ;; Nominal (defined by wcs info in file) SST FOV center
    ;; coordinates in arc sec.
    xc_nominal = mean(wcs[iframe].hpln)
    yc_nominal = mean(wcs[iframe].hplt)

    xc_dif = mean(wcs[-1].hpln) - xc_nominal
    yc_dif = mean(wcs[-1].hplt) - yc_nominal

    ;; Nominal subimage FOV center coordinates in arc sec. Center of
    ;; original FOV plus offset from center of the XROI-selected
    ;; subfield.
    xs_nominal = xc_nominal + ((x1+x0)/2. - dims_sst[0]/2.) * sst_image_scale
    ys_nominal = yc_nominal + ((y1+y0)/2. - dims_sst[1]/2.) * sst_image_scale

    ;; Where is this in HMI image pixel cordinates?
    xc_nominal_hmipix = xs_nominal / hmi_image_scale - dims_hmi[0]/2.
    yc_nominal_hmipix = ys_nominal / hmi_image_scale - dims_hmi[1]/2.
  endif
  
  ;; Now we identify a point that can berecognized in both images. The
  ;; mouse-click in the SST image represents the nominal position of
  ;; that point, while the mouse-click in the HMI image represents the
  ;; actual position. (To be refined later with cross correlation.)
  
  print
  print, inam + ' : Please click on a position in the SST subimage, that you can identify in the HMI image.'
  wshow, wid_subim_sst
  wset, wid_subim_sst
  cursor, xp_click, yp_click, /device

  if sst_spatial_ok then begin
    xp_nominal = xs_nominal + (xp_click-Nx/2.) * sst_image_scale
    yp_nominal = ys_nominal + (yp_click-Ny/2.) * sst_image_scale
  endif
  
  print
  print, inam + ' : Click on the same position in the HMI image.'
  wshow, wid_hmi_before
  wset, wid_hmi_before
  cursor, xp_hmipix, yp_hmipix, /device

  xp = (xp_hmipix - dims_hmi[0]/2.) * hmi_image_scale
  yp = (yp_hmipix - dims_hmi[1]/2.) * hmi_image_scale

  if sst_spatial_ok then begin
    ;; Shift in arc sec defined by mouse clicks
    xp_shift = xp_nominal - xp
    yp_shift = yp_nominal - yp

    print
    print, 'Clicked point is shifted by (' $
         + strtrim(round(xp_shift), 2) + ',' $
           + strtrim(round(yp_shift), 2) + ') arcsec.' 
    print, 'Dimensions of SST image are (' $
           + strtrim(round(dims_sst[0]*sst_image_scale), 2) + ',' $
           + strtrim(round(dims_sst[1]*sst_image_scale), 2) + ') arcsec.'
    print
  endif
  
  ;; Now read out a subim from im_hmi_before that (approximately)
  ;; matches the sst image.

  ;; The clicked SST image point (xp_sstpix, yp_sstpix) is then "the
  ;; same" as clicked (xp_hmipix, yp_hmipix) point in the HMI image.
  ;; Each SST pixel corresponds to self.image_scale/hmi_image_scale
  ;; HMI pixels. So the center point (Nx/2,Ny/2) of the SST image is
  ;; (xp_sstpix, yp_sstpix) - (Nx/2,Ny/2) away from the clicked
  ;; point. That should correspond to (xp_hmipix, yp_hmipix) -
  ;; ((xp_sstpix, yp_sstpix) - (Nx/2,Ny/2)) * self.image_/hmi_image_scale in the HMI image.
;  xc_hmipix = round(xp_hmipix - (xp_sstpix-Nx/2.)*sst_image_scale/hmi_image_scale)
;  yc_hmipix = round(yp_hmipix - (yp_sstpix-Ny/2.)*sst_image_scale/hmi_image_scale)
;  ;; The size of the HMI subimage should match the SST subimage. 
;  xsz_hmipix = round(Nx * sst_image_scale/hmi_image_scale)
;  ysz_hmipix = round(Ny * sst_image_scale/hmi_image_scale)
;  subim_hmi_before = congrid(red_pic_at_coord(im_hmi_before $
;                                              , round(xc_hmipix), round(yc_hmipix) $
;                                              , xsz_hmipix, ysz_hmipix) $
;                             , Nx, Ny, cub=-0.5)

  ;; Pixel coordinates of the FOV of subim_sst in the HMI image (all
  ;; pixels, for interpolation):
  xx_hmipix = xp_hmipix + (findgen(Nx) - xp_click) * sst_image_scale/hmi_image_scale
  yy_hmipix = yp_hmipix + (findgen(Ny) - yp_click) * sst_image_scale/hmi_image_scale

  ;; Center of alignment subimage
  xs_hmipix = mean(xx_hmipix)
  ys_hmipix = mean(yy_hmipix)
  

  ;; Draw alignment subimage, with and without shift
  wshow, wid_hmi_before
  wset, wid_hmi_before
  if sst_spatial_ok then begin
    red_draw_lines, /minmax, /device, color='yellow' $ ; Relative to nominal position
                    , xx_hmipix + xp_shift/hmi_image_scale $
                    , yy_hmipix + yp_shift/hmi_image_scale 
    print
    print, 'The position of the SST subimage is marked in the HMI image window.'
    print, 'In yellow based on the nominal pointing information and'
    print, 'in green according to the HMI image mouse click.'
    print
  endif
  red_draw_lines, /minmax, /device, color='green' $ ; Relative to clicked position
                  , xx_hmipix $
                  , yy_hmipix

   
  ;; Crop HMI images to that FOV
  subim_hmi_before = interpolate(im_hmi_before, xx_hmipix, yy_hmipix, cub=-0.5, /grid)
  subim_hmi_after  = interpolate(im_hmi_after,  xx_hmipix, yy_hmipix, cub=-0.5, /grid)

  print
  print, 'I will now blink two images: 0) the SST subimage and 1) the corresponding'
  print, 'subimage of the HMI image magnified to the image scale of the SST image.'
  print, 'Their alignment will be estimated by use of a model fit and the result will'
  print, 'also be blinked (at that point with the SST image downgraded to the resolution'
  print, 'of the HMI image). Then this will be repeated for one more HMI image'
  print, '(collected after the SST image). Please strike RETURN when asked to.'
  print
  
  red_show, subim_hmi_before, w = 1
  blink, [wid_subim_sst,1]
  wdelete, 1, wid_subim_sst, wid_sst

  ;; We need to substitute NaNs 
  indx = where(~finite(subim_sst))
  subim_sst[indx] = median(subim_sst)
  ;; Estimate rotation with mpfit
  p_before = red_rot_magn_align(subim_hmi_before, subim_sst, /displ, maxiter = 10)
  p_after  = red_rot_magn_align(subim_hmi_after,  subim_sst, /displ, maxiter = 10)

;  p_between = red_rot_magn_align(subim_hmi_before, subim_hmi_after, /displ)
  
  print
  print, 'Config file rotation [deg]:   ', self.rotation
  print, 'Original cube rotation [deg]: ', rotation
  print, 'Measured rotational misalignment [deg]: ', p_before[0]
  print, 'New measured rotation parameter [deg]:  ', rotation - p_before[0]
  print
  read, 'Do you want to use the measured rotation and make a new cube [yN]? ', s
  if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
    print
    print, 'Change the rotation parameter in config.txt to the measured value.'
    print, 'Then make a new cube in a new IDL session and rerun fitscube_wcs_improve_spatial'
    print, 'with that (or just use the rotation keyword if you want to iterate a few times).'
    print
    wdelete, 1, wid_subim_sst, wid_sst
    return
  endif

  ;; If we are past this point, we can assume the rotation angle is fine!

  ;; First calculate the "before" and "after shifts
  shifts_subim_before = -p_before[2:3] * sst_image_scale
  shifts_subim_after  = -p_after[2:3]  * sst_image_scale

  ;; Draw "before" alignment subimage with final coordinates
  wshow, wid_hmi_before
  wset, wid_hmi_before
  red_draw_lines, /minmax, /device, color='darkgreen' $
                  , xx_hmipix - shifts_subim_before[0]/hmi_image_scale $
                  , yy_hmipix - shifts_subim_before[1]/hmi_image_scale
;  red_draw_lines, /minmax, /device, color='orange' $
;                  , xx_hmipix - shifts_before[0]/hmi_image_scale $
;                  , yy_hmipix - shifts_before[1]/hmi_image_scale 

  print
  print, 'The position of the SST subimage is now marked in the HMI image window'
  print, 'in dark green. This time according to the HMI image mouse click and subsequent'
  print, 'alignment fitting. Compare to the blue box in the full size SST image.'
  print
  
  if sst_spatial_ok then begin
    shifts_before = [xp_shift, yp_shift] + shifts_subim_before
    shifts_after  = [xp_shift, yp_shift] + shifts_subim_after

    ;; Draw shifted FOV box
    wshow, wid_hmi_before
    wset, wid_hmi_before
    red_draw_lines, color='green', /device, /minmax $
                  , (wcs[iframe].hpln - shifts_before[0])/hmi_image_scale + dims_hmi[0]/2 $
                  , (wcs[iframe].hplt - shifts_before[1])/hmi_image_scale + dims_hmi[1]/2
  endif else begin

    ;; Draw the box based on absolute numbers rather than shift!

    ;; Center of the subfield in SST pixel coordinates:
    xs_sstpix = (x0+x1)/2.
    ys_sstpix = (y0+y1)/2.

    ;; Center of the SST image in SST pixel coordinates:
    xc_sstpix = dims_sst[0]/2.
    yc_sstpix = dims_sst[1]/2.
    
    ;; Center of the SST image in HMI pixel coordinates:

    xc_hmipix = xs_hmipix + (xc_sstpix-xs_sstpix) * sst_image_scale/hmi_image_scale
    yc_hmipix = ys_hmipix + (yc_sstpix-ys_sstpix) * sst_image_scale/hmi_image_scale

    xc = (xc_hmipix - dims_hmi[0]/2.) * hmi_image_scale
    yc = (yc_hmipix - dims_hmi[1]/2.) * hmi_image_scale

    hpln = xc + [[-1, 1], [-1, 1]] * dims_sst[0]/2. * sst_image_scale
    hplt = yc + [[-1, -1], [1, 1]] * dims_sst[1]/2. * sst_image_scale

    ;; Draw FOV box
    wshow, wid_hmi_before
    wset, wid_hmi_before
    red_draw_lines, color='green', /device, /minmax $
                    , hpln/hmi_image_scale + dims_hmi[0]/2 $
                    , hplt/hmi_image_scale + dims_hmi[1]/2
  
  endelse
  


  print
  print, 'The estimated position of the full SST image is now marked in the HMI image window'
  print, 'in green. Compare to the full size SST image.'
  print
  
;  wshow, wid_hmi_after
;  wset, wid_hmi_after
;  red_draw_lines, color='green', /device, /minmax $
;                  , wcs[iframe].hpln/hmi_image_scale + dims_hmi[0]/2 $
;                  , wcs[iframe].hplt/hmi_image_scale + dims_hmi[1]/2

  if sst_spatial_ok then begin

    ;; Interpolate to get shifts at the time of the SST image
    shifts_subim = ( shifts_subim_before*(htimes[1]-time_sst) $
                     + shifts_subim_after*(time_sst-htimes[0]) ) / (htimes[1]-htimes[0])
    shifts = [xp_shift, yp_shift] + shifts_subim
    
    print, inam + ' : Detected X displacement: ' + strtrim(shifts[0], 2) + ' arc sec.'
    print, inam + ' : Detected Y displacement: ' + strtrim(shifts[1], 2) + ' arc sec.'

    ;; Modify the WCS info (for all frames)
    wcs.hpln -= shifts[0]
    wcs.hplt -= shifts[1]

    if keyword_set(log) then begin
      lfile = 'wcs_improve_spatial.log'
;    if file_test(lfile) then openu, llun, lfile, /get_lun else openw, llun, lfile, /get_lun
;    printf, llun, date_beg+' '+strtrim(x_shift, 2)+' '+strtrim(y_shift, 2)
;    free_lun, llun
      ;; For some reason the log file kept getting overwritten in spite
      ;; of being opened with openu. So we'll do it differently for now.
      spawn, 'echo '+ date_beg+' '+strtrim(shifts[0], 2)+' '+strtrim(shifts[1], 2) + ' >> ' + lfile
    endif
  endif else begin

    ;; Set the WCS info (for all frames)
    for iwcs = 0, n_elements(wcs)-1 do begin
      wcs[iwcs].hpln = hpln
      wcs[iwcs].hplt = hplt
    end
  endelse

  read, 'Do you want to update the WCS info in the cube file [yN]? ', s

  if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
    
    ;; Write the updated WCS info back into the file. 
    red_fitscube_addwcs, filename, wcs $
                         , csyer_spatial_value = 5. $ ; 5 arcsec, minor rotation error may remain
                         , csyer_spatial_comment = 'Aligned with HMI images' $
                         , dimensions = fxpar(h, 'NAXIS*') $
                         , /update
  
    hdr = headfits(filename)
    self -> headerinfo_addstep, hdr $
                                , prstep = 'ALIGNMENT-SPATIAL' $
                                , prpara = prpara $
                                , prproc = inam $
                                , prref = file_basename(hnames)
    red_fitscube_newheader, filename, hdr
  endif

end

if 1 then begin
  cd, '/scratch/mats/2016.09.19/CHROMIS-jan19'
  a = chromisred(/dev)
  a -> fitscube_wcs_improve_spatial, 'cubes_wb/wb_3950_2016-09-19T09:28:36_scans=0-3_test_corrected_im.fits'
endif else begin
  cd, '/scratch/mats/2016.09.19/CRISP-aftersummer'
  a = crispred(/dev)
  a -> fitscube_wcs_improve_spatial, 'cubes_wb/wb_6302_2016-09-19T09:30:20_scans=0-2_corrected_im.fits'
end


end
