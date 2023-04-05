; docformat = 'rst'

;+
; Bypass MOMFBD restoration for bad quality data.
;
; Instead, prepare for making cubes by destretching and summing raw
; data. 
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
;   cfgfile : in, type=string
;
;      Path to a momfbd config file.
; 
; 
; :Keywords:
; 
;    clips : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment. 
;
;    overwrite : in, optional, type=boolean
;
;       Don't care if output is already on disk, overwrite with new
;       versions.
;
;    nostretch : in, optional, type=boolean
;
;       Skip destretching for these data. Might be useful if the data
;       quality is particularly low.
;
;    tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;    global_gain : in, optional, type=boolean
;
;       Set this keyword to not use time-dependent gains (for CRISP).
;
;    no_display : in, optional, type=boolean
;
;       Set this keyword not to show alignment process.
;
;    outdir_tag : in, optional, type=string
;
;       String to be add to putput directory name - 'results'+outdir_tag,
;       default - '_bypass'.
;
;    nthreads : in, optional, type=integer
;
;       Number of threads to be used in processing.
;
;    align_interactive : in, optional, type=boolean
;
;      Set this keyword to define the alignment FOV by use of the XROI GUI.
;
;    roi_align : in, optional, type=integer array
;
;      Coordinates of a region to be used for alignment.
;
; :History:
; 
;   2022-11-29 : MGL. First version.
; 
;   2023-01-11 : MGL. New keyword nostretch.
;
;   2023-04-03 : OA. New keywords global_gain, no_display, outdir_tag, nthreads,
;                    align_interactive, roi_align
; 
;-
pro red::bypass_momfbd, cfgfile $
                        , clips = clips $
                        , overwrite = overwrite $
                        , nostretch = nostretch $
                        , tiles = tiles $
                        , global_gain = global_gain $
                        , no_display = no_display $
                        , outdir_tag = outdir_tag $
                        , nthreads = nthreads $
                        , align_interactive = align_interactive $
                        , roi_align = roi_align

  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if ~keyword_set(outdir_tag) then outdir_tag = '_bypass'
  
  if n_elements(tiles) eq 0 then tiles =   [ 8, 16, 16, 32, 32]  
  if n_elements(clips) eq 0 then clips = 2*[16, 16,  8,  4,  2]  

  ;; Base default tiles and clips on num_points and max_local_shift?
  ;; Num_points gives final number of tiles. Max_local_shift gives
  ;; max(clips)? 

  red_make_prpara, prpara, cfgfile
  red_make_prpara, prpara, clips
  red_make_prpara, prpara, tiles  
  
  if n_elements(nthreads) eq 0 then nthreads = 10
  
  cfg = redux_readcfg(cfgfile)
  cfgdir = file_dirname(cfgfile)

  if file_test(cfgdir+'/fov_mask.fits') then fov_mask = readfits(cfgdir+'/fov_mask.fits')

  Nobjects = n_elements(redux_cfggetkeyword(cfg, 'OBJECT*'))

  trace = redux_cfggetkeyword(cfg, 'TRACE')

  image_nums = rdx_str2ints(redux_cfggetkeyword(cfg, 'IMAGE_NUMS'))
  Nexposures = long(n_elements(image_nums))

  ;; FOV from sim_x,sim_y and num_points.
  num_points = long(redux_cfggetkeyword(cfg, 'NUM_POINTS'))
  margin = num_points/8
  sim_xy = redux_cfggetkeyword(cfg, 'SIM_XY', count = cnt)
  if cnt gt 0 then begin
    sim_xy = rdx_str2ints(sim_xy)
    indx = indgen(n_elements(sim_xy)/2)*2
    indy = indx+1
    sim_x = sim_xy[indx]
    sim_y = sim_xy[indy]   
  endif else begin
    sim_x = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_X'))
    sim_y = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_Y'))
  endelse
  x0 = min(sim_x) + margin - num_points/2 
  x1 = max(sim_x) - margin + num_points/2 - 1
  y0 = min(sim_y) + margin - num_points/2 
  y1 = max(sim_y) - margin + num_points/2 - 1
  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1

  
  ;;                       ----- Anchor images -----  
  
  ;; Calculate stretch vectors for all raw WB (anchor) images. Also
  ;; make a sum of the stretched images. Save as the anchor object.
  ;; Disregard any diversity images, use only the focused WB images.

  anchor_output_file = redux_cfggetkeyword(cfg, 'OBJECT0.OUTPUT_FILE')
  anchor_gain_file = redux_cfggetkeyword(cfg, 'OBJECT0.CHANNEL0.GAIN_FILE')  
  anchor_dark_file = redux_cfggetkeyword(cfg, 'OBJECT0.CHANNEL0.DARK_TEMPLATE')  
  image_data_dir = redux_cfggetkeyword(cfg, 'OBJECT0.CHANNEL0.IMAGE_DATA_DIR')
  filename_template = redux_cfggetkeyword(cfg, 'OBJECT0.CHANNEL0.FILENAME_TEMPLATE')
  back_gain = redux_cfggetkeyword(cfg, 'OBJECT0.CHANNEL0.BACK_GAIN', count = DoBackscatter)
  if DoBackscatter then begin
    ;; Use the presence of the back_gain keyword only as a boolean.
    ;; Load the backscatter gain and psf the usual way.
    split_template = strsplit(filename_template, '_', /extract)
    detector = split_template[0]
    pref = split_template[2]
    self -> loadbackscatter, detector, pref, bgt, Psft
  endif

  filename_parts = strsplit(filename_template, '%', /extract)
  wfiles = file_search(image_data_dir+'/'+filename_parts[0]+'*', count = Nfiles)
  self -> extractstates, wfiles, wstates

  nframes_perfile = wstates.nframes

  wcube=red_readdata_multiframe(wfiles)
  dims = size(wcube, /dim)
  wcube1 = fltarr(dims)
  Nx = dims[0]
  Ny = dims[1]
  Nexposures = dims[2]
    
  anchor_dark = red_readdata(anchor_dark_file)
  anchor_gain = red_readdata(anchor_gain_file)
  mask = anchor_gain eq 0
  if n_elements(fov_mask) ne 0 then mask *= fov_mask

  for iexposure = 0, Nexposures-1 do begin
    red_progressbar, iexposure, Nexposures, 'Processing raw anchor images'
    im = wcube[*, *, iexposure]
    im -= anchor_dark
    if DoBackscatter then begin
;      im = rdx_descatter(im, bgt, Psft, verbose = verbose, nthreads = nthreads)
      im = red_cdescatter(im, bgt, Psft,nthreads = nthreads)
    endif
    im *= anchor_gain
    wcube[*, *, iexposure] = rdx_fillpix(im, nthreads=nthreads, mask = mask)
  endfor                        ; iexposure
  
  anchor_bg = median(wcube[*,*,0])

  xc = Nx/2
  yc = Ny/2
  xbd = round(20/float(self.image_scale))  < round(Nx*0.9)
  ybd = round(20/float(self.image_scale))  < round(Ny*0.9)
  align_size = [xbd, ybd]

  if n_elements(roi_align) ne 0 then begin
     xc = roi_align[0]
     yc = roi_align[1]
     align_size[0] = roi_align[2]
     align_size[1] = roi_align[3]
  endif else begin
     if keyword_set(align_interactive) then begin
        
        print
        print, 'Use the XROI GUI to either modify an initial alignment ROI/FOV or define a new one from scratch.'
        print, 'Select Quit in the File menu. The last ROI is used.'
        print

        ;; Define default roi
        X_in = xc + [-1,  1, 1, -1]*xbd/2 
        Y_in = yc + [-1, -1, 1,  1]*ybd/2
        roiobject_in = OBJ_NEW('IDLgrROI', X_in, Y_in)
        roiobject_in -> setproperty, name = 'Default'
        
        ;; Fire up the XROI GUI.
;    dispim = bytscl(red_histo_opt(total(wcube, 3)))
        dispim = bytscl(red_histo_opt(wcube[*,*,Nexposures/2]))
        xroi, dispim, regions_in = [roiobject_in], regions_out = roiobject, /block $
              , tools = ['Translate-Scale', 'Rectangle'] $
              , title = 'Modify or define alignment ROI'
        roiobject[-1] -> getproperty, roi_xrange = roi_x
        roiobject[-1] -> getproperty, roi_yrange = roi_y

        obj_destroy, roiobject_in
        obj_destroy, roiobject

        xc = round(mean(roi_x))
        yc = round(mean(roi_y))

        align_size = [round(roi_x[1])-round(roi_x[0]), round(roi_y[1])-round(roi_y[0])]

        roi_align = intarr(4)
        roi_align[0] = xc
        roi_align[1] = yc
        roi_align[2] = align_size[0]
        roi_align[3] = align_size[1]
     endif
  endelse

  np = 3
  shift = red_aligncube(wcube, np, xc = xc, yc = yc, xbd = align_size[0], ybd = align_size[1] $
                        , nthreads=nthreads, no_display = no_display)
  for iexposure = 0, Nexposures-1 do begin
    red_progressbar, iexposure, Nexposures, 'Shifting anchor images'
    wcube1[*, *, iexposure] $
       = red_rotation(wcube[*, *, iexposure], 0., shift[0,iexposure], shift[1,iexposure] $
                      , background = anchor_bg $
                      , nthreads=nthreads)
  endfor                        ; iexposure

  if ~keyword_set(nostretch) then begin
    for iexposure = 1, Nexposures-1 do begin
      red_progressbar, iexposure, Nexposures, 'Calculating stretch vectors for anchor images.'
;      grid = rdx_cdsgridnest(wcube1[*, *, iexposure-1], wcube1[*, *, iexposure] $
;                             , tiles, clips, nthreads=nthreads)
      grid = red_dsgridnest(wcube1[*, *, iexposure-1], wcube1[*, *, iexposure] $
                             , tiles, clips, nthreads=nthreads, mask=fov_mask)
      if iexposure eq 1 then begin
        gdims = size(grid, /dim)
        grids = dblarr(Nexposures, 2, gdims[1], gdims[2])
      endif
      grids[iexposure, *, *, *] = grid
    endfor                      ; iexposure
    
    for itile = 0, gdims[1]-1 do begin  
      red_progressbar, itile, gdims[1], 'Preparing displacement grids'
      for jtile = 0, gdims[2]-1 do begin
        ;; Note that we apply a cumulative displacement grid  
        grids[*, 0, itile, jtile] = total(grids[*, 0, itile, jtile], /cumulative)      
        grids[*, 1, itile, jtile] = total(grids[*, 1, itile, jtile], /cumulative)      
        ;; Subtract average displacements, don't want to stretch
        ;; to first frame
        grids[*, 0, itile, jtile] = grids[*, 0, itile, jtile] - mean(grids[*, 0, itile, jtile])
        grids[*, 1, itile, jtile] = grids[*, 1, itile, jtile] - mean(grids[*, 1, itile, jtile])      
      endfor
    endfor
  endif
  
  for iexposure = 0, Nexposures-1 do begin
    if keyword_set(nostretch) then begin
      ;; red_progressbar, iexposure, Nexposures, 'Shifting the anchor images.'
      ;; imm = red_rotation(wcube[*,*,iexposure] $
      ;;                    , 0.0, shift[0,iexposure], shift[1,iexposure] $
      ;;                    , background = anchor_bg $
      ;;                    , nthreads=nthreads, nearest=nearest)
      if n_elements(fov_mask) ne 0 then wcube1[*,*,iexposure] *= fov_mask
      ;;wcube1[*,*,iexposure] = imm
    endif else begin
      red_progressbar, iexposure, Nexposures, 'Shifting and stretching the anchor images.'
      imm = red_rotation(wcube[*,*,iexposure] $
                         , 0.0, shift[0,iexposure], shift[1,iexposure] $
                         , stretch_grid = reform(grids[iexposure,*,*,*]) $
                         , background = anchor_bg $
                         , nthreads=nthreads, nearest=nearest)
      if n_elements(fov_mask) ne 0 then imm *= fov_mask
      wcube1[*,*,iexposure] = imm      
    endelse
  endfor                        ; iexposure
  
  anchorim = mean(wcube1, dim = 3)
  red_missing, anchorim, /inplace, missing_type_wanted = 'nan'

;; Add the string "_bypass" to the "results" outdir, otherwise
;; generate the same filenames as momfbd would do (except with a
;; .fits extension rather than .momfbd). Add the .fitsheader file
;; from the regular output directory to the file header (but
;; remove/change some keywords that have to do with momfbd
;; processing).

  header = headfits(cfgdir+'/'+anchor_output_file + '.fitsheader')

  red_headerinfo_deletestep, header, /last
  if keyword_set(nostretch) then begin
    self -> headerinfo_addstep, header $
                                , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT' $
                                , prpara = prpara $
                                , prproc = inam
  endif else begin
    self -> headerinfo_addstep, header $
                                , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                                , prpara = prpara $
                                , prproc = inam
  endelse
  
  anchor_output_file = cfgdir+'/'+red_strreplace(anchor_output_file, 'results', 'results'+outdir_tag) + '.fits'
  
  if keyword_set(nostretch) then begin
    fxaddpar, header, 'FILENAME', file_basename(anchor_output_file), 'Aligned raw data'
  endif else begin
    fxaddpar, header, 'FILENAME', file_basename(anchor_output_file), 'Destretched raw data'
  endelse
  fxaddpar, header, 'CROP_X0', x0, 'Crop llx pixel'  
  fxaddpar, header, 'CROP_X1', x1, 'Crop urx pixel'  
  fxaddpar, header, 'CROP_Y0', y0, 'Crop lly pixel'  
  fxaddpar, header, 'CROP_Y1', y1, 'Crop ury pixel'  
  file_mkdir, file_dirname(anchor_output_file)
  red_writedata, anchor_output_file, anchorim[x0:x1,y0:y1], header = header, overwrite = overwrite 
                                ; , extra_info = { name:wfiles, Nx_orig:dims[0], Ny_orig:dims[1] $
                                ;                 , x0:x0, y0:y0, x1:x1, y1:y1 }



  ;;                       ----- All other images -----
  
  for iobject = 1, Nobjects-1 do begin
    
    ;; For each NB object, read the cube of raw NB images, and correct
    ;; for dark and flat (and possibly backscatter). Apply the align map
    ;; transformation, and stretch with the appropriate WB stretch
    ;; vectors. Sum the NB images and save as FITS files.

    ;; If TRACE, then do the same thing with the WB images that are
    ;; simultaneous to the NB images (except for the align map
    ;; transformation).

    object_string = 'OBJECT'+strtrim(iobject, 2)
    output_file = redux_cfggetkeyword(cfg, object_string+'.OUTPUT_FILE')
    gain_file = redux_cfggetkeyword(cfg, object_string+'.CHANNEL0.GAIN_FILE')  
    dark_file = redux_cfggetkeyword(cfg, object_string+'.CHANNEL0.DARK_TEMPLATE')  
    image_data_dir = redux_cfggetkeyword(cfg, object_string+'.CHANNEL0.IMAGE_DATA_DIR')
    filename_template = redux_cfggetkeyword(cfg, object_string+'.CHANNEL0.FILENAME_TEMPLATE')
    align_map = float(reform(strsplit(redux_cfggetkeyword(cfg, object_string+'.CHANNEL0.ALIGN_MAP') $
                                      , ',', /extract), 3, 3))
    back_gain = redux_cfggetkeyword(cfg, object_string+'.CHANNEL0.BACK_GAIN', count = DoBackscatter)
    if DoBackscatter then begin
      ;; Use the presence of the back_gain keyword only as a boolean.
      ;; Load the backscatter gain and psf the usual way.
      split_template = strsplit(filename_template, '_', /extract)
      if detector ne split_template[0] or pref ne split_template[2] then begin
        detector = split_template[0]
        pref = split_template[2]
        self -> loadbackscatter, detector, pref, bgt, Psft
      endif
    endif
    
    filename_parts = strsplit(filename_template, '%', /extract)
    files = file_search(image_data_dir+'/'+filename_parts[0]+'*'+strmid(filename_parts[1],3,strlen(filename_parts[1])-3), count = Nfiles)
    self -> extractstates, files, states

    if trace then toutput_file = red_strreplace(output_file, states[0].detector, wstates[0].detector)    
    
    ;; Find the wcube frame indices for this object
    match2, image_nums, states.framenumber, suba, file_indx
    undefine, indx
    for iindx = 0, n_elements(file_indx)-1 do begin
      if file_indx[iindx] eq 0 then begin
        istart = 0
      endif else begin
        istart = round(total(nframes_perfile[0:file_indx[iindx]-1]))
      endelse
      red_append, indx, lindgen(nframes_perfile[file_indx[iindx]]) + istart
    endfor                      ; iindx    

    if keyword_set(global_gain) then begin
      self -> get_calib, states[0] $
                         , darkdata = dark, darkstatus = darkstatus $
                         , cgaindata = gain, cgainstatus = gainstatus
      if darkstatus ne 0 then stop
      if gainstatus ne 0 then stop
    endif else begin
      dark = red_readdata(dark_file)
      gain = red_readdata(gain_file)
    endelse

    gain = rdx_img_transform(invert(align_map), gain, /preserve) >0
    mask = gain eq 0
    if n_elements(fov_mask) ne 0 then mask *= fov_mask
    
    cube = red_readdata_multiframe(files)
    dims = size(cube, /dim)
    cube1 = fltarr(dims)
    Nframes = dims[2]

    if Nframes ne n_elements(indx) then stop

    for iframe = 0, Nframes-1 do begin
      red_progressbar, iframe, Nframes, 'Processing raw images for object ' + strtrim(iobject, 2) $
                       + ' of ' + strtrim(Nobjects, 2)
      im = cube[*, *, iframe]
      im -= dark
      if DoBackscatter then begin
;        im = rdx_descatter(im, bgt, Psft, nthreads = nthreads)
        im = red_cdescatter(im, bgt, Psft, nthreads = nthreads)
      endif
      im = rdx_img_transform(invert(align_map), im, /preserve) >0
      im *= gain
      im = rdx_fillpix(im, nthreads=nthreads, mask = mask)
      if max(im) eq min(im) then stop

      cube[*, *, iframe] = im
    endfor                      ; iframe
    bg = median(cube[*,*,0])
   
    for iframe = 0, Nframes-1 do begin
      if keyword_set(nostretch) then begin
        red_progressbar, iframe, Nframes, 'Shifting the images.'
      endif else begin
        red_progressbar, iframe, Nframes, 'Shifting and stretching the images.'
      endelse
      iexposure = indx[iframe]
      xshift = shift[0,iexposure]
      yshift = shift[1,iexposure]
      ;;if ~keyword_set(nostretch) then grid = reform(grids[iexposure,*,*,*])
     
      if keyword_set(nostretch) then begin
         imm = red_rotation(cube[*,*,iframe] $
                            , 0.0, xshift, yshift $
                            , background = bg $
                            , nthreads=nthreads, nearest=nearest)
        if  n_elements(fov_mask) ne 0 then imm *= fov_mask
        cube1[*,*,iframe] = imm
      endif else begin
        imm = red_rotation(cube[*,*,iframe] $
                           , 0.0, xshift, yshift $
                           , stretch_grid = reform(grids[iexposure,*,*,*]) $
                           , background = bg $
                           , nthreads=nthreads, nearest=nearest)
         if n_elements(fov_mask) ne 0 then imm *= fov_mask
         cube1[*,*,iframe] = imm
      endelse
    endfor                      ; iexposure

    im = mean(cube1, dim = 3)
;    if Nfiles gt 1 then im = mean(cube1, dim = 3) else im = reform(cube1)
    red_missing, im, /inplace, missing_type_wanted = 'median'

    header = headfits(cfgdir+'/'+output_file + '.fitsheader')
    red_headerinfo_deletestep, header, /last

    if keyword_set(nostretch) then begin
      self -> headerinfo_addstep, header $
                                  , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT' $
                                  , prpara = prpara $
                                  , prproc = inam
      fxaddpar, header, 'FILENAME', file_basename(output_file), 'Aligned raw data'
    endif else begin
      self -> headerinfo_addstep, header $
                                  , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                                  , prpara = prpara $
                                  , prproc = inam
      fxaddpar, header, 'FILENAME', file_basename(output_file), 'Destretched raw data'
    endelse
    
    output_file = cfgdir+'/'+red_strreplace(output_file, 'results', 'results'+outdir_tag) + '.fits'
    fxaddpar, header, 'CROP_X0', x0, 'Crop llx pixel'  
    fxaddpar, header, 'CROP_X1', x1, 'Crop urx pixel'  
    fxaddpar, header, 'CROP_Y0', y0, 'Crop lly pixel'  
    fxaddpar, header, 'CROP_Y1', y1, 'Crop ury pixel'  
    red_writedata, output_file, im[x0:x1,y0:y1], header = header, overwrite = overwrite
    
    if trace then begin

      theader = headfits(cfgdir+'/'+toutput_file + '.fitsheader')
      toutput_file = cfgdir+'/'+red_strreplace(toutput_file, 'results', 'results'+outdir_tag) + '.fits'
      if file_test(toutput_file) and ~keyword_set(overwrite) then continue
      tcube1 = wcube1[*, *, indx]
      tim = mean(tcube1, dim = 3)
;      if n_elements(indx) gt 1 then tim = mean(tcube1, dim = 3) else tim = reform(tcube1)
      red_missing, tim, /inplace, missing_type_wanted = 'median' 
      red_headerinfo_deletestep, theader, /last
      if keyword_set(nostretch) then begin
        self -> headerinfo_addstep, theader $
                                    , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT' $
                                    , prpara = prpara $
                                    , prproc = inam
        fxaddpar, theader, 'FILENAME', file_basename(toutput_file), 'Aligned raw data'
      endif else begin
        self -> headerinfo_addstep, theader $
                                  , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                                  , prpara = prpara $
                                  , prproc = inam
        fxaddpar, theader, 'FILENAME', file_basename(toutput_file), 'Destretched raw data'
      endelse
      fxaddpar, theader, 'CROP_X0', x0, 'Crop llx pixel'  
      fxaddpar, theader, 'CROP_X1', x1, 'Crop urx pixel'  
      fxaddpar, theader, 'CROP_Y0', y0, 'Crop lly pixel'  
      fxaddpar, theader, 'CROP_Y1', y1, 'Crop ury pixel'  
      red_writedata, toutput_file, tim[x0:x1,y0:y1], header = theader, overwrite = overwrite
    endif
    
  endfor                        ; iobject

  
end 


;cd, '/scratch/mats/2016.09.19/CHROMIS-jan19'
;
;
;a = chromisred("config.txt", /dev, /no)
;root_dir = "/data/2016/2016.09/2016.09.19/"
;nthreads=20
;
;cfgfile = 'momfbd/09:28:36/3950/cfg/momfbd_reduc_3950_00004.cfg'
;a -> bypass_momfbd, cfgfile, /over
;cfgfile = 'momfbd/09:28:36/3950/cfg/momfbd_reduc_3950_00005.cfg'
;a -> bypass_momfbd, cfgfile, /over
;
;
;stop

cd, '/scratch/mats/2016.09.19/CRISP-aftersummer'

a = crispred(/dev, /no)

cfgfile = 'momfbd_nopd/09:30:20/6302/cfg/momfbd_reduc_6302_00001.cfg'
cfgfile = 'momfbd_nopd/09:30:20/6302/cfg/momfbd_reduc_6302_00002.cfg'

;cfgfile = 'momfbd_nopd/09:30:20/8542/cfg/momfbd_reduc_8542_00003.cfg'

a -> bypass_momfbd, cfgfile, /over

stop

a -> setproperty, 'filetype', 'FITS'
a -> make_scan_cube, 'momfbd_nopd/09:30:20/6302/cfg/results_bypass/', odir = 'cubes_bypass/', scanno = 2, /over

;a -> make_stokes_cubes, 'momfbd_nopd/09:30:20/6302/cfg/results_bypass/', 1, /redemodulate, /over

stop

a -> setproperty, 'filetype', 'FITS'
a -> make_wb_cube, 'momfbd_nopd/09:30:20/6302/cfg/results_bypass/', nametag = 'bypass', scanno = '1,2'
a -> make_nb_cube, 'cubes_wb/wb_6302_2016-09-19T09:30:20_09:30:20=1,2_bypass_corrected_im.fits'

end

