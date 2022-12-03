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
;    tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
; 
; :History:
; 
;   2022-11-29 : MGL. First version.
; 
;-
pro red::bypass_momfbd, cfgfile $
                        , clips = clips $
                        , overwrite = overwrite $
                        , tiles = tiles

                        
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  outdir_tag = '_bypass'
  
  if n_elements(tiles) eq 0 then tiles =   [ 8, 16, 16, 32, 32]  ;, 64, 64]
  if n_elements(clips) eq 0 then clips = 2*[16, 16,  8,  4,  2]  ;,  2,  1  ]

  ;; Base default tiles and clips on num_points and max_local_shift?
  ;; Num_points gives final number of tiles. Max_local_shift gives
  ;; max(clips)? 

  red_make_prpara, prpara, cfgfile
  red_make_prpara, prpara, clips
  red_make_prpara, prpara, tiles  
  
  if n_elements(nthreads) eq 0 then nthreads = 10
  
  cfg = redux_readcfg(cfgfile)
  cfgdir = file_dirname(cfgfile)

  Nobjects = n_elements(redux_cfggetkeyword(cfg, 'OBJECT*'))

  trace = redux_cfggetkeyword(cfg, 'TRACE')

  image_nums = rdx_str2ints(redux_cfggetkeyword(cfg, 'IMAGE_NUMS'))
  Nexposures = long(n_elements(image_nums))

  ;; FOV from sim_x,sim_y and num_points.
  num_points = long(redux_cfggetkeyword(cfg, 'NUM_POINTS'))
  margin = num_points/8
  sim_xy = redux_cfggetkeyword(cfg, 'SIM_XY', count = cnt)
  if cnt gt 0 then begin
    sim_xy = reform(rdx_str2ints(sim_xy), 2, cnt/2)
    sim_x = sim_xy[0, *]    
    sim_y = sim_xy[1, *]    
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
  
  anchor_dark = red_readdata(anchor_dark_file)
  anchor_gain = red_readdata(anchor_gain_file)
  dims = size(anchor_dark, /dim)
  wcube = fltarr([dims, Nexposures])  
  wcube1 = fltarr([dims, Nexposures])  
  for iexposure = 0, Nexposures-1 do begin
    red_progressbar, iexposure, Nexposures, 'Reading raw anchor images'
    file_num = image_nums[iexposure]
    filename = filename_parts[0] + string(file_num, format = '(i07)')
    if filename_parts[1] ne '07d' then filename += '.fits' ; Add .fits extension
    im = red_readdata(image_data_dir+'/'+filename)
    im -= anchor_dark
    if DoBackscatter then begin
      im = rdx_descatter(im, bgt, Psft, verbose = verbose, nthreads = nthreads)
    endif
    im *= anchor_gain
    mask = anchor_gain eq 0
    wcube[*, *, iexposure] = rdx_fillpix(im, nthreads=nthreads, mask = mask)
  endfor                        ; iexposure
  
  anchor_bg = median(wcube)

  np = 3
  shift = red_aligncube(wcube, np, xc = Nx/2, yc = Ny/2, xbd = 500, ybd = 500, nthreads=nthreads)
  for iexposure = 0, Nexposures-1 do begin
    red_progressbar, iexposure, Nexposures, 'Shifting anchor images'
    wcube1[*, *, iexposure] $
       = red_rotation(wcube[*, *, iexposure], 0., shift[0,iexposure], shift[1,iexposure] $
                      , background = anchor_bg, nthreads=nthreads)
  endfor                        ; iexposure
  
  for iexposure = 1, Nexposures-1 do begin
    red_progressbar, iexposure, Nexposures, 'Calculating stretch vectors for anchor images.'
    grid = rdx_cdsgridnest(wcube1[*, *, iexposure-1], wcube1[*, *, iexposure] $
                           , tiles, clips, nthreads=nthreads)
    if iexposure eq 1 then begin
      gdims = size(grid, /dim)
      grids = dblarr(Nexposures, 2, gdims[1], gdims[2])
    endif
    grids[iexposure, *, *, *] = grid
  endfor                        ; iexposure

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
  
  tighttv,[wcube[*,*,0],wcube1[*,*,0]]
  for iexposure = 0, Nexposures-1 do begin
    red_progressbar, iexposure, Nexposures, 'Stretching the anchor images.'
    wcube1[*,*,iexposure] = red_rotation(wcube[*,*,iexposure] $
                                         , 0.0, shift[0,iexposure], shift[1,iexposure] $
                                         , stretch_grid = reform(grids[iexposure,*,*,*]) $
                                         , background = anchor_bg $
                                         , nthreads=nthreads, nearest=nearest)
  endfor                        ; iexposure
  anchorim = mean(wcube1, dim = 3)
  

;;for iexposure = 0, Nexposures-1 do tvscl,[wcube[*,*,iexposure],wcube1[*,*,iexposure]]


;; Add the string "_bypass" to the "results" outdir, otherwise
;; generate the same filenames as momfbd would do (except with a
;; .fits extension rather than .momfbd). Add the .fitsheader file
;; from the regular output directory to the file header (but
;; remove/change some keywords that have to do with momfbd
;; processing).

  header = headfits(cfgdir+'/'+anchor_output_file + '.fitsheader')

  red_headerinfo_deletestep, header, /last
  self -> headerinfo_addstep, header $
                              , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                              , prpara = prpara $
                              , prproc = inam

  anchor_output_file = cfgdir+'/'+red_strreplace(anchor_output_file, 'results', 'results'+outdir_tag) + '.fits'
;  anchor_output_file = cfgdir+'/'+red_strreplace(anchor_output_file, 'results', 'results'+outdir_tag) + '.momfbd'
  
  fxaddpar, header, 'FILENAME', file_basename(anchor_output_file), 'Destretched raw data'
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
    files = file_search(image_data_dir+'/'+filename_parts[0]+'*', count = Nfiles)
    self -> extractstates, files, states

    if trace then toutput_file = red_strreplace(output_file, states[0].detector, wstates[0].detector)    
    

    match2, image_nums, states.framenumber, suba, indx    

    dark = red_readdata(dark_file)
    gain = red_readdata(gain_file)
    
    cube = fltarr([dims, Nfiles])  
    cube1 = fltarr([dims, Nfiles])
    
    for ifile = 0, Nfiles-1 do begin
      red_progressbar, ifile, Nfiles, 'Reading raw images for object ' + strtrim(iobject, 2) $
                       + ' of ' + strtrim(Nobjects, 2)
;      file_num = image_nums[ifile] ; <-------------------------
;      filename = filename_parts[0] + string(file_num, format = '(i07)')
;      if filename_parts[1] ne '07d' then filename += '.fits' ; Add .fits extension
      im = red_readdata(files[ifile])
      im -= dark
      if DoBackscatter then begin
        im = rdx_descatter(im, bgt, Psft, verbose = verbose, nthreads = nthreads)
      endif
      im *= gain
      mask = gain eq 0
      im = rdx_fillpix(im, nthreads=nthreads, mask = mask)
      im = rdx_img_transform(invert(align_map), im, /preserve) >0
      if max(im) eq min(im) then stop
      red_missing, im, missing_type_wanted = 'median', /inplace
      cube[*, *, ifile] = im
                                ;cube[*, *, ifile] = red_rotate(im, self.direction)
    endfor                      ; ifile
    bg = median(cube)

    
    for ifile = 0, Nfiles-1 do begin
      red_progressbar, ifile, Nfiles, 'Stretching the images.'
      iexposure = indx[ifile]
      xshift = shift[0,iexposure]
      yshift = shift[1,iexposure]
      grid = reform(grids[iexposure,*,*,*])

;      print, wfiles[iexposure]
;      print, files[ifile]
;      tighttv, wcube[*, *, iexposure]    
;      tvscl, cube[*, *, ifile]
;
;
;      stop
      cube1[*,*,ifile] = red_rotation(cube[*,*,ifile] $
                                      , 0.0, xshift, yshift $
                                      , stretch_grid = grid $
                                      , nthreads=nthreads, nearest=nearest)
    endfor                      ; iexposure
    im = mean(cube1, dim = 3)
    red_missing, im, missing_type_wanted = 'nan', /inplace

    header = headfits(cfgdir+'/'+output_file + '.fitsheader')

    red_headerinfo_deletestep, header, /last
    self -> headerinfo_addstep, header $
                                , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                                , prpara = prpara $
                                , prproc = inam

    output_file = cfgdir+'/'+red_strreplace(output_file, 'results', 'results'+outdir_tag) + '.fits'
    fxaddpar, header, 'FILENAME', file_basename(output_file), 'Destretched raw data'
    fxaddpar, header, 'CROP_X0', x0, 'Crop llx pixel'  
    fxaddpar, header, 'CROP_X1', x1, 'Crop urx pixel'  
    fxaddpar, header, 'CROP_Y0', y0, 'Crop lly pixel'  
    fxaddpar, header, 'CROP_Y1', y1, 'Crop ury pixel'  
    red_writedata, output_file, im[x0:x1,y0:y1], header = header, overwrite = overwrite


    
    if trace then begin
    
;      tcube = fltarr([dims, Nfiles])  
;      tcube1 = fltarr([dims, Nfiles])
;      for ifile = 0, Nfiles-1 do begin
;        red_progressbar, ifile, Nfiles, 'Reading raw trace images for object'+strtrim(iobject, 2)
;        file = files[ifile]     ; <----------------------------------------------------
;        im = red_readdata(file, direction = self.direction)
;        im -= anchor_dark
;        im *= anchor_gain
;        mask = anchor_gain eq 0
;        tcube[*, *, ifile] = rdx_fillpix(im, nthreads=nthreads, mask = mask)
;      endfor                    ; ifile
      tcube1 = wcube1[*, *, indx]
      tim = mean(tcube1, dim = 3)
      red_missing, tim, missing_type_wanted = 'nan', /inplace
      theader = headfits(cfgdir+'/'+toutput_file + '.fitsheader')
      red_headerinfo_deletestep, theader, /last
      self -> headerinfo_addstep, theader $
                                  , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                                  , prpara = prpara $
                                  , prproc = inam
      toutput_file = cfgdir+'/'+red_strreplace(toutput_file, 'results', 'results'+outdir_tag) + '.fits'
      fxaddpar, theader, 'FILENAME', file_basename(toutput_file), 'Destretched raw data'
      fxaddpar, theader, 'CROP_X0', x0, 'Crop llx pixel'  
      fxaddpar, theader, 'CROP_X1', x1, 'Crop urx pixel'  
      fxaddpar, theader, 'CROP_Y0', y0, 'Crop lly pixel'  
      fxaddpar, theader, 'CROP_Y1', y1, 'Crop ury pixel'  
      red_writedata, toutput_file, tim[x0:x1,y0:y1], header =theader, overwrite = overwrite
    endif
    
  endfor                        ; iobject

  
end 

cd, '/scratch/mats/2016.09.19/CRISP-aftersummer'

a = crispred(/dev, /no)

cfgfile = 'momfbd_nopd/09:30:20/6302/cfg/momfbd_reduc_6302_00001.cfg'
cfgfile = 'momfbd_nopd/09:30:20/6302/cfg/momfbd_reduc_6302_00002.cfg'

cfgfile = 'momfbd_nopd/09:30:20/8542/cfg/momfbd_reduc_8542_00003.cfg'

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

