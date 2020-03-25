; docformat = 'rst'

;+
; Check orientation and position of the FOV of a fitscube file by
; comparison with an SDO/HMI intensity image collected close in time.
; 
; Use the iframe, ituning, istokes, and/or iscan keywords to select
; the fitscube frame to compare.
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
;       The name of a WB or NB fitscube file.
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
;-
pro red::compare_hmi, filename $
                      , iframe = iframe $
                      , iscan = iscan $
                      , istokes = istokes $
                      , ituning = ituning $
                      , log = log


  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Download HMI images
  self -> download, /hmi

  ;; Read the SST data
  h = headfits(filename)
  red_fitscube_getframe, filename, im_sst $
                         , iframe = iframe $
                         , iscan = iscan $
                         , istokes = istokes $
                         , ituning = ituning
  dims_sst = size(im_sst, /dim)
  date_beg = fxpar(h, 'DATE-BEG')
  time_beg = red_time2double((strsplit(date_beg, 'T', /extract))[1])

  ;; SST FOV WCS info, array coordinates in arc sec
  red_fitscube_getwcs, filename, coordinates=wcs

  ;; SST FOV center coordinates in arc sec
  xc_sst = mean(wcs[iframe].hpln) 
  yc_sst = mean(wcs[iframe].hplt)

  red_show, im_sst, w = 0 
  
  ;; Find the HMI image closest in time
  hnames=file_search('downloads/HMI/*Ic_flat_4k.jpg')
  htimes = float(strmid(file_basename(hnames), 9, 2)) * 3600.
  mn = min(abs(time_beg - htimes),ml)
  hname = hnames[ml]
  read_jpeg, hname, im_hmi
  im_hmi = reform(im_hmi[1, *, *])          ; use green channel
  hmi_image_scale = 0.5                     ; arc sec / pixel
  dims_hmi = size(im_hmi, /dim)

 
  ;; Display the HMI image
  scrollwindow, xs = dims_hmi[0], ys = dims_hmi[1], wid = wid_hmi, /free, sizefraction = 0.75
  tvscl, im_hmi
  stop
  ;; Draw axes
  cgPolygon, [0, dims_hmi[0], dims_hmi[0]/2, dims_hmi[0]/2, dims_hmi[0]/2] $
             , [dims_hmi[1]/2, dims_hmi[1]/2, dims_hmi[1]/2, 0, dims_hmi[1]] $
             , color='yellow', /DEVICE
  ;; Draw FOV box
  cgPolygon, round(wcs[iframe].hpln[[0, 1, 3, 2, 0]]/hmi_image_scale+dims_hmi[0]/2) $
             , round(wcs[iframe].hplt[[0, 1, 3, 2, 0]]/hmi_image_scale+dims_hmi[1]/2) $
             , color='yellow', /DEVICE
  
  print, 'Please use scroll bars to center the HMI image on the SST FOV.'
  print, 'It should be close to the drawn rectangle.'
  print, 'We will worry about the alignment later but first:'
  s = ''
  read, 'Do the orientations (90 deg rotations and mirroring) match in the two images [yN]? ', s
  if strupcase(strmid(s, 0, 1)) ne 'Y' then begin

    ;; Help determining the direction parameter, ask user to make new
    ;; cube, then exit.

    ;; The direction is available in the class object but you can
    ;; override it when making the WB cube, so better get it from the
    ;; file header.
    pos = where(strmatch(fxpar(h,'PRPROC*'),'*make_wb_cube*'), Nmatch)
    if Nmatch eq 0 then stop
    prpara = json_parse((fxpar(h,'PRPARA*'))[pos])
    if prpara.haskey('DIRECTION') then begin
      ;; Direction in the header
      direction = prpara['DIRECTION']
    endif else begin
      ;; Direction from the class object, i.e., from the config file.
      direction = self.direction
    endelse
    
    ;; First invert the current "direction" of the SST data, then
    ;; determine the "direction" needed to get the matching
    ;; orientaiton.
    dims_sst_small = round(dims_sst*float(self.image_scale)/hmi_image_scale)
    im_sst_orig = congrid(red_rotate(im_sst, -direction) $
                          , dims_sst_small[0], dims_sst_small[1])
    maxdim_small = max(dims_sst_small)
    im_sst_orig = red_centerpic(im_sst_orig, sz = maxdim_small, z = median(im_sst_orig))
    ;;scrollwindow, xs = maxdim_small*8, ys = maxdim_small, wid = wid_sst_small, /free
    window, xs = maxdim_small*8, ys = maxdim_small
    for i = 0, 7 do begin
      tvscl, red_rotate(im_sst_orig, i), i
      cgtext, i*maxdim_small+10, 10, strtrim(i, 2), /dev, color = 'yellow'
    endfor                      ; i

    ;; Change this to click in window.
    repeat read, 'Which panel is oriented like the HMI image [0-7]? ', s $
    until long(s) ge 0 and long(s) le 7
    
    if long(s) eq direction then begin
      print, 'This is the orientation of the cube, no change.'
    endif else begin
      print, 'Please edit config.txt, change or add the following line:'
      print, 'direction='+s
      print, 'Then make a new WB cube and continue from that.'
      retall
    endelse
    
  endif 

  print, 'Will now try to measure the misalignment with cross correclation.'

  ;; SST FOV in the HMI image
  xc_hmi_pix = dims_hmi[0]/2. + xc_sst/hmi_image_scale
  yc_hmi_pix = dims_hmi[1]/2. + yc_sst/hmi_image_scale

  xsz_hmi_pix = round(dims_sst[0] * float(self.image_scale)/hmi_image_scale)
  ysz_hmi_pix = round(dims_sst[1] * float(self.image_scale)/hmi_image_scale)

  subim_hmi = congrid(red_pic_at_coord(im_hmi $
                                       , round(xc_hmi_pix), round(yc_hmi_pix) $
                                       , xsz_hmi_pix, ysz_hmi_pix) $
                      , dims_sst[0], dims_sst[1], cub=-0.5)
  red_show, subim_hmi, w = 1

  shifts = red_alignoffset(red_centerpic(im_sst,sz=1024) $
                           , red_centerpic(subim_hmi,sz=1024),/win)

  xc_hmi_pix -= shifts[0] * float(self.image_scale)/hmi_image_scale
  yc_hmi_pix -= shifts[1] * float(self.image_scale)/hmi_image_scale
  subim_hmi = congrid(red_pic_at_coord(im_hmi $
                                       , round(xc_hmi_pix), round(yc_hmi_pix) $
                                       , xsz_hmi_pix, ysz_hmi_pix) $
                      , dims_sst[0], dims_sst[1], cub=-0.5)
  red_show, subim_hmi, w = 1
  blink, [0, 1]

  read, 'Did the alignment work [yN]? ', s
  wdelete, 1

  if strupcase(strmid(s, 0, 1)) ne 'Y' then begin

    print, 'Please click on a recognicable structure in the SST image.'
    wshow, 0
    wset, 0
    cursor, xp_sst_pix, yp_sst_pix, /device

    print, 'Click on the same structure in the HMI image.'
    wshow, wid_hmi
    wset, wid_hmi
    cursor, xp_hmi_pix, yp_hmi_pix, /device

    ;; Now read out a subim from im_hmi that matches the sst image.
    xl = round(xp_hmi_pix - xp_sst_pix*float(self.image_scale)/hmi_image_scale)
    xh = round(xp_hmi_pix + (dims_sst[0]-xp_sst_pix)*float(self.image_scale)/hmi_image_scale)
    yl = round(yp_hmi_pix - yp_sst_pix*float(self.image_scale)/hmi_image_scale)
    yh = round(yp_hmi_pix + (dims_sst[1]-yp_sst_pix)*float(self.image_scale)/hmi_image_scale)
    xc_hmi_pix = (xl+xh)/2
    yc_hmi_pix = (yl+yh)/2
    subim_hmi = congrid(red_pic_at_coord(im_hmi $
                                         , round(xc_hmi_pix), round(yc_hmi_pix) $
                                         , xsz_hmi_pix, ysz_hmi_pix) $
                        , dims_sst[0], dims_sst[1], cub=-0.5)
    red_show, subim_hmi, w = 1

    shifts = red_alignoffset(red_centerpic(im_sst,sz=1024) $
                             , red_centerpic(subim_hmi,sz=1024),/win)

    xc_hmi_pix -= shifts[0] * float(self.image_scale)/hmi_image_scale
    yc_hmi_pix -= shifts[1] * float(self.image_scale)/hmi_image_scale
    subim_hmi = congrid(red_pic_at_coord(im_hmi $
                                         , round(xc_hmi_pix), round(yc_hmi_pix) $
                                         , xsz_hmi_pix, ysz_hmi_pix) $
                        , dims_sst[0], dims_sst[1], cub=-0.5)
    red_show, subim_hmi, w = 1
    blink, [0, 1]
    
    read, 'Did the alignment work this time [yN]? ', s
    if strupcase(strmid(s, 0, 1)) ne 'Y' then begin
      stop
    endif
   
  endif

  ;; Calculate total shifts in arc sec
  xc_hmi = (xc_hmi_pix - dims_hmi[0]/2.) * hmi_image_scale 
  yc_hmi = (yc_hmi_pix - dims_hmi[1]/2.) * hmi_image_scale
  x_shift = xc_hmi - xc_sst
  y_shift = yc_hmi - yc_sst
  
  print, 'Detected X displacement: ' + strtrim(x_shift, 2) + ' arc sec.'
  print, 'Detected Y displacement: ' + strtrim(y_shift, 2) + ' arc sec.'

  ;; Modify the WCS info (for all frames)
  wcs.hpln += x_shift
  wcs.hplt += y_shift

  ;; Confirm that the new coordinates are good.
  wshow, wid_hmi
  wset, wid_hmi
  ;; Draw new FOV box
  cgPolygon, round(wcs[iframe].hpln[[0, 1, 3, 2, 0]]/hmi_image_scale+dims_hmi[0]/2) $
             , round(wcs[iframe].hplt[[0, 1, 3, 2, 0]]/hmi_image_scale+dims_hmi[1]/2) $
             , color='green', /DEVICE
  
  ;; If we want to also offer to measure the rotation angle vs the HMI
  ;; image, this is where we could do it. The correction as a rotation
  ;; matrix could be specified in the PCm_n keywords of the fitscube
  ;; file. But for that we need to first update red_fitscube_getwcs so
  ;; that it can use that information. (We probably have to do that
  ;; anyway if we want to be able to get the spatial coordinates right
  ;; in un-rotated scan cubes and wb/nb cubes where an offset angle is
  ;; used, or where the average rotation angle is subtracted to save
  ;; space.)
  
  ;; Or we skip measuring (small) rotation misalignments in this step,
  ;; leaving that for software that align data from multiple
  ;; instruments.


  if keyword_set(log) then begin
    lfile = 'compare_hmi.log'
    if file_test(lfile) then openu, llun, lfile, /get_lun else openw, llun, lfile, /get_lun
    printf, llun, date_beg+' '+strtrim(x_shift, 2)+' '+strtrim(y_shift, 2)
    free_lun, llun
  endif
  
  read, 'Do you want to update the WCS info in the cube file [yN]? ', s

  if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
    
    ;; Write the updated WCS info back into the file. 
    red_fitscube_addwcs, filename, wcs, dimensions = fxpar(h, 'NAXIS*'), /update

  endif

end

