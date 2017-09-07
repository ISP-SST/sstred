; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
; 
;    xbd  : 
;   
;   
;   
;    ybd  : 
;   
;   
;   
;    np  : 
;   
;   
;   
;    clip  : 
;   
;   
;   
;    tile  : 
;   
;   
;   
;    tstep  : 
;   
;   
;   
;    scale  : 
;   
;   
;   
;    ang  : 
;   
;   
;   
;    shift  : 
;   
;   
;   
;    square : 
;   
;   
;   
;    fitsoutput : in, optional, type=boolean
;   
;      Set this to get output in a FITS file.
;   
;   
;   
;    negang  : 
;   
;   
;   
;    crop : 
;   
;    momfbddir :  in, optional, type=string, default='momfbd'
;   
;       Top directory of output tree.
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-07-24 : MGL. Limited tstep to length of scan.
; 
;   2013-08-27 : MGL. Added support for logging. Let the subprogram
;                find out its own name.
;
;   2013-09-11 : MGL. Use red_lp_write rather than lp_write.
;
;   2014-01-14 : PS. Code cleanup. Use self.filetype.
;
;   2014-01-15 : PS. Proper FITS header parsing. Support EXT_TIME for
;                all formats
;
;   2014-11-29 : JdlCR. added support for fullframe cubes (aka,
;                despite rotation and shifts, the entire FOV is inside
;                the image.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding SolarNet
;                keywords.
;
;   2016-10-28 : MGL. Adapt to new pipeline. 
;
;   2016-11-01 : MGL. New output file names including date and scan
;                selection. Changed default clip and tile. Improved
;                WB lightcurve plot. 
;
;   2016-11-03 : MGL. Removed keywords ext_time, ext_date, ang, and
;                shift. 
;
;   2016-11-30 : MGL. New keyword fitsoutput.
;
;   2017-01-30 : JdlCR. Allow to specify a mean offset angle to, e.g.,
;                create the cube aligned with the solar north. The
;                angle is in radians.
;
;   2017-03-20 : MGL. New keyword momfbddir.
;
;   2017-03-24 : MGL. Add prpara info.
;
;   2017-04-27 : MGL. New defaults for clip and tile.
;
;   2017-06-05 : MGL. Works with /fitsoutput now. Add 'BUNIT' and
;                'BTYPE' to FITS header.
;
;   2017-09-07 : MGL. Changed red_fitsaddpar --> red_addfitskeyword. 
;
;-
pro chromis::polish_tseries, blur = blur $
                             , clip = clip $
                             , crop = crop $
                             , fitsoutput = fitsoutput $
                             , fullframe = fullframe $
                             , momfbddir = momfbddir $
                             , negang = negang $
                             , np = np $
                             , offset_angle = offset_angle $
                             , scale = scale $
                             , square = square $
                             , tile = tile $
                             , timefiles = timefiles $
                             , tstep = tstep $
                             , xbd = xbd $
                             , ybd = ybd 
;                             , ang = ang $
;                             , ext_date = ext_date $
;                             , ext_time = ext_time $
;                             , shift = shift $

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if keyword_set(fitsoutput) then begin
    ;; Allow only in developer mode for now
    if ~self.developer_mode then begin
      print, inam + ' : fitsoutput allowed only in developer mode'
      retall
    endif
  endif
  
  ;; Make prpara
  if n_elements(blur        ) ne 0 then red_make_prpara, prpara, 'blur'         , blur         
  if n_elements(clip        ) ne 0 then red_make_prpara, prpara, 'clip'         , clip         
  if n_elements(crop        ) ne 0 then red_make_prpara, prpara, 'crop'         , crop         
  if n_elements(fitsoutput  ) ne 0 then red_make_prpara, prpara, 'fitsoutput'   , fitsoutput   
  if n_elements(fullframe   ) ne 0 then red_make_prpara, prpara, 'fullframe'    , fullframe    
  if n_elements(momfbddir   ) ne 0 then red_make_prpara, prpara, 'momfbddir'    , momfbddir    
  if n_elements(negang      ) ne 0 then red_make_prpara, prpara, 'negang'       , negang       
  if n_elements(np          ) ne 0 then red_make_prpara, prpara, 'np'           , np           
  if n_elements(offset_angle) ne 0 then red_make_prpara, prpara, 'offset_angle' , offset_angle 
  if n_elements(scale       ) ne 0 then red_make_prpara, prpara, 'scale'        , scale        
  if n_elements(square      ) ne 0 then red_make_prpara, prpara, 'square'       , square       
  if n_elements(tile        ) ne 0 then red_make_prpara, prpara, 'tile'         , tile         
  if n_elements(timefiles   ) ne 0 then red_make_prpara, prpara, 'timefiles'    , timefiles    
  if n_elements(tstep       ) ne 0 then red_make_prpara, prpara, 'tstep'        , tstep        
  if n_elements(ybd         ) ne 0 then red_make_prpara, prpara, 'ybd'          , ybd          
  if n_elements(xbd         ) ne 0 then red_make_prpara, prpara, 'xbd'          , xbd          

  ;; Default keywords
  if n_elements(momfbddir) eq 0 then momfbddir = 'momfbd' 
  
  ;; Camera/detector identification
  self->getdetectors
  wbindx = where(strmatch(*self.cameras,'Chromis-W'))
  wbcamera = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx = where(strmatch(*self.cameras,'Chromis-N')) 
  nbcamera = (*self.cameras)[nbindx[0]]
  nbdetector = (*self.detectors)[nbindx[0]]
  ;; Should be generalized to multiple NB cameras if CHROMIS gets
  ;; polarimetry. We don't need to identify any PD cameras for
  ;; restored data.

  ;; Find timestamp subdirs
  search_dir = self.out_dir +'/'+momfbddir+'/'
  timestamps = file_basename(file_search(search_dir + '*' $
                                         , count = Ntimestamps, /test_dir))
  if Ntimestamps eq 0 then begin
    print, inam + ' : No timestamp sub-directories found in: ' + search_dir
    return
  endif

  ;; Select timestamp folders
  selectionlist = strtrim(indgen(Ntimestamps), 2)+ '  -> ' + timestamps
  tmp = red_select_subset(selectionlist $
                          , qstring = inam + ' : Select timestamp directory ID:' $
                          , count = Ntimestamps, indx = sindx)
  if Ntimestamps eq 0 then begin
    print, inam + ' : No timestamp sub-folders selected.'
    return                      ; Nothing more to do
  endif
  timestamps = timestamps[sindx]
  print, inam + ' : Selected -> '+ strjoin(timestamps, ', ')

  ;; Loop over timestamp directories
  for itimestamp = 0L, Ntimestamps-1 do begin

    timestamp = timestamps[itimestamp]
    datestamp = self.isodate+'T'+timestamp

    ;; Find prefilter subdirs
    search_dir = self.out_dir +'/'+momfbddir+'/'+timestamp+'/'
    prefilters = file_basename(file_search(search_dir + '*' $
                                           , count = Nprefs, /test_dir))
    if Nprefs eq 0 then begin
      print, inam + ' : No prefilter sub-directories found in: ' + search_dir
      continue                  ; Next timestamp
    endif
    
    ;; Select prefilter folders
    selectionlist = strtrim(indgen(Nprefs), 2)+ '  -> ' + prefilters
    tmp = red_select_subset(selectionlist $
                            , qstring = inam + ' : Select prefilter directory ID:' $
                            , count = Nprefs, indx = sindx)
    if Nprefs eq 0 then begin
      print, inam + ' : No prefilter sub-folders selected.'
      continue                  ; Go to next timestamp
    endif
    prefilters = prefilters[sindx]
    print, inam + ' : Selected -> '+ strjoin(prefilters, ', ')

    ;; Loop over WB prefilters
    for ipref = 0L, Nprefs-1 do begin

      search_dir = self.out_dir + '/'+momfbddir+'/' + timestamp $
                   + '/' + prefilters[ipref] + '/cfg/results/'
      case self.filetype of
        'ANA': extension = '.f0'
        'MOMFBD': extension = '.momfbd'
        'FITS': extension = '.fits'
      endcase
      files = file_search(search_dir + '*'+self->getdetector('Chromis-W')+'*'+extension, count = Nfiles)

      ;; Find the global WB images and the number of scans.
;      self -> selectfiles, files = files, states = states $
;                           , cam = wbcamera, ustat = '' $
;                           , sel = windx, count = Nscans
      ;; We have no special state (or absence of state) to identify
      ;; the global WB images but we do know that their exposure times
      ;; are much larger than the ones corresponding to the individual
      ;; NB states.
      self -> extractstates, files, states
      windx = where(states.EXPOSURE gt mean(states.EXPOSURE))
      wstates = states[windx]
      wfiles = files[windx]

      ;; Select scan numbers
      selectionlist = wstates[uniq(wstates.scannumber, sort(wstates.scannumber))].scannumber
      tmp = red_select_subset(strtrim(selectionlist, 2) $
                              , qstring = inam + ' : Select scans:' $
                              , count = Nscans, indx = scanindx)

      uscans = selectionlist[scanindx]
      wstates = wstates[scanindx]
      wfiles = wfiles[scanindx]

      time = strarr(Nscans)
      date = strarr(Nscans)

      ;; Read headers to get obs_time and load the images into a cube
      for iscan = 0L, Nscans -1 do begin
        
        red_progressbar, iscan, Nscans, 'Read headers and load the images into a cube', clock = clock

        tmp = red_readdata(wfiles[iscan], head = hdr)

        if keyword_set(timefiles) then begin
          
          timefile = red_strreplace(wfiles[iscan], '.f0', '.time')
          spawn, 'cat '+timefile, ts
          ts_split = strsplit(ts, 'T', /extract)
          date[iscan] = ts_split[0]
          time[iscan] = ts_split[1]

        endif else begin

          ;; This part needs work! Only used the /timefiles option so
          ;; far. 

          date_avg = fxpar(hdr, 'DATE-AVG', count = hasdateavg)
          if hasdateavg then begin
            date_avg_split = strsplit(date_avg, 'T', /extract, count = Nsplit)
            ddate = date_avg_split[0]
            if Nsplit gt 1 then ttime = date_avg_split[1] else undefine, ttime
          endif else undefine, ddate, ttime

          if n_elements(ddate) eq 0 then begin
            print, inam+' : No date and time information for scan '+strtrim(uscans[iscan], 2)
            stop
          endif else begin
            date[iscan] = ddate
            time[iscan] = ttime
          endelse

        endelse

        IF n_elements(crop) NE 4 THEN crop = [0,0,0,0]
        
        IF iscan EQ 0 THEN BEGIN
          dim = size(tmp, /dimension)
          dimim = red_getborder(tmp, x0, x1, y0, y1, square = square)
          x0 += crop[0]
          x1 -= crop[1]
          y0 += crop[2]
          y1 -= crop[3]
          nx = x1 - x0 + 1
          ny = y1 - y0 + 1
          cub = fltarr(nx, ny, Nscans)
        endif
        
        ;; Why fillpix? Isn't this done in the momfbding? /MGL
        cub[*, *, iscan] = red_fillpix((temporary(tmp))[x0:x1, y0:y1], nthreads = 4L)
        
      endfor                    ; iscan

      if (keyword_set(fullframe)) then cub1 = cub

      ;; Get derotation angles
;      if(~keyword_set(ang)) then begin
        ang = red_lp_angles(time, date)
        mang = median(ang)
        ang -= mang
        if(keyword_set(negang)) then ang = -ang
        if(n_elements(offset_angle)) then ang += offset_angle
;      endif else begin
;
;        print, inam + ' : Using external angles'
;
;        if(n_elements(ang) NE Nscans) then begin
;          print, inam + ' : Error, the number of angles (' + red_stri(n_elements(ang)) $
;                 + ')!= number of images (' + red_stri(Nscans) + ')'
;          stop
;        endif
;      endelse
      
      ;; De-rotate images in the cube
      for iscan = 0L, Nscans -1 do begin
        red_progressbar, iscan, Nscans, inam+' : De-rotating images.', clock = clock
        cub[*,*,iscan] = red_rotation(cub[*,*,iscan], ang[iscan])
      endfor                    ; iscan
      
      ;; Align cube
;      if(~keyword_set(shift)) then begin
        if(~keyword_set(np)) then begin
          np = 0L
          read, np, prompt = inam +' : Please introduce the factor to recompute the reference image: '
        endif

        print, inam + ' : aligning images ... ', format = '(A, $)'
        shift = red_aligncube(cub, np, xbd = xbd, ybd = ybd, cubic = cubic, /aligncube)
        print, 'done'
;      endif else begin
;        print, inam + ' : Using external shifts'
;
;        if(n_elements(shift[0,*]) NE Nscans) then begin
;          print, inam + ' : Error, incorrect number of elements in shift array'
;          return
;        endif 
;
;        for iscan = 0L, Nscans - 1 do begin
;          red_progressbar, iscan, Nscans, inam+' : Applying shifts to images.', clock = clock
;          cub[*,*,iscan] = red_shift_im(cub[*,*,iscan], shift[0,iscan], shift[1,iscan])
;        endfor                  ; iscan
;
;      endelse


      if(keyword_set(fullframe)) then begin

        ;; Get maximum angle and maximum shift in each direction
        maxangle = max(abs(ang))
        mdx0 = reform(min(shift[0,*]))
        mdx1 = reform(max(shift[0,*]))
        mdy0 = reform(min(shift[1,*]))
        mdy1 = reform(max(shift[1,*]))
        ff = [maxangle, mdx0, mdx1, mdy0, mdy1]

        ;; Recreate cube
        dum = red_rotation(cub1[*,*,0], ang[0], shift[0,0], shift[1,0], full=ff)
        nd = size(dum,/dim)
        nx = nd[0]
        ny = nd[1]
        cub = fltarr([nd, Nscans])
        cub[*,*,0] = temporary(dum)

        for iscan=1, Nscans-1 do begin
          red_progressbar, clock = clock, iscan, Nscans $
                           , inam+' : Making full-size cube, de-rotating and shifting.'
          cub[*,*,iscan] = red_rotation(cub1[*,*,iscan], ang[iscan], shift[0,iscan] $
                                        , shift[1,iscan], full=ff)
        endfor                   ; iscan
        
        ;; Note that the fullsize option changes the FOV of the array.
        ;; This may make it necessary to adjust some header info.
 
      endif else ff = 0
      
      ;; De-stretch
      if(~keyword_set(clip)) then clip = [12, 6, 3, 1]
      if(~keyword_set(tile)) then tile = [10, 20, 30, 40]
      if(~keyword_set(scale)) then scale = 1.0 / float(self.image_scale)
      if(~keyword_set(tstep)) then begin
        dts = dblarr(Nscans)
        for iscan = 0L, Nscans - 1 do dts[iscan] = red_time2double(time[iscan])
        tstep = fix(round(180. / median(abs(dts[0:Nscans-2] - dts[1:*])))) <Nscans-1
      endif

      print, inam + ' : Using the following parameters for de-stretching the time-series: '
      print, '   tstep [~3 m. (?)]= ', tstep
      print, '   scale [pixels / arcsec] = ', scale
      print, '   tile = ['+strjoin(string(tile, format='(I3)'),',')+']'
      print, '   clip = ['+strjoin(string(clip, format='(I3)'),',')+']'

      grid = red_destretch_tseries(cub, scale, tile, clip, tstep)

      for iscan = 0L, Nscans - 1 do begin
        red_progressbar, iscan, Nscans, inam+' : Applying the stretches.', clock = clock
        cub[*,*,iscan] = red_stretch(cub[*,*,iscan], reform(grid[iscan,*,*,*]))
      endfor                    ; iscan

      ;; Measure time-dependent intensity variation (sun move's in the Sky)
      tmean = fltarr(Nscans)
      for ik = 0, Nscans-1 do tmean[ik] = median(cub[*,*,ik])
      
      cgplot, uscans, tmean, xtitle = 'Scan number', ytitle = 'Mean WB intensity', psym=-1

      ;; Prepare for making output file names
      midpart = prefilters[ipref] + '_' + datestamp + '_scans=' $
                + red_collapserange(uscans, ld = '', rd = '')
      if keyword_set(blur) then midpart += '_blurred'

      ;; Save angles, shifts and de-stretch grids
      odir = self.out_dir + '/calib_tseries/'
      file_mkdir, odir
      ofil = 'tseries_'+midpart+'_calib.sav'
      print, inam + ' : saving calibration data -> ' + odir + ofil
      save, file = odir + ofil $
            , tstep, clip, tile, scale, ang, shift, grid, time, date $
            , wfiles, tmean, crop, mang, x0, x1, y0, y1, ff, nd

      ;; Sample saved quantities with comments about their use in make_crispex:
      ;;
      ;; TSTEP           LONG      =            5          [Not used]
      ;; CLIP            INT       = Array[4]              [Not used]
      ;; TILE            INT       = Array[4]              [Not used]
      ;; SCALE           FLOAT     =       26.3852         [Not used]
      ;; ANG             DOUBLE    = Array[5]              [Yes]
      ;; SHIFT           DOUBLE    = Array[2, 5]           [Yes]
      ;; GRID            FLOAT     = Array[5, 2, 40, 25]   [Yes]
      ;; TIME            STRING    = Array[5]              [Not used]
      ;; DATE            STRING    = Array[5]              [Not used]
      ;; WFILES          STRING    = Array[5]              [Yes]
      ;; TMEAN           FLOAT     = Array[5]              [Yes]
      ;; CROP            INT       = Array[4]              [Not used]
      ;; MANG            DOUBLE    =        6.2653633      [Not used]
      ;; X0              LONG      =            0          [Yes]
      ;; X1              LONG      =         1831          [Yes]
      ;; Y0              LONG      =            0          [Yes]
      ;; Y1              LONG      =         1151          [Yes]
      ;; FF              INT       =        0              [Yes] ( = [maxangle, mdx0, mdx1, mdy0, mdy1] if fullframe)
      ;; ND              UNDEFINED = <Undefined>           [Yes] (2d array if fullframe)
      ;;
      ;; (Some of) these quantities should really go into the metadata
      ;; of the FITS file for two reasons: documentation and use in
      ;; make_crispex. The latter so we can avoid using the save file.
      ;;
      ;; Already in prpara: clip, crop, fullframe (same as ff?), grid.
      ;; So we could parse prpara to get those. Can make_crispex
      ;; calculate x0,x1,y0,y1 from naxis1 and naxis2 and any of the
      ;; other parameters?

      
      ;; Normalize intensity
      me = mean(tmean)
      for iscan = 0L, Nscans - 1 do cub[*,*,iscan] *= (me / tmean[iscan])

      if keyword_set(blur) then cub = smooth(cub, [29, 29, 1], /edge_wrap)

      if keyword_set(fitsoutput) then begin
        ;; Save WB results as a fits file
        ofil = 'wb_'+midpart+'_corrected_im.fits'
        print, inam + ' : saving WB corrected cube -> ' + odir + ofil

        ;; Add the wavelength and Stokes dimensions
        dims = [nx, ny, 1, 1, Nscans]
        cub = reform(cub, dims, /overwrite) 

        ;; Make header. Start with header from last input file, it has
        ;; info about MOMFBD processing.
        check_fits, cub, hdr, /update                          ; Get dimensions right
        red_fitsaddkeyword, hdr, 'DATE', red_timestamp(/iso) $     ; DATE with time
                        , 'Creation UTC date of FITS header'   ;
        red_fitsaddkeyword, hdr, 'BITPIX', 16 $                    ; Because we round() before saving. 
                        , 'Number of bits per data pixel'      ;
        red_fitsaddkeyword, hdr, 'FILENAME', ofil, anchor = 'DATE' ; New file name

        if keyword_set(blur) then begin
          red_fitsaddkeyword, hdr, before='DATE', 'COMMENT', 'Intentionally blurred version'
        endif

        ;; Add info to headers
        red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', 'DU', 'Units in array'
        red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

        ;; Add info about this step
        self -> headerinfo_addstep, hdr $
                                    , prstep = 'Prepare WB science data cube' $
                                    , prpara = prpara $
                                    , prproc = inam
        

        ;; Make time tabhdu extension with Nscans rows
        s_array = lonarr(Nscans)
        s_array[0] = wstates.scannumber
        t_array = dblarr(1, Nscans)
        t_array[0] = red_time2double(time)                       ; In [s] since midnight
;        w_array = fltarr(1)
        ;;w_array[0] = wstates.tun_wavelength
        w_array = replicate(float(prefilters[0])/10., 1, Nscans) ; In [nm]
;        tabhdu = {EXTNAME-WCS-TABLES: {TIME-TABULATION: {val:t_array  $
;                                                         , comment:'time-coordinates'}, $
;                                       WAVE-TABULATION: {val:w_array $
;                                                         , comment:'wavelength-coordinates'}, $
;                                       comment: ' For storing tabulated keywords'}}
;        tabhdu = {tabulations: {time: {val:t_array $
;                                       , comment:'time-coordinates'}, $
;                                wavelength: {val:w_array $
;                                             , comment:'wavelength-coordinates'}, $
;                                scannumber: {val:s_array $
;                                             , comment:'scannumbers'}, $
;                                comment: ' For storing tabulated keywords'}}
        
        ;; TIME reference value, all times are seconds since midnight.
        dateref = fxpar(hdr, 'DATEREF', count = count)
        if count eq 0 then begin
          ;; One should really check also for the existence of MJDREF
          ;; and JDREF but within the pipeline we can be sure we don't
          ;; use them.
          dateref = self.isodate+'T00:00:00.000000' ; Midnight
          red_fitsaddkeyword, hdr, 'DATEREF', dateref, 'Reference time in ISO-8601', after = 'DATE'
        endif

        help, round(cub)
        print, 'n_elements:', n_elements(cub)
        ;; Write the file
;        red_fits_createfile, odir + ofil, hdr, lun, fileassoc;, tabhdu = tabhdu
        self -> fitscube_initialize, odir + ofil, hdr, lun, fileassoc, dims ; $
;                                     , wcs_time_coordinate = t_array $
;                                     , wcs_wave_coordinate = w_array $
;                                     , scannumber = s_array


        for iscan = 0, Nscans-1 do begin
          self -> fitscube_addframe, fileassoc, round(cub[*, *, 0, 0, iscan]) $
                                     , iscan = iscan
        endfor                  ; iscan
        free_lun, lun
        print, inam + ' : Wrote file '+odir + ofil

;        ;; Experimental extra tabulated data:
;        Ntuning = 1
;        temp_array = fltarr(Ntemp)
;        r0_array   = fltarr(Nr0)
;        temp_array[0] = cos(s_array/max(s_array)) ; Really function of scannumber
;        r0_array[0]   = sin(s_array/max(s_array)) ; Really function of time, i.e., of tuning and scannumber
;        
;
;        red_fitsaddkeyword, hdr, 'TABULATD', 'TABULATIONS;ATMOS_R0,AMB_TEMP'
;
;        fxbhmake,bdr,1,'TABULATIONS','For storing tabulated keywords'
;        fxbaddcol, 1, bdr, r0_array, 'ATMOS_R0', TUNIT = 'm', 'Table of atmospheric r0'
;        fxbadd, bdr, '1CTYP1', 'TIME-TAB', 'Time since DATEREF, tabulated for ATMOS_R0'
;        fxbadd, bdr, '1CNAM1', 'Time since DATEREF for ATMOS_R0'
;        fxbadd, bdr, '1S1_0', 'TABULATIONS', 'Extension w/tabulations for ATMOS_R0'
;        fxbadd, bdr, '1S1_1', 'TIME-ATMOS_R0', 'TTYPE of col. w/TIME for ATMOS_R0'
;        fxbaddcol, 2, bdr, r0_time, 'TIME-ATMOS_R0', 'Tabulations of TIME for ATMOS_R0', tunit = 's'
;        
;        
;        fxbaddcol, 3, bdr, temp_array, 'AMB_TEMP', TUNIT = 'K', 'Table of ambient temperature'
;        fxbcreate, lun, filename, bdr, extension_no
;        fxbwrite, lun, r0_array, 1, 1
;        fxbwrite, lun, temp_array, 2, 1
;        fxbfinish, lun
;    

        self -> fitscube_addwcs, odir + ofil, w_array, t_array
;                                     , scannumber = s_array
;
;
;cccc = red_readdata(odir+ofil, head = hhhh)
;stop
;        
;fxbopen, tlun, odir + ofil, 'WCS-TAB', bhhh
;fxbread, tlun, wavetime_tab, 'WAVE+TIME'
;fxbclose, tlun
;
;
;
;stop
      endif else begin

        ;; Save WB results as lp_cube
        ofil = 'wb_'+midpart+'_corrected.icube'
        print, inam + ' : saving WB corrected cube -> ' + odir + ofil
        red_lp_write, fix(round(temporary(cub))), odir + ofil

      endelse


    endfor                      ; ipref
  endfor                        ; itimestamp


end
