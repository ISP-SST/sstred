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
;    negang  : 
;   
;   
;   
;    crop : 
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
;
;-
pro chromis::polish_tseries, xbd = xbd $
                             , ybd = ybd $
                             , np = np $
                             , clip = clip $
                             , tile = tile $
                             , tstep = tstep $
                             , scale = scale $
;                             , ang = ang $
;                             , shift = shift $
                             , square = square $
                             , negang = negang $
                             , crop = crop $
;                             , ext_time = ext_time $
                             , fullframe = fullframe $
;                             , ext_date = ext_date $
                             , timefiles = timefiles $
                             , fitsoutput = fitsoutput $
                             , blur = blur
  

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

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
  search_dir = self.out_dir +'/momfbd/'
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
    search_dir = self.out_dir +'/momfbd/'+timestamp+'/'
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

      search_dir = self.out_dir + '/momfbd/' + timestamp $
                   + '/' + prefilters[ipref] + '/cfg/results/'
      case self.filetype of
        'ANA': extension = '.f0'
        'MOMFBD': extension = '.momfbd'
        'FITS': extension = '.fits'
      endcase
      files = file_search(search_dir + '*'+extension, count = Nfiles)      
      
      ;; Find the global WB images and the number of scans
      self -> selectfiles, files = files, states = states $
                           , cam = wbcamera, ustat = '' $
                           , sel = windx, count = Nscans $
                           , complement = complement, Ncomplement = Ncomplement
      wstates = states[windx]
      wfiles = files[windx]

      ;; Select can numbers
      selectionlist = strtrim(wstates[uniq(wstates.scannumber, sort(wstates.scannumber))].scannumber, 2)
      tmp = red_select_subset(selectionlist $
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

        tmp = red_readdata(wfiles[iscan], head = h)

        if keyword_set(timefiles) then begin
          
          timefile = strreplace(wfiles[iscan], '.f0', '.time')
          spawn, 'cat '+timefile, ts
          ts_split = strsplit(ts, 'T', /extract)
          date[iscan] = ts_split[0]
          time[iscan] = ts_split[1]

        endif else begin

          ;; This part needs work! Only used the /timefiles option so
          ;; far. 

          date_ave = fxpar(h, 'DATE-AVE', count = hasdateave)
          if hasdateave then begin
            date_ave_split = strsplit(date_ave, 'T', /extract, count = Nsplit)
            ddate = date_ave_split[0]
            if Nsplit gt 1 then ttime = date_ave_split[1] else undefine, ttime
          endif else undefine, ddate, ttime

          if n_elements(ddate) eq 0 then begin
            print, inam+' : No date and time information for scan '+strtrim(uscans[iscan], 2)
            stop
          endif else begin
            date[iscan] = ddate
            time[iscan] = ttime
          endelse

;          if n_elements(ext_time) gt 0 then begin
;            time[iscan] = ext_time[iscan] 
;          endif else begin
;            if n_elements(ttime) gt 0 then time[iscan] = ttime
;          endelse
;          if(n_elements(ext_date) ne 0) then begin
;            date[iscan] = ext_date
;          endif else begin
;            date[iscan] = ddate
;          endelse
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
        
        cub[*, *, iscan] = red_fillpix((temporary(tmp))[x0:x1, y0:y1], nthreads = 4L)
        
      endfor                    ; iscan

      if (keyword_set(fullframe)) then cub1 = cub

      ;; Get derotation angles
;      if(~keyword_set(ang)) then begin
        ang = red_lp_angles(time, date)
        mang = median(ang)
        ang -= mang
        if(keyword_set(negang)) then ang = -ang
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
        cub = fltarr([nd, Nscans])
        cub[*,*,0] = temporary(dum)
        for iscan=1, Nscans-1 do begin
          red_progressbar, clock = clock, iscan, Nscans $
                           , inam+' : Making full-size cube, de-rotating and shifting.'
          cub[*,*,iscan] = red_rotation(cub1[*,*,iscan], ang[iscan], shift[0,iscan] $
                                        , shift[1,iscan], full=ff)
        endfor                   ; iscan
 
      endif else ff = 0
      
      ;; De-stretch
      if(~keyword_set(clip)) then clip = [8, 4, 1]
      if(~keyword_set(tile)) then tile = [8, 16, 32]
      if(~keyword_set(scale)) then scale = 1.0 / float(self.image_scale)
      if(~keyword_set(tstep)) then begin
        dts = dblarr(Nscans)
        for iscan = 0L, Nscans - 1 do dts[iscan] = red_time2double(time[iscan])
        tstep = fix(round(180. / median(abs(dts[0:Nscans-2] - dts[1:*])))) <Nscans
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
      tmean = total(total(cub,1),1) / float(nx) / float(ny)
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

      ;; Normalize intensity
      me = mean(tmean)
      for iscan = 0L, Nscans - 1 do cub[*,*,iscan] *= (me / tmean[iscan])

      if keyword_set(blur) then cub = smooth(cub, [29, 29, 1], /edge_wrap)

      if keyword_set(fitsoutput) then begin
        ;; Save WB results as a fits file
        ofil = 'wb_'+midpart+'_corrected_im.fits'
        print, inam + ' : saving WB corrected cube -> ' + odir + ofil

        ;; Add the wavelength dimension
        cub = reform(cub, nx, ny, 1, Nscans, /overwrite) 

        ;; Make header
        dims = size(cub, /dim)
        type=(size(fix(1)))[1]
        mkhdr,hdr,type,dims,/extend
        ;; Keywords to describe data cube
        sxaddpar, hdr, before='DATE', 'SOLARNET', 0.5
        sxaddpar, hdr, before='DATE', 'BTYPE',    'Intensity'
        sxaddpar, hdr, before='DATE', 'BUNIT',    'DU' ; Digital unit?
        sxaddpar, hdr, before='DATE', 'CDELT1',   float(self.image_scale), ' [arcsec] x-coordinate increment'
        sxaddpar, hdr, before='DATE', 'CDELT2',   float(self.image_scale), ' [arcsec] y-coordinate increment'
        sxaddpar, hdr, before='DATE', 'CDELT3',   1.,               ' [m] wavelength-coordinate axis increment'
        sxaddpar, hdr, before='DATE', 'CDELT4',   1.,               ' [s] time-coordinate axis increment'
        sxaddpar, hdr, before='DATE', 'CTYPE1',   'x',              ' [arcsec]'
        sxaddpar, hdr, before='DATE', 'CTYPE2',   'y',              ' [arcsec]'
        sxaddpar, hdr, before='DATE', 'CTYPE3',   'wavelength',     ' [m]'
        sxaddpar, hdr, before='DATE', 'CTYPE4',   'time',           ' [s]'
        sxaddpar, hdr, before='DATE', 'CUNIT1',   'arcsec'
        sxaddpar, hdr, before='DATE', 'CUNIT2',   'arcsec'
        sxaddpar, hdr, before='DATE', 'CUNIT3',   'm'
        sxaddpar, hdr, before='DATE', 'CUNIT4',   's'
        sxaddpar, hdr, before='DATE', 'COMMENT',  "Index order is (x,y,lambda,t)"
        ;; DATE_OBS should be getting the value including decimals
        ;; from the raw data headers, not just integer seconds as
        ;; here:
        sxaddpar, hdr, before='DATE', 'DATE_OBS', datestamp,        ' '
        ;; IS STARTOBS needed?
        sxaddpar, hdr, before='DATE', 'STARTOBS', datestamp,        ' '
        ;; Keywords from wstates
        sxaddpar, hdr, before='DATE', 'DETECTOR', wstates[0].detector
        sxaddpar, hdr, before='DATE', 'CAMERA',   wstates[0].camera
        sxaddpar, hdr, before='DATE', 'IS_WB',    wstates[0].is_wb
        if wstates[0].cam_settings ne '' then $
           sxaddpar, hdr, before='DATE', 'CAM_SETTINGS',  wstates[0].cam_settings
        if wstates[0].pf_wavelength ne 0.0 then $
           sxaddpar, hdr, before='DATE', 'PF_WAVELENGTH', wstates[0].pf_wavelength
        if wstates[0].exposure ne 0.0d then $
           sxaddpar, hdr, before='DATE', 'EXPOSURE',      wstates[0].exposure
        ;; Keyword info from other variables
        sxaddpar, hdr, before='DATE', 'FILTER1',  prefilters[ipref]

        if keyword_set(blur) then begin
          sxaddpar, hdr, before='DATE', 'COMMENT', 'Intentionally blurred version'
        endif

        ;; Make time tabhdu extension with Nscans rows
        s_array = lonarr(1, 1, 1, Nscans)
        s_array[0] = wstates.scannumber
        t_array = dblarr(1, 1, 1, Nscans)
        t_array[0] = red_time2double(time);-red_time2double(timestamp)
        w_array = fltarr(1, 1, 1, 1)
        ;;w_array[0] = wstates.tun_wavelength
        w_array[0] = float(prefilters[0])*1e-10
        tabhdu = {tabulations: {time: {val:t_array $
                                       , comment:'time-coordinates'}, $
                                wavelength: {val:w_array $
                                             , comment:'wavelength-coordinates'}, $
                                scannumber: {val:s_array $
                                             , comment:'scannumbers'}, $
                                comment: ' For storing tabulated keywords'}}
stop
        red_fits_createfile, odir + ofil, hdr, lun, fileassoc, tabhdu = tabhdu
        for iscan = 0, Nscans-1 do fileassoc[iscan] = round(cub[*, *, *, iscan])
        free_lun, lun
        print, inam + ' : Wrote file '+odir + ofil

;        ;; Write the file
;        red_writedata, odir + ofil, fix(round(temporary(cub))) $
;                       , header = hdr, tabhdu = tabstruct, /over
 

rcub = readfits(odir + ofil, rhdr)
rdcub = red_readdata(odir + ofil, head = rdhdr)
help, cub
print, hdr, format = '(a0)'
help, rcub
print, rhdr, format = '(a0)'
;help, rdcub
;print, rdhdr, format = '(a0)'
stats, round(cub)-rcub
stop
      endif else begin

        ;; Save WB results as lp_cube
        ofil = 'wb_'+midpart+'_corrected.icube'
        print, inam + ' : saving WB corrected cube -> ' + odir + ofil
        red_lp_write, fix(round(temporary(cub))), odir + ofil

      endelse


    endfor                      ; ipref
  endfor                        ; itimestamp


end
