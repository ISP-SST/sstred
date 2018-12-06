; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
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
; :history:
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
;   2014-01-14 : PS  Code cleanup.  Use self.filetype.  
;   2014-01-15 : PS  Proper FITS header parsing.  Support EXT_TIME for
;                all formats
;
;   2014-11-29 : JdlCR, added support for fullframe cubes (aka,
;                despite rotation and shifts, the entire FOV is inside
;                the image
;-
pro red::polish_tseries, xbd = xbd, ybd = ybd, np = np, clip = clip, $
                         tile = tile, tstep = tstep, scale = scale, $
                         ang = ang, shift = shift, square=square, $
                         negang = negang, crop=crop, ext_time = ext_time, $
                         fullframe = fullframe, ext_date = ext_date, offset_angle = offset_angle
  

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Logging
  help, /obj, self, output = selfinfo 
  red_writelog, selfinfo = selfinfo

  ;; Get time stamp
  data_dir = file_search(self.out_dir + '/momfbd/*', /test_directory, count = nf)
  if(nf GT 1) then begin
     print, inam + ' : Available folders: '
     for ii = 0L, nf - 1 do print, '  '+red_stri(ii)+' -> '+ data_dir[ii]
     idx = 0L
     read, idx, prompt='Select folder ID: '
     data_dir = data_dir[idx]
  endif 
  time_stamp = file_basename(data_dir)

  ;; Search reduction folders
  fold = file_search(data_dir+'/*', /test_directory, count = ct)
  if(ct eq 0) then begin
     print, inam + ' : Error, no subfolders were found in '+data_dir
     return
  endif

  ;; Select one of them if ct >1
  iread = 0
  if(ct gt 1) then begin
     print, inam + ' : reduction subfolders:'
     for ii = 0L, ct-1 do print, string(ii,format='(I2)') +' -> '+ file_basename(fold[ii])
     read, iread, prompt = 'Please, choose one state (enter the ID number): '
  endif
  fold = fold[iread]
  pref = file_basename(fold)
  print, inam + ' : selected state -> '+ pref
  
    ;; Extensions
  case self.filetype of
      'ANA': exten = '.f0'
      'MOMFBD': exten = '.momfbd'
      'FITS': exten =  '.fits'
      ELSE: begin
          print, inam+' : WARNING -> could not determine a file type for the output'
          exten = ''
      END
  endcase

  ;; Search files
  self->getcamtags, dir = self.data_dir
  wfiles = file_search(fold+'/cfg/results/'+self.camwbtag+'.?????.'+pref+exten, count = ct)
  if(ct eq 0) then begin
     print, inam + ' : Error, no WB files found in -> '+ fold+'/cfg/results/'
     stop
  endif

  time = strarr(ct)
  date = strarr(ct)

    ;; Read headers to get obs_time and load the images into a cube
  FOR ii = 0L, ct -1 DO BEGIN
      
      CASE self.filetype OF
          'ANA': BEGIN
              fzread, tmp, wfiles[ii], h
              dum = strsplit(h, ' =', /extract)
              IF n_elements(ext_time) GT 0 THEN $
                time[ii] = ext_time[ii] $
              ELSE $
                time[ii] = dum[1]
              date[ii] = dum[3]
          END
          'MOMFBD': BEGIN
              tmp = red_mozaic((dum =  momfbd_read(wfiles[ii])))
              IF n_elements(ext_time) GT 0 THEN $
                time[ii] = ext_time[ii] $
              ELSE $
                time[ii] = dum.time + ''
              date[ii] = strmid(dum.date, 0, 10) + ''
          END
          'FITS': BEGIN
              tmp = readfits(wfiles[ii], h, /SILENT)
              IF n_elements(ext_time) GT 0 THEN $
                time[ii] = ext_time[ii] $
              ELSE BEGIN
                  idx = where(strpos(h, 'TIME-OBS') GE 0)
                  ok = execute('time[ii] = '+(strsplit(strmid(h(idx), 9), '/', /extr))[0])
              ENDELSE
              idx = where(strpos(h, 'DATE-OBS') GE 0)
              ok = execute('date[ii] = '+(strsplit(strmid(h(idx), 9), '/', /extr))[0])
           END
      ENDCASE
      if(n_elements(ext_date) ne 0) then date[ii] = ext_date

      IF n_elements(crop) NE 4 THEN crop = [0,0,0,0]
      
      IF ii EQ 0 THEN BEGIN
          dim = size(tmp, /dimension)
          dimim = red_getborder(tmp, x0, x1, y0, y1, square = square)
          x0 += crop[0]
          x1 -= crop[1]
          y0 += crop[2]
          y1 -= crop[3]
          nx = x1 - x0 + 1
          ny = y1 - y0 + 1
          cub = fltarr(nx, ny, ct)
      ENDIF

      
      cub[*, *, ii] = red_fillpix((temporary(tmp))[x0:x1, y0:y1], nthreads = 4L)
  ENDFOR
  if (keyword_set(fullframe)) then cub1 = cub

  ;; Get derotation angles
  if(~keyword_set(ang)) then begin
     ang = red_lp_angles(time, date)
     mang = median(ang)
     ang -= mang
     if(n_elements(offset_angle)) then ang += offset_angle
     if(keyword_set(negang)) then ang = -ang
  endif else begin

     print, inam + ' : Using external angles'

     if(n_elements(ang) NE ct) then begin
        print, inam + ' : Error, the number of angles ('+red_stri(n_elements(ang))+')!= number of images ('+red_stri(ct)+')'
        stop
     endif
  endelse
  
  ;; De-rotate images in the cube
  print, inam+' : de-rotating WB images ... ', format = '(A,$)'
  for ii = 0L, ct -1 do cub[*,*,ii] = red_rotation(cub[*,*,ii], ang[ii])
  print, 'done'
  
  ;; Align cube
  if(~keyword_set(shift)) then begin
     if(~keyword_set(np)) then begin
        np = 0L
        read, np, prompt = inam +' : Please introduce the factor to recompute the reference image: '
     endif

     print, inam + ' : aligning images ... ', format = '(A, $)'
     shift = red_aligncube(cub, np, xbd = xbd, ybd = ybd, cubic = cubic, /aligncube)
     print, 'done'
  endif else begin
     print, inam + ' : Using external shifts'

     if(n_elements(shift[0,*]) NE ct) then begin
        print, inam + ' : Error, incorrect number of elements in shift array'
        return
     endif 

     for ii = 0L, ct - 1 do cub[*,*,ii] = red_shift_im(cub[*,*,ii], reform(shift[0,ii]), reform(shift[1,ii]))
;     FOR ii = 0L, ct-1 DO cub[*, *, i] = $
;       poly_2d(cub[*, *, ii], [-shift[0, ii], 0, 1., 0.], [-shift[1, ii], 1., 0., 0.], 2, CUBIC = -0.5)
  endelse


  if(keyword_set(fullframe)) then begin

     ; Get maximum angle and maximum shift in each direction
     maxangle = max(abs(ang))
     mdx0 = reform(min(shift[0,*]))
     mdx1 = reform(max(shift[0,*]))
     mdy0 = reform(min(shift[1,*]))
     mdy1 = reform(max(shift[1,*]))
     ff = [maxangle, mdx0, mdx1, mdy0, mdy1]

     ;; Recreate cube
     dum = red_rotation(cub1[*,*,0], ang[0], shift[0,0], shift[1,0], full=ff)
     nd = size(dum,/dim)
     cub = fltarr([nd, ct])
     cub[*,*,0] = temporary(dum)
     for ii=1, ct-1 do cub[*,*,ii] = red_rotation(cub1[*,*,ii], ang[ii], shift[0,ii], shift[1,ii], full=ff)
     
  endif else begin
     ff = 0
     ;for ii=0, ct-1 do cub[*,*,ii] = red_rotation(cub[*,*,ii], 0.0, sh[0,ii], sh[1,ii])
  endelse
  
  ;; De-stretch
  if(~keyword_set(clip)) then clip = [12,4,2,1]
  if(~keyword_set(tile)) then tile = [6,8,14,24]
  if(~keyword_set(scale)) then scale = 1.0 / float(self.image_scale)
  if(~keyword_set(tstep)) then begin
     dts = dblarr(ct)
     for ii = 0L, ct - 1 do dts[ii] = red_time2double(time[ii])
     tstep = fix(round(180. / median(abs(dts[0:ct-2] - dts[1:*])))) <ct
  endif

  print, inam + ' : Using the following parameters for de-stretching the time-series: '
  print, '   tstep [~3 m. (?)]= ', tstep
  print, '   scale [pixels / arcsec] = ', scale
  print, '   tile = ['+strjoin(string(tile, format='(I3)'),',')+']'
  print, '   clip = ['+strjoin(string(clip, format='(I3)'),',')+']'

  print, inam + " : computing destretch-grid (using LMSAL's routines) ... ", format = '(A, $)'
  grid = red_destretch_tseries(cub, scale, tile, clip, tstep)
  for ii = 0L, ct - 1 do cub[*,*,ii] = red_stretch(cub[*,*,ii], reform(grid[ii,*,*,*]))
  print, 'done'

  ;; Measure time-dependent intensity variation (sun move's in the Sky)
  tmean = total(total(cub,1),1) / float(nx) / float(ny)
  plot, tmean, xtitle = 'Time Step', ytitle = 'Mean WB intensity', psym=-1

  ;; Save angles, shifts and de-stretch grids
  odir = self.out_dir + '/calib_tseries/'
  file_mkdir, odir
  ofil = 'tseries.'+pref+'.'+time_stamp+'.calib.sav'
  print, inam + ' : saving calibration data -> ' + odir + ofil
  save, file = odir + ofil, tstep, clip, tile, scale, ang, shift, grid, time, date, wfiles, tmean, crop, mang, x0, x1, y0, y1, ff, nd

  ;; Normalize intensity
  me = mean(tmean)
  for ii = 0L, ct - 1 do cub[*,*,ii] *= (me / tmean[ii])

  ;; Save WB results as lp_cube
  ofil = 'wb.'+pref+'.'+time_stamp+'.corrected.icube'
  print, inam + ' : saving WB corrected cube -> ' + odir + ofil
                                ;fcwrite, fix(round(temporary(cub))), odir + ofil, ' '
  red_lp_write, fix(round(temporary(cub))), odir + ofil

  return
end
