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
; 
;-
pro red::polish_Tseries, xbd = xbd, ybd = ybd, np = np, clip = clip, $
                         tile = tile, tstep = tstep, scale = scale, $
                         ang = ang, shift = shift, square=square, negang = negang, $
                         crop=crop
  nostokes = 0B
  if(keyword_set(no_stokes)) then nostokes = 1B
  inam = 'red::polish_tseries : '
  ;;
  ;; Get time stamp
  ;;
  data_dir = file_search(self.out_dir + '/momfbd/*', count = nf)
  if(nf GT 1) then begin
     print, inam + 'Available folders: '
     for ii = 0L, nf - 1 do print, '  '+red_stri(ii)+' -> '+ data_dir[ii]
     idx = 0L
     read, idx, prompt='Select folder ID: '
     data_dir = data_dir[idx]
  endif 
  time_stamp = strsplit(data_dir, '/', /extract)
  time_stamp = time_stamp[n_elements(time_stamp)-1]
  

                                ;
                                ; Search reduction folders
                                ;
  fold = file_search(data_dir+'/*', /test_directory, count = ct)
  if(ct eq 0) then begin
     print, inam + 'Error, no subfolders were found in '+data_dir
     return
  endif
                                ;
                                ; select one of them if ct >1
                                ;
  iread = 0
  if(ct gt 1) then begin
     print, inam + 'reduction subfolders:'
     for ii = 0L, ct-1 do print, string(ii,format='(I2)') +' -> '+ file_basename(fold[ii])
     read, iread, prompt = 'Please, choose one state (enter the ID number): '
  endif
  fold = fold[iread]
  pref = file_basename(fold)
  print, inam + 'selected state -> '+ pref
                                ;
                                ; search files
                                ;
  self->getcamtags, dir = self.data_dir
  wfiles = file_search(fold+'/cfg/results/'+self.camwbtag+'.?????.'+pref+'.{momfbd,f0}', count = ct)
  if(ct eq 0) then begin
     print, inam + 'Error, no WB files found in -> '+ fold+'/cfg/results/'
     stop
  endif
                                ;
  time = strarr(ct)
  date = strarr(ct)
                                ;
                                ; Read headers to get obs_time and load the images into a cube
                                ;
  for ii = 0L, ct -1 do begin
     dum = strsplit(wfiles[ii],'.',/extract)
     ex =  dum[n_elements(dum)-1]
     print, ex
     if(ex eq 'f0') then begin 
        fzread, tmp, wfiles[ii], h
        if(n_elements(crop) ne 4) then crop = [0,0,0,0]

        if(ii eq 0) then begin
           dim = size(tmp, /dimension)
           dimim = red_getborder(tmp, x0, x1, y0, y1, square=square)
           x0 += crop[0]
           x1 -= crop[1]
           y0 += crop[2]
           y1 -= crop[3]
           nx = x1 - x0 + 1
           ny = y1 - y0 + 1
           cub = fltarr(nx, ny, ct)
        endif

        cub[*,*,ii] = red_fillzero((temporary(tmp))[x0:x1, y0:y1])
        dum = strsplit(h, ' =',/extract)
        time[ii] = dum[1]
        date[ii] = dum[3]
     endif 
     if(ex eq 'momfbd') then begin
        dum = momfbd_read(wfiles[ii])
        tmp = red_mozaic(dum)
        if(n_elements(crop) ne 4) then crop = [0,0,0,0] ; Added by MGL 2013-06-04

        if(ii eq 0) then begin
           dim = size(tmp, /dimension)
           dimim = red_getborder(tmp, x0, x1, y0, y1, square=square)
           x0 += crop[0]
           x1 -= crop[1]
           y0 += crop[2]
           y1 -= crop[3]
           nx = x1 - x0 + 1
           ny = y1 - y0 + 1
           cub = fltarr(nx, ny, ct)
        endif
        
        cub[*,*,ii] = red_fillzero((temporary(tmp))[x0:x1, y0:y1])
        time[ii] = dum.time + ''
        date[ii] = strmid(dum.date, 0, 10) + ''
     endif
  endfor
                                ;
                                ; get derotation angles
                                ;
  if(~keyword_set(ang)) then begin
     ang = red_lp_angles(time, date)
     ang -= median(ang)
     if(keyword_set(negang)) then ang = -ang
  endif else begin
                                ;
     print, inam + 'Using external angles'
                                ;
     if(n_elements(ang) NE ct) then begin
        print, inam + 'Error, the number of angles ('+red_stri(n_elements(ang))+')!= number of images ('+red_stri(ct)+')'
        stop
     endif
  endelse

                                ;
                                ; de-rotate images in the cube
                                ;
  print, inam+'de-rotating WB images ... ', format = '(A,$)'
  for ii = 0L, ct -1 do cub[*,*,ii] = red_rotation(cub[*,*,ii], ang[ii])
  print, 'done'
  
                                ;
                                ; align cube
                                ;
  if(~keyword_set(shift)) then begin
     if(~keyword_set(np)) then begin
        np = 0L
        read, np, prompt = inam +'Please introduce the factor to recompute the reference image: '
     endif
                                ;
     print, inam + 'aligning images ... ', format = '(A, $)'
     shift = red_aligncube(cub, np, xbd = xbd, ybd = ybd, cubic = cubic, /aligncube)
     print, 'done'
  endif else begin
     print, inam + 'Using external shifts'
                                ;
     if(n_elements(shift[0,*]) NE ct) then begin
        print, inam + 'Error, incorrect number of elements in shift array'
        return
     endif 
                                ;
     for ii = 0L, ct - 1 do cub[*,*,ii] = red_shift_im(cub[*,*,ii], reform(shift[0,ii]), reform(shift[1,ii]))
                                ;
  endelse
  
                                ;
                                ; De-stretch
                                ;
  if(~keyword_set(clip)) then clip = [12,4,2,1]
  if(~keyword_set(tile)) then tile = [6,8,8,14]
  if(~keyword_set(scale)) then scale = 1.0 / 0.0592
  if(~keyword_set(tstep)) then begin
     dts = dblarr(ct)
     for ii = 0L, ct - 1 do dts[ii] = red_time2double(time[ii])
     tstep = fix(round(180. / median(abs(dts[0:ct-2] - dts[1:*]))))
  endif
                                ;
  print, inam + 'Using the following parameters for de-stretching the time-series: '
  print, '   tstep [~3 m. (?)]= ', tstep
  print, '   scale [pixels / arcsec] = ', scale
  print, '   tile = ['+strjoin(string(tile, format='(I3)'),',')+']'
  print, '   clip = ['+strjoin(string(clip, format='(I3)'),',')+']'
                                ;
  print, inam + "computing destretch-grid (using T. Berger's routines) ... ", format = '(A, $)'
  grid = red_destretch_tseries(cub, scale, tile, clip, tstep)
  for ii = 0L, ct - 1 do cub[*,*,ii] = stretch(cub[*,*,ii], reform(grid[ii,*,*,*]))
  print, 'done'
  
                                ;
                                ; Measure time-dependent intensity variation (sun move's in the Sky)
                                ;
  tmean = total(total(cub,1),1) / float(nx) / float(ny)
  plot, tmean, xtitle = 'Time Step', ytitle = 'Mean WB intensity', psym=-1

                                ;
                                ; Save angles, shifts and de-stretch grids
                                ;
  odir = self.out_dir + '/calib_tseries/'
  file_mkdir, odir
  ofil = 'tseries.'+pref+'.'+time_stamp+'.calib.sav'
  print, inam + 'saving calibration data -> ' + odir + ofil
  save, file = odir + ofil, tstep, clip, tile, scale, ang, shift, grid, time, date, wfiles, tmean, crop
  

                                ;
                                ; Normalize intensity
                                ;
  me = mean(tmean)
  for ii = 0L, ct - 1 do cub[*,*,ii] *= (me / tmean[ii])


                                ;
                                ; Save WB results as lp_cube
                                ;
  ofil = 'wb.'+pref+'.'+time_stamp+'.corrected.icube'
  print, inam + 'saving WB corrected cube -> ' + odir + ofil
                                ;fcwrite, fix(round(temporary(cub))), odir + ofil, ' '
  lp_write, fix(round(temporary(cub))), odir + ofil
                                ;
  return
end
